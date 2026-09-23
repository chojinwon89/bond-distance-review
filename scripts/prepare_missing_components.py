#!/usr/bin/env python3
"""Prepare isolated jobs for audited website gaps, retaining all source inputs."""
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib,json,re,shutil
from pathlib import Path
import numpy as np
from ase.constraints import FixAtoms
from ase.io import read,write
from add_missing_dft_images import incar_values
from component_store import reference_candidates
from extract_kestrel_singlepoint import infer_functional
from prepare_completion_jobs import sha

FUNCS={'pbe':'PBE','pbe_d3':'PBE_D3','r2scan':'r2scan','beef_vdw':'beef_vdw'}


def potential_blocks(path):
    blocks=re.findall(r'.*?End of Dataset[^\n]*\n?',path.read_text(),re.S)
    result={}
    for block in blocks:
        match=re.search(r'TITEL\s*=\s*([^\n]+)',block)
        element=re.search(r'VRHFIN\s*=\s*([A-Z][a-z]?)\s*:',block)
        if not match or not element:raise ValueError('Unrecognized potential block: '+str(path))
        result[match[1].strip()]=(element[1],block)
    if not result:raise ValueError('No POTCAR datasets: '+str(path))
    return result


def bottom_constraints(atoms,template,metal):
    """Transfer a verified lower-layer constraint convention, not atom indices."""
    left=[i for i,a in enumerate(template) if a.symbol==metal]
    right=[i for i,a in enumerate(atoms) if a.symbol==metal]
    if len(left)!=len(right) or not np.allclose(template.cell,atoms.cell,atol=1e-5,rtol=0):
        raise ValueError('Surface template cell/count mismatch')
    fixed=set()
    for c in template.constraints:
        if not isinstance(c,FixAtoms):raise ValueError('Unsupported constraint type')
        fixed.update(c.get_indices())
    normal=np.cross(atoms.cell[0],atoms.cell[1]);normal/=np.linalg.norm(normal)
    old_z=template.positions@normal;new_z=atoms.positions@normal
    left.sort(key=lambda i:old_z[i]);right.sort(key=lambda i:new_z[i])
    if not fixed or set(left[:len(fixed)])!=fixed or len(fixed)>=len(right):
        raise ValueError('Template does not constrain only lower metal layers')
    n=len(fixed)
    if new_z[right[n]]-new_z[right[n-1]]<0.5:raise ValueError('Constraint boundary splits a layer')
    atoms.set_constraint(FixAtoms(indices=right[:n]))
    return n


def settings_for(original,functional,role,robust=True):
    v=dict(original)
    if infer_functional(v)!=functional:raise ValueError('Functional template mismatch')
    v.update(ISTART='0',ICHARG='2',LWAVE='.FALSE.',LCHARG='.FALSE.')
    if robust:
        v.pop('IALGO',None);v.update(ALGO='Normal',NELM=str(max(300,int(v.get('NELM','150')))))
    if role=='spe':v.update(NSW='0',IBRION='-1');v.pop('EDIFFG',None)
    else:v.update(NSW='1000',IBRION='1',POTIM='0.2',EDIFFG=original.get('EDIFFG','-0.05'),EDIFF=str(min(float(original.get('EDIFF','1e-5')),1e-6)))
    return v


def validate_inputs(directory):
    atoms=read(directory/'POSCAR',format='vasp')
    if not np.isfinite(atoms.positions).all() or abs(np.linalg.det(atoms.cell))<1e-5:raise ValueError('Invalid coordinates')
    distances=atoms.get_all_distances(mic=True);np.fill_diagonal(distances,np.inf)
    if len(atoms)>1 and distances.min()<0.5:raise ValueError('Near-coincident atoms')
    labels=(directory/'POSCAR').read_text().splitlines()[5].split()
    blocks=potential_blocks(directory/'POTCAR')
    # Our generated structures have one block per species; native inputs must too.
    if [v[0] for v in blocks.values()]!=labels:raise ValueError('POSCAR/POTCAR element ordering mismatch')
    if not (directory/'KPOINTS').stat().st_size:raise ValueError('Empty KPOINTS')
    settings=incar_values((directory/'INCAR').read_text())
    if infer_functional(settings)=='beef_vdw' and not (directory/'vdw_kernel.bindat').is_file():raise ValueError('Missing vdW kernel')
    return atoms


def prepare(root,project,campaign,plans,records,live_snapshot):
    if campaign.exists() or campaign.parent!=project/'completion_runs':raise ValueError('New isolated campaign required')
    bydir={r['calculation']['directory']:r for r in records if r['cluster']=='perlmutter'}
    groups=defaultdict(list)
    for r in records:
        if r['role']=='complex':groups[(r['surface'],r['molecule'],r['calculation']['functional'])].append(r)
    local=[r for r in records if r['cluster']=='perlmutter' and all((Path(r['calculation']['directory'])/n).is_file() for n in ['INCAR','POSCAR','POTCAR','KPOINTS'])]
    donors=defaultdict(list)
    for r in local:
        for element,title in r['calculation'].get('potentials',{}).items():donors[title].append(Path(r['calculation']['directory'])/'POTCAR')
    cached_blocks={}
    def block(title):
        for path in donors[title]:
            if path not in cached_blocks:cached_blocks[path]=potential_blocks(path)
            if title in cached_blocks[path]:return cached_blocks[path][title][1],path
        raise ValueError('No local potential donor for '+title)
    geometry=json.loads((root/'functional_geometry.json').read_text())['systems']
    drafts=[];ready=[]
    for plan in plans:
        s,m,f,mode=plan['key'];metal=re.match(r'[A-Z][a-z]?',s)[0]
        if plan['status']!='new complex':ready.append(plan);continue
        role=mode;source=Path(plan['source']) if plan['source'] else None
        if source:
            settings=incar_values((source/'INCAR').read_text())
            if infer_functional(settings)!=f:raise ValueError('Need a same-functional source for '+str(plan['key']))
            geo=source/'POSCAR'
            # Failed electronic steps can corrupt the subsequent ionic move.
            # Reuse saved ionic geometry only when the source was not invalid/SCF-failed.
            src=bydir[str(source)]['calculation']
            if role=='relaxed' and src['status'] not in ('invalid','electronic_unconverged') and (source/'CONTCAR').exists():
                try:read(source/'CONTCAR',format='vasp');geo=source/'CONTCAR'
                except Exception:pass
            atoms=read(geo,format='vasp')
            pot=(source/'POTCAR').read_text();pot_sources=[source/'POTCAR']
            note='Original ML/input POSCAR retained for SPE; existing native input/last usable geometry for relaxation.'
        else:
            geo=Path(geometry[s+'|'+m]['structures']['mlip']['source_path']);atoms=read(geo)
            nmetal=atoms.get_chemical_symbols().count(metal)
            existing=sorted(groups[(s,m,f)],key=lambda r:r['cluster']!='perlmutter')
            expected=next((r['calculation'].get('potentials',{}) for r in existing if r['calculation'].get('potentials')), {})
            templates=[r for r in local if r['surface']==s and r['calculation']['functional']==f and r['calculation'].get('cell')
                       and r['calculation']['composition'].get(metal)==nmetal
                       and np.allclose(atoms.cell,r['calculation']['cell'],atol=1e-5,rtol=0)
                       and (not expected.get(metal) or r['calculation']['potentials'].get(metal)==expected[metal])]
            templates.sort(key=lambda r:(r['calculation']['status']!='converged',r['role']!='complex',len(r['calculation']['directory']),r['calculation']['directory']))
            if not templates:raise ValueError('No matching surface template for '+str(plan['key']))
            template=templates[0];source=Path(template['calculation']['directory']);settings=incar_values((source/'INCAR').read_text())
            if any(k in settings for k in ['MAGMOM','NELECT','NUPDOWN','LDAUL','LDAUU','LDAUJ']):raise ValueError('Structure-dependent template settings need review')
            nfixed=bottom_constraints(atoms,read(source/'POSCAR',format='vasp'),metal)
            order=sorted(range(len(atoms)),key=lambda i:(atoms[i].symbol!=metal,atoms[i].symbol,i));atoms=atoms[order]
            species=list(dict.fromkeys(atoms.get_chemical_symbols()));parts=[];pot_sources=[]
            for element in species:
                title=expected.get(element) or template['calculation']['potentials'].get(element)
                if not title:
                    candidates=[r for r in local if r['role']=='molecule' and element in r['calculation'].get('potentials',{})]
                    candidates.sort(key=lambda r:(r['calculation']['status']!='converged',r['calculation']['functional']!=f,r['calculation']['directory']))
                    title=candidates[0]['calculation']['potentials'][element]
                text,path=block(title);parts.append(text);pot_sources.append(path)
            pot=''.join(parts)
            note=f'New DFT validation from published MLIP geometry; {nfixed} lower-layer metal atoms fixed using a template with the same cell and metal count. Not an asserted reproduction of unavailable Kestrel inputs.'
        drafts.append(dict(role=role,surface=s,molecule=m,functional=f,source=source,geometry_source=geo,atoms=atoms,pot=pot,pot_sources=pot_sources,settings=settings,targets=[plan['key']],note=note,why=plan['why'],copy_geometry=bool(plan['source'])))
    # Missing reference requests are deduplicated by their actual physical setup.
    needs=[]
    for plan in plans:
        if plan['missing_refs']:
            comp=next(r for r in records if r['calculation']['directory']==plan['completed_directory'])
            for role in plan['missing_refs']:needs.append((role,comp,plan['key']))
    # Also check references for newly prepared complexes that have missing energies.
    for d in drafts:
        if 'energy' not in d['why']:continue
        potentials={re.search(r'VRHFIN\s*=\s*([A-Z][a-z]?)\s*:',b)[1]:re.search(r'TITEL\s*=\s*([^\n]+)',b)[1].strip() for b in re.findall(r'.*?End of Dataset[^\n]*\n?',d['pot'],re.S)}
        pseudo=dict(snapshot_id='planned',surface=d['surface'],molecule=d['molecule'],calculation=dict(directory=str(d['source']),functional=d['functional'],composition=dict(Counter(d['atoms'].get_chemical_symbols())),potentials=potentials,cell=d['atoms'].cell.tolist()))
        for role in reference_candidates(pseudo,records,allow_unverified=True)['unresolved_roles']:needs.append((role,pseudo,d['targets'][0]))
    refs={}
    for role,comp,target in needs:
        s,m,f,_=target;metal=re.match(r'[A-Z][a-z]?',s)[0];c=comp['calculation']
        composition={k:v for k,v in c['composition'].items() if (k==metal)==(role=='slab')}
        potmap={k:c['potentials'][k] for k in composition}
        key=json.dumps([role,s if role=='slab' else m,f,composition,c.get('cell') if role=='slab' else None,potmap],sort_keys=True)
        if key in refs:refs[key]['targets'].append(target);continue
        if role=='molecule':
            candidates=[r for r in local if r['role']=='molecule' and r['molecule']==m and r['calculation']['functional']==f and r['calculation'].get('composition')==composition and all(r['calculation'].get('potentials',{}).get(k)==v for k,v in potmap.items())]
            candidates.sort(key=lambda r:(len(r['calculation']['directory']),r['calculation']['directory']))
            if not candidates:raise ValueError('No gas reference inputs for '+key)
            source=Path(candidates[0]['calculation']['directory']);settings=incar_values((source/'INCAR').read_text());geo=source/'POSCAR'
            if candidates[0]['calculation']['status']!='electronic_unconverged' and (source/'CONTCAR').is_file():
                try:read(source/'CONTCAR',format='vasp');geo=source/'CONTCAR'
                except Exception:pass
            atoms=read(geo,format='vasp');pot=(source/'POTCAR').read_text();pot_sources=[source/'POTCAR'];copy_geometry=True
        else:
            directory=c.get('directory');source=Path(directory) if directory else None
            if not source or not (source/'POTCAR').exists():raise ValueError('No local complex for slab reference '+key)
            settings=incar_values((source/'INCAR').read_text());geo=source/'POSCAR';a=read(geo,format='vasp');atoms=a[[i for i,x in enumerate(a) if x.symbol==metal]]
            pot,pot_source=block(potmap[metal]);pot_sources=[pot_source];copy_geometry=False
        refs[key]=dict(role=role,surface=s,molecule=m,functional=f,source=source,geometry_source=geo,atoms=atoms,pot=pot,pot_sources=pot_sources,settings=settings,targets=[target],why=['reference'],copy_geometry=copy_geometry,
                       note='Recover missing reference with matching atom count and recorded potentials; gas identity retained; slab cell and constraints inherited from complex.')
    drafts=list(refs.values())+drafts
    drafts.sort(key=lambda d:({'molecule':0,'slab':1,'spe':2,'relaxed':3}[d['role']],-len(d['targets']),d['surface'],d['molecule'],d['functional']))
    campaign.mkdir(parents=True);(campaign/'logs').mkdir();jobs=[]
    for i,d in enumerate(drafts):
        role=d['role'];f=d['functional'];s=d['surface'];m=d['molecule'];source=d['source'];tag=hashlib.sha256(json.dumps(d['targets'],sort_keys=True).encode()).hexdigest()[:8]
        label=m if role=='molecule' else s+'_'+tag if role=='slab' else m+'_'+s
        directory=campaign/role/label/FUNCS[f]
        # Multiple potential variants are explicit nested groups for gas references.
        if directory.exists():directory=directory/('variant_'+tag)
        directory.mkdir(parents=True)
        atoms=d['atoms'];settings=settings_for(d['settings'],f,role)
        settings['SYSTEM']=label
        if d['copy_geometry']:shutil.copyfile(d['geometry_source'],directory/'POSCAR')
        else:write(directory/'POSCAR',atoms,format='vasp',direct=True,vasp5=True)
        (directory/'POTCAR').write_text(d['pot'])
        (directory/'INCAR').write_text(''.join(k+' = '+v+'\n' for k,v in settings.items()))
        for name in ['KPOINTS','vdw_kernel.bindat']:
            if (source/name).is_file():shutil.copyfile(source/name,directory/name)
        if f=='beef_vdw' and not (directory/'vdw_kernel.bindat').exists():raise ValueError('Template lacks vdW kernel')
        validate_inputs(directory)
        sources=set([d['geometry_source'],source/'INCAR',source/'KPOINTS',*d['pot_sources']])
        if (source/'vdw_kernel.bindat').exists():sources.add(source/'vdw_kernel.bindat')
        hashes={name:sha(directory/name) for name in ['POSCAR','INCAR','POTCAR','KPOINTS']}
        if (directory/'vdw_kernel.bindat').exists():hashes['vdw_kernel.bindat']=sha(directory/'vdw_kernel.bindat')
        changes={k:dict(before=d['settings'].get(k),after=settings.get(k)) for k in sorted(set(settings)|set(d['settings'])) if settings.get(k)!=d['settings'].get(k)}
        job=dict(array_index=i,role=role,system=label,surface=s,molecule=m,functional=f,directory=str(directory),source_directory=str(source),geometry_source=str(d['geometry_source']),
                 source_files_sha256={str(p):sha(p) for p in sorted(sources)},staged_sha256=hashes,changes=changes,natoms=len(atoms),cell=atoms.cell.tolist(),constraints=[x.todict() for x in atoms.constraints],
                 target_cells=d['targets'],target_systems=sorted({x[1]+'_'+x[0] for x in d['targets']}),note=d['note'],status='prepared')
        (directory/'completion_inputs.json').write_text(json.dumps(job,indent=2,default=lambda x:x.tolist())+'\n')
        jobs.append(job)
    manifest=dict(schema_version=1,created_at=datetime.now(timezone.utc).isoformat(),project_root=str(project),jobs=jobs,
                  live_snapshot_sha256=sha(live_snapshot),plan_summary=dict(Counter(p['status'] for p in plans)),
                  policy='Fill blank energy entries and absent functional geometries. Preserve originals and existing 169-job restart campaign. New unavailable-native geometries start from archived MLIP structures with cell/count-matched lower-layer constraints and recorded PAW variants. SPE remains fixed geometry. Gas/slab references deduplicated. No calculation is labelled converged until output verification.',
                  kestrel_note='Archived Kestrel components retained; live Kestrel authentication unavailable.')
    (campaign/'manifest.json').write_text(json.dumps(manifest,indent=2,default=lambda x:x.tolist())+'\n')
    (campaign/'joblist.txt').write_text(''.join(j['directory']+'\n' for j in jobs))
    (campaign/'gap_plan.json').write_text(json.dumps(plans,indent=2)+'\n')
    shutil.copyfile(live_snapshot,campaign/'preflight_components.jsonl.gz')
    shutil.copyfile(root/'scripts/completion_array.slurm',campaign/'run_array.slurm')
    print(json.dumps(dict(jobs=len(jobs),roles=dict(Counter(j['role'] for j in jobs)),campaign=str(campaign)),indent=2))
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','project','campaign','plans','records','live-snapshot']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.root.resolve(),a.project.resolve(),a.campaign.resolve(),json.loads(a.plans.read_text()),json.loads(a.records.read_text()),a.live_snapshot)
