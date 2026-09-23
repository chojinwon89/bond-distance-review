"""Archive actual coordinates and atom-resolved adsorption geometry for all functionals."""
import csv, hashlib, itertools, json, math, re
from collections import Counter, defaultdict
from pathlib import Path
import networkx as nx
import numpy as np
from ase.io import read
from ase.data import covalent_radii, atomic_numbers
from ase.geometry import find_mic
from bs4 import BeautifulSoup
from component_store import latest, read_store
from molecule_names import canonical, equivalent_names, ALIASES

ROOT=Path(__file__).resolve().parents[1]
FUNCS=['pbe','pbe_d3','r2scan','beef_vdw']


def graph(atoms, indices):
    g=nx.Graph()
    for i in indices:g.add_node(i,element=atoms[i].symbol)
    for i,j in itertools.combinations(indices,2):
        d=atoms.get_distance(i,j,mic=True)
        if 0.35<d<1.25*(covalent_radii[atoms[i].number]+covalent_radii[atoms[j].number]):g.add_edge(i,j)
    return g


def correspondence(reference, target, metal):
    """Map heavy-atom molecular graphs; resolve symmetry by centered RMSD, no index assumption."""
    a=[i for i,x in enumerate(reference) if x.symbol not in (metal,'H')]
    b=[i for i,x in enumerate(target) if x.symbol not in (metal,'H')]
    if not a or len(a)!=len(b):return {},'No heavy-atom correspondence'
    ga,gb=graph(reference,a),graph(target,b)
    matcher=nx.algorithms.isomorphism.GraphMatcher(ga,gb,node_match=lambda x,y:x['element']==y['element'])
    best=None;count=0
    ar=reference.positions[a];ar=find_mic(ar-ar[0],reference.cell,reference.pbc)[0];ar-=ar.mean(axis=0)
    for mapping in itertools.islice(matcher.isomorphisms_iter(),256):
        br=target.positions[[mapping[i] for i in a]];br=find_mic(br-br[0],target.cell,target.pbc)[0];br-=br.mean(axis=0)
        u,_,vt=np.linalg.svd(br.T@ar);fix=np.eye(3);fix[-1,-1]=np.linalg.det(u@vt)
        score=float(np.sum((br@u@fix@vt-ar)**2));count+=1
        item=(score,tuple(mapping[i] for i in a),mapping)
        if best is None or item[:2]<best[:2]:best=item
    if best is None:return {},'Molecular connectivity differs; atoms are not paired'
    mapping=dict(best[2])
    # Hydrogens are paired only within the same mapped covalent attachment.
    from scipy.optimize import linear_sum_assignment
    full_a=graph(reference,[i for i,x in enumerate(reference) if x.symbol!=metal])
    full_b=graph(target,[i for i,x in enumerate(target) if x.symbol!=metal])
    for i,j in list(mapping.items()):
        ha=[h for h in full_a[i] if reference[h].symbol=='H' and len(full_a[h])==1]
        hb=[h for h in full_b[j] if target[h].symbol=='H' and len(full_b[h])==1]
        if not ha or len(ha)!=len(hb):continue
        va=reference.get_distances(i,ha,mic=True,vector=True);vb=target.get_distances(j,hb,mic=True,vector=True)
        va/=np.linalg.norm(va,axis=1)[:,None];vb/=np.linalg.norm(vb,axis=1)[:,None]
        rows,cols=linear_sum_assignment(np.sum((va[:,None,:]-vb[None,:,:])**2,axis=2))
        mapping.update({ha[a]:hb[b] for a,b in zip(rows,cols)})
    return mapping,f'Heavy-atom graph match; symmetry choice from {count} mappings; attached H paired by bond direction'+(' (search capped)' if count==256 else '')


def measure(atoms, metal):
    ads=[i for i,x in enumerate(atoms) if x.symbol!=metal];slab=[i for i,x in enumerate(atoms) if x.symbol==metal]
    if not ads or not slab:raise ValueError('Missing adsorbate or metal atoms')
    if not np.isfinite(atoms.positions).all():raise ValueError('Nonfinite coordinates')
    labels={i:atoms[i].symbol+str(k+1) for k,i in enumerate(ads)}
    contacts=[]
    for i in ads:
        distances=atoms.get_distances(i,slab,mic=True);j=slab[int(np.argmin(distances))]
        cutoff=1.25*(covalent_radii[atoms[i].number]+covalent_radii[atoms[j].number])
        contacts.append(dict(atom=i,label=labels[i],element=atoms[i].symbol,metal_atom=j,
            metal_label=metal+str(j+1),distance_A=round(float(min(distances)),6),contact_candidate=bool(min(distances)<=cutoff)))
    contacts.sort(key=lambda c:c['distance_A'])
    heavy=[i for i in ads if atoms[i].symbol!='H'];g=graph(atoms,ads)
    angles=[];torsions=[]
    for center in g:
        neighbors=sorted(g[center]);heavy_neighbors=[i for i in neighbors if atoms[i].symbol!='H']
        if atoms[center].symbol=='C' and len(heavy_neighbors)>=2:neighbors=heavy_neighbors
        for a,b in itertools.combinations(neighbors,2):
            ids=[a,center,b]
            try:value=float(atoms.get_angle(*ids,mic=True))
            except (ZeroDivisionError,ValueError):continue
            angles.append(dict(atoms=ids,label='–'.join(labels[i] for i in ids),degrees=round(value,4)))
    for b,c in sorted(g.edges()):
        if atoms[b].symbol=='H' or atoms[c].symbol=='H':continue
        left=sorted(set(g[b])-{c});right=sorted(set(g[c])-{b})
        left=[i for i in left if atoms[i].symbol!='H'] or left[:1]
        right=[i for i in right if atoms[i].symbol!='H'] or right[:1]
        for a in left:
            for d in [i for i in right if i!=a]:
                ids=[a,b,c,d]
                try:value=float(atoms.get_dihedral(*ids,mic=True));value=(value+180)%360-180
                except (ZeroDivisionError,ValueError):continue
                if math.isfinite(value):torsions.append(dict(atoms=ids,label='–'.join(labels[i] for i in ids),degrees=round(value,4)))
    return dict(adsorbate_indices=ads,labels=labels,contacts=contacts,angles=angles,torsions=torsions,
        nearest_contact_A=contacts[0]['distance_A'],nearest_heavy_contact_A=next((c['distance_A'] for c in contacts if c['element']!='H'),None))


def serialize(atoms, path, cluster, status, snapshot=None, metal=None):
    if metal is None:raise ValueError('Surface metal must be explicit')
    ads=[i for i,a in enumerate(atoms) if a.symbol!=metal]
    center=atoms.positions[ads[0]]
    positions=find_mic(atoms.positions-center,atoms.cell,atoms.pbc)[0]
    return dict(display_positions=np.round(positions,7).tolist(),molecular_bonds=list(graph(atoms,ads).edges()),symbols=atoms.get_chemical_symbols(),positions=np.round(atoms.positions,7).tolist(),cell=np.round(atoms.cell.array,7).tolist(),
        pbc=atoms.pbc.tolist(),source_path=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cluster=cluster,status=status,snapshot_id=snapshot)


def build(root=ROOT):
    page=(root/'dft_comparison.html').read_text();soup=BeautifulSoup(page,'html.parser')
    cache_path=root/'data/adsorption_results/coordinate_sources.json'
    coordinate_cache=json.loads(cache_path.read_text()) if cache_path.exists() else {}
    current=latest(read_store(root/'dft_component_store.jsonl.gz'));groups=defaultdict(list)
    for r in current:
        c=r['calculation']
        if r['role']=='complex' and (c.get('nsw') or 0)>0:groups[(r['surface'],canonical(r['molecule']),c['functional'])].append(r)
    with (root/'dft_cluster_selection.csv').open() as f:
        selected={(r['surface'],r['molecule'],r['functional']):r for r in csv.DictReader(f) if r['mode']=='relaxed'}
    data={};sources=[];errors=[];counts=Counter()
    for card in soup.select('.g'):
        surface,molecule=card['data-surf'],canonical(card['data-mol']);metal=re.match(r'[A-Z][a-z]?',surface)[0];key=surface+'|'+molecule
        record=dict(surface=surface,molecule=molecule,metal=metal,structures={});data[key]=record
        base='published base' in card.select_one('.m').get_text()
        names=[f'{surface}_{card["data-mol"]}.cif',f'{surface}_{card["data-mol"]}_sevennet_omni.cif'] if base else [f'{surface}_{card["data-mol"]}_sevennet_omni.cif',f'{surface}_{card["data-mol"]}.cif']
        names.insert(0,f'{surface}_{card["data-mol"]}_5m.cif') if base else names.append(f'{surface}_{card["data-mol"]}_5m.cif')
        for alias in equivalent_names(molecule,ALIASES):
            for suffix in ['_sevennet_omni.cif','_5m.cif','.cif']:
                name=f'{surface}_{alias}{suffix}'
                if name not in names:names.append(name)
        ml=None
        for name in names:
            path=root/name
            if not path.exists():continue
            try:
                ml=read(path);metrics=measure(ml,metal)
                record['structures']['mlip']=dict(serialize(ml,path,'published MLIP','relaxed',metal=metal),**metrics,model='SevenNet-OMNI' if 'sevennet_omni' in name else 'published GOAD 5-model baseline' if '_5m' in name else 'published GOAD baseline')
                break
            except Exception as e:errors.append(dict(path=str(path),reason=str(e)))
        for functional in FUNCS:
            target=selected.get((surface,molecule,functional),{});candidates=groups.get((surface,molecule,functional),[])
            candidates.sort(key=lambda r:(r['calculation']['status']!='converged',r['cluster']!='perlmutter',r['calculation']['directory']!=target.get('complex_directory'),len(r['calculation']['directory']),r['calculation']['directory']))
            result=None
            for candidate in candidates:
                c=candidate['calculation'];path=Path(c['directory'])/'CONTCAR'
                if not path.is_file() or not path.stat().st_size:
                    cached=coordinate_cache.get(candidate['cluster']+'|'+c['directory'])
                    if not cached:continue
                    path=root/cached['path']
                    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=cached['sha256']:continue
                try:
                    atoms=read(path,format='vasp')
                    if dict(Counter(atoms.get_chemical_symbols()))!=c['composition']:raise ValueError('CONTCAR composition does not match archived calculation')
                    metrics=measure(atoms,metal)
                    result=dict(serialize(atoms,path,candidate['cluster'],c['status'],candidate['snapshot_id'],metal=metal),**metrics,
                        matches_energy_complex=c['directory']==target.get('complex_directory'),functional=functional)
                    if ml is not None:
                        mapping,note=correspondence(ml,atoms,metal);result['ml_to_dft_atom_map']=mapping;result['atom_mapping_note']=note
                        pairs=[];lookup={x['atom']:x for x in result['contacts']}
                        for contact in record['structures']['mlip']['contacts']:
                            if contact['element']!='H' and contact['atom'] in mapping:
                                dft=lookup[mapping[contact['atom']]]
                                pairs.append(dict(ML_atom=contact['label'],DFT_atom=dft['label'],ML_distance_A=contact['distance_A'],DFT_distance_A=dft['distance_A'],delta_A=round(dft['distance_A']-contact['distance_A'],6),ML_metal=contact['metal_label'],DFT_metal=dft['metal_label']))
                        result['paired_contacts']=pairs
                        surface_angles=[];ml_graph=graph(ml,[i for i in mapping if ml[i].symbol!='H'])
                        ml_contacts={x['atom']:x for x in record['structures']['mlip']['contacts']}
                        for anchor in ml_graph:
                            for neighbor in sorted(ml_graph[anchor]):
                                mid=ml_contacts[anchor]['metal_atom'];did=lookup[mapping[anchor]]['metal_atom']
                                try:
                                    va=float(ml.get_angle(mid,anchor,neighbor,mic=True))
                                    vb=float(atoms.get_angle(did,mapping[anchor],mapping[neighbor],mic=True))
                                except (ValueError,ZeroDivisionError):continue
                                surface_angles.append(dict(label=metal+str(mid+1)+'–'+record['structures']['mlip']['labels'][anchor]+'–'+record['structures']['mlip']['labels'][neighbor],
                                    DFT_label=metal+str(did+1)+'–'+result['labels'][mapping[anchor]]+'–'+result['labels'][mapping[neighbor]],
                                    ML_degrees=round(va,4),DFT_degrees=round(vb,4),delta_degrees=round(vb-va,4)))
                        result['paired_surface_angles']=surface_angles
                        for kind in ['angles','torsions']:
                            paired=[]
                            for entry in record['structures']['mlip'][kind]:
                                if not all(i in mapping for i in entry['atoms']):continue
                                ids=[mapping[i] for i in entry['atoms']]
                                try:v=float(atoms.get_angle(*ids,mic=True) if kind=='angles' else atoms.get_dihedral(*ids,mic=True))
                                except (ValueError,ZeroDivisionError):continue
                                if not math.isfinite(v):continue
                                if kind=='torsions':v=(v+180)%360-180
                                delta=v-entry['degrees']
                                if kind=='torsions':delta=(delta+180)%360-180
                                paired.append(dict(label=entry['label'],DFT_label='–'.join(result['labels'][i] for i in ids),ML_degrees=entry['degrees'],DFT_degrees=round(v,4),delta_degrees=round(delta,4)))
                            result['paired_'+kind]=paired
                    record['structures'][functional]=result;counts['DFT '+c['status']]+=1
                    sources.append(dict(surface=surface,molecule=molecule,functional=functional,source=str(path),cluster=candidate['cluster'],status=c['status'],matches_energy_complex=result['matches_energy_complex'],sha256=result['source_sha256']))
                    break
                except Exception as e:errors.append(dict(path=str(path),reason=str(e)))
            if result is None:counts['DFT missing coordinates']+=1
        counts['MLIP available' if ml is not None else 'MLIP missing']+=1
    folder=root/'data/adsorption_results';folder.mkdir(exist_ok=True,parents=True)
    # Content-addressed coordinate snapshots cannot be overwritten by a later refresh.
    payload=json.dumps(dict(schema_version=1,systems=data),separators=(',',':'),allow_nan=False).encode()
    digest=hashlib.sha256(payload).hexdigest();archive=folder/'geometry_snapshots';archive.mkdir(exist_ok=True)
    import gzip
    path=archive/(digest+'.json.gz')
    if not path.exists():path.write_bytes(gzip.compress(payload,mtime=0))
    (root/'functional_geometry.json').write_bytes(payload)
    report=dict(structures=dict(counts),source_count=len(sources),coordinate_snapshot=str(path.relative_to(root)),sha256=digest,
        distance_convention='ASE minimum-image distances using source cell and periodicity. Contacts are geometric candidates, not a bond-order calculation.',
        atom_labels='Molecular labels use 1-based adsorbate order; metal labels use 1-based source atom index. Cross-model heavy-atom pairing uses molecular graph isomorphism, with symmetry resolved by centered rigid-fit RMSD; attached hydrogens use bond-direction matching within a mapped attachment; changed connectivity is not paired.',
        angles='Molecular bond paths, including heteroatom-bound H and terminal H torsions when heavy-only paths are unavailable; angle at the middle atom. Torsions in [-180,180) degrees; wrapped DFT-minus-ML differences.',
        kestrel_access='Direct SSH was reachable but authentication was unavailable. Archived exports and locally present coordinates only.',errors=errors,sources=sources)
    (root/'functional_geometry_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    rows=[]
    for sys in data.values():
        for f,s in sys['structures'].items():
            for c in s['contacts']:rows.append(dict(surface=sys['surface'],molecule=sys['molecule'],method=f,status=s['status'],source_path=s['source_path'],**c))
    with (root/'functional_geometry_contacts.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # External enhancement leaves the seven energy cells stable for the data pipeline.
    if 'functional_geometry.js' not in page:
        page=page.replace('</head>','<link rel="stylesheet" href="functional_geometry.css?v=static-ase">\n</head>')
        page=page.replace('</body>','<script src="functional_geometry.js?v=static-ase" defer></script>\n</body>')
    page=re.sub(r'<img(?![^>]*\bloading=)', '<img loading="lazy"',page)
    (root/'dft_comparison.html').write_text(page)
    print(json.dumps(dict(counts),indent=2),flush=True)


if __name__=='__main__':build()
