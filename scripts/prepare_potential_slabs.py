"""Stage missing matched-potential slab references from audited local complexes."""
import csv,json,hashlib,re,shutil
from pathlib import Path
from datetime import datetime,timezone
from ase.io import read,write
from prepare_completion_jobs import sha,incar_text,metal_potential
from add_missing_dft_images import incar_values

ROOT=Path(__file__).resolve().parents[1]
PROJECT=Path('/pscratch/sd/j/jcho5/VASP')
CAMPAIGN=PROJECT/'completion_runs/matched_potential_slabs_20260915'

def prepare():
    if (CAMPAIGN/'manifest.json').exists():return
    with (ROOT/'dft_potential_mismatches.csv').open() as f:rows=list(csv.DictReader(f))
    groups={};unavailable=[]
    for row in sorted(rows,key=lambda r:(r['complex_cluster']!='perlmutter',r['complex_directory'])):
        if int(row['compatible_slab_count']):continue
        source=Path(row['complex_directory'])
        if row['complex_cluster']!='perlmutter' or not (source/'POTCAR').exists():
            unavailable.append(row);continue
        metal=re.match(r'[A-Z][a-z]?',row['surface'])[0]
        atoms=read(source/'POSCAR',format='vasp');atoms=atoms[[i for i,a in enumerate(atoms) if a.symbol==metal]]
        settings=incar_values((source/'INCAR').read_text())
        electronic={k:v for k,v in settings.items() if k not in ['SYSTEM','NSW','IBRION','POTIM','ISTART','LWAVE','LCHARG','EDIFFG','NWRITE','NELM','NELMIN']}
        pot=metal_potential((source/'POTCAR').read_text(),metal)
        key=json.dumps(dict(surface=row['surface'],functional=row['functional'],cell=atoms.cell.round(6).tolist(),natoms=len(atoms),
            constraints=[c.todict() for c in atoms.constraints],potential_sha256=hashlib.sha256(pot.encode()).hexdigest(),
            kpoints=(source/'KPOINTS').read_text(),settings=electronic),sort_keys=True,default=lambda v:v.tolist())
        if key not in groups:groups[key]=dict(source=source,atoms=atoms,pot=pot,row=row,targets=set())
        groups[key]['targets'].add(next(part for part in reversed(source.parts) if part.startswith(row['molecule']+'_') or part.startswith('CH3OCH3_')) if row['molecule']=='DME' else row['molecule']+'_'+row['surface'])
    CAMPAIGN.mkdir(parents=True,exist_ok=False);(CAMPAIGN/'logs').mkdir();jobs=[]
    for key,g in sorted(groups.items()):
        row=g['row'];source=g['source'];tag=hashlib.sha256(key.encode()).hexdigest()[:10]
        directory=CAMPAIGN/'slab'/(row['surface']+'_'+tag)/row['functional'];directory.mkdir(parents=True)
        write(directory/'POSCAR',g['atoms'],format='vasp',direct=True,vasp5=True)
        (directory/'POTCAR').write_text(g['pot'])
        for name in ['KPOINTS','vdw_kernel.bindat']:
            if (source/name).exists():shutil.copyfile(source/name,directory/name)
        changes=dict(SYSTEM=row['surface'],ISTART='0',NSW='500',IBRION='2',POTIM='0.3',EDIFFG='-0.05',NELM='200',IALGO='38',LWAVE='.FALSE.',LCHARG='.FALSE.')
        (directory/'INCAR').write_text(incar_text((source/'INCAR').read_text(),changes))
        job=dict(role='slab',system=row['surface']+'_'+tag,surface=row['surface'],functional=row['functional'],directory=str(directory),source_directory=str(source),
            target_systems=sorted(g['targets']),changes=changes,natoms=len(g['atoms']),cell=g['atoms'].cell.tolist(),status='prepared',
            metal_potential=json.loads(row['complex_potentials'])[re.match(r'[A-Z][a-z]?',row['surface'])[0]],
            staged_sha256={name:sha(directory/name) for name in ['POSCAR','INCAR','POTCAR','KPOINTS']})
        (directory/'completion_inputs.json').write_text(json.dumps(job,indent=2)+'\n');jobs.append(job)
    manifest=dict(created_at=datetime.now(timezone.utc).isoformat(),project_root=str(PROJECT),jobs=jobs,
        policy='Metal-only subset of original complex POSCAR; constraints and cell preserved. Exact metal POTCAR block and k mesh copied from complex. Deduplicate by surface, cell, atom count, constraints, potential hash, functional and electronic settings. No vdW parameter changes; consistent with existing complex calculations. Electronic Davidson solver, same original EDIFF. Originals untouched.',
        archived_only_candidates_without_local_inputs=len(unavailable))
    (CAMPAIGN/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (CAMPAIGN/'joblist.txt').write_text(''.join(j['directory']+'\n' for j in jobs))
    worker=(ROOT/'scripts/completion_array.slurm').read_text().replace('#SBATCH -n 32','#SBATCH -N 1\n#SBATCH -n 32').replace('#SBATCH -t 06:00:00','#SBATCH -t 02:00:00')
    (CAMPAIGN/'worker.slurm').write_text(worker)
    (ROOT/'dft_potential_slab_recovery.json').write_text(json.dumps(manifest,indent=2)+'\n');print('Prepared',len(jobs),'slab references at',CAMPAIGN)

if __name__=='__main__':prepare()
