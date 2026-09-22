#!/usr/bin/env python3
"""Prepare two isolated DME gas-reference recoveries; leave original runs intact."""
import json
import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil
from collections import Counter
from ase.io import read
from prepare_completion_jobs import sha, incar_text
from extract_kestrel_singlepoint import component

ROOT=Path('/pscratch/sd/j/jcho5/VASP')
CAMPAIGN=ROOT/'completion_runs/dme_references_20260915'


def main():
    global CAMPAIGN
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',type=Path,default=CAMPAIGN)
    parser.add_argument('--ialgo',choices=['38','48'],default='38')
    args=parser.parse_args();CAMPAIGN=args.campaign
    if (CAMPAIGN/'manifest.json').exists():
        print('Already prepared:',CAMPAIGN);return
    geometry=ROOT/'vasp_mol/CH3OCH3/beef_vdw'
    assert component(geometry)['status']=='converged'
    atoms=read(geometry/'CONTCAR',format='vasp')
    assert Counter(atoms.get_chemical_symbols())=={'C':2,'O':1,'H':6}
    CAMPAIGN.mkdir(parents=True,exist_ok=False);(CAMPAIGN/'logs').mkdir();jobs=[]
    for folder,functional in [('PBE','pbe'),('PBE_D3','pbe_d3')]:
        source=ROOT/'vasp_mol/CH3OCH3'/folder
        assert component(source)['status']!='converged'
        directory=CAMPAIGN/'molecule/DME'/folder;directory.mkdir(parents=True)
        shutil.copyfile(geometry/'CONTCAR',directory/'POSCAR')
        for name in ['POTCAR','KPOINTS']:shutil.copyfile(source/name,directory/name)
        changes=dict(IBRION='1',POTIM='0.2',NSW='500',NELM='200',EDIFF='1E-06',ISTART='0',LWAVE='.FALSE.',LCHARG='.FALSE.',IALGO=args.ialgo)
        (directory/'INCAR').write_text(incar_text((source/'INCAR').read_text(),changes))
        targets=sorted(p.name for p in (ROOT/'dft_jobs').glob('CH3OCH3_*') if p.is_dir())
        job=dict(role='molecule',molecule='DME',system='DME',functional=functional,surface='',directory=str(directory),
            source_directory=str(source),geometry_source=str(geometry/'CONTCAR'),geometry_source_sha256=sha(geometry/'CONTCAR'),
            target_systems=targets,natoms=len(atoms),cell=atoms.cell.tolist(),changes=changes,
            staged_sha256={name:sha(directory/name) for name in ['INCAR','POSCAR','POTCAR','KPOINTS']},status='prepared')
        (directory/'completion_inputs.json').write_text(json.dumps(job,indent=2)+'\n');jobs.append(job)
    manifest=dict(created_at=datetime.now(timezone.utc).isoformat(),project_root=str(ROOT),jobs=jobs,
        policy='Two gas-only PBE/PBE+D3 retries from converged BEEF DME coordinates; no BEEF energy substitution. Original potentials, cutoff, cell, k points and functional preserved. RMM-DIIS ionic optimization replaces failed conjugate-gradient line search; electronic IALGO='+args.ialgo+'; tighter electronic tolerance; original force threshold retained.')
    (CAMPAIGN/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (CAMPAIGN/'joblist.txt').write_text(''.join(j['directory']+'\n' for j in jobs))
    worker=(Path(__file__).parent/'completion_array.slurm').read_text()
    worker=worker.replace('#SBATCH -n 32','#SBATCH -n 16').replace('#SBATCH -t 06:00:00','#SBATCH -t 02:00:00')
    worker=worker.replace('#SBATCH -n 16','#SBATCH -N 1\n#SBATCH -n 16')
    worker=worker.replace('-n 32 -c 2 --mem=120G','-n 16 -c 2 --mem=60G')
    worker=worker.replace("if record['role']=='slab':","if record['role'] in ('slab','molecule'):")
    (CAMPAIGN/'worker.slurm').write_text(worker)
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
