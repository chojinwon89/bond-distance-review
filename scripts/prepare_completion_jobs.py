#!/usr/bin/env python3
"""Stage isolated jobs for five specific missing page entries; never edit source jobs."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil

from ase.io import read, write
from add_missing_dft_images import incar_values
from extract_kestrel_singlepoint import component


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def incar_text(original, overrides):
    values = incar_values(original)
    values.update(overrides)
    return '\n'.join(f'{k} = {v}' for k, v in values.items()) + '\n'


def metal_potential(text, metal):
    blocks = re.findall(r'.*?End of Dataset[^\n]*\n?', text, re.S)
    matches = [b for b in blocks if re.search(r'VRHFIN\s*=\s*'+metal+r'\s*:', b)]
    if len(matches) != 1:
        raise ValueError('Cannot identify exactly one metal POTCAR block')
    return matches[0]


def prepare(project, campaign):
    if (campaign/'manifest.json').exists():
        return json.loads((campaign/'manifest.json').read_text())
    campaign.mkdir(parents=True, exist_ok=False)
    jobs = []
    specs = [('slab', 'C2H4_Cu100', 'beef_vdw', ['C2H4_Cu100','C2H6_Cu100']),
             ('slab', 'C2H4_Ir100', 'r2scan', ['C2H4_Ir100']),
             ('spe', 'C2H6_Ir100', 'r2scan', ['C2H6_Ir100']),
             ('spe', 'C2H6_Rh100', 'beef_vdw', ['C2H6_Rh100'])]
    for role, system, functional, targets in specs:
        source = project/'dft_jobs'/system/'singlepoint'/functional
        result = component(source)
        if role == 'spe' and result['status'] != 'electronic_unconverged':
            raise ValueError('Retry requires a finished electronically unconverged SPE: '+str(source))
        if role == 'slab' and result['status'] != 'converged':
            raise ValueError('Reference source SPE must be completed: '+str(source))
        surface = re.search(r'_([A-Z][a-z]?\d+)', system).group(1)
        directory = campaign/role/system/functional
        directory.mkdir(parents=True, exist_ok=False)
        atoms = read(source/'POSCAR', format='vasp')
        overrides = {'ISTART':'0', 'NELM':'300', 'LWAVE':'.FALSE.', 'LCHARG':'.FALSE.'}
        if role == 'slab':
            metal = re.match(r'[A-Z][a-z]?', surface).group()
            atoms = atoms[[i for i,s in enumerate(atoms.get_chemical_symbols()) if s==metal]]
            write(directory/'POSCAR', atoms, format='vasp', direct=True, vasp5=True)
            (directory/'POTCAR').write_text(metal_potential((source/'POTCAR').read_text(),metal))
            overrides.update(SYSTEM=surface, NSW='1000', IBRION='2', EDIFFG='-0.05')
        else:
            for name in ['POSCAR','POTCAR']:shutil.copyfile(source/name, directory/name)
            overrides.update(NSW='0', IBRION='-1')
        for name in ['KPOINTS','vdw_kernel.bindat']:
            if (source/name).exists():shutil.copyfile(source/name,directory/name)
        (directory/'INCAR').write_text(incar_text((source/'INCAR').read_text(),overrides))
        hashes = {name:sha(directory/name) for name in ['POSCAR','INCAR','POTCAR','KPOINTS']}
        if (directory/'vdw_kernel.bindat').exists():hashes['vdw_kernel.bindat']=sha(directory/'vdw_kernel.bindat')
        if role=='spe':
            record=json.loads((source/'spe_inputs.json').read_text())
            record.update(destination=str(directory),staged_sha256=hashes,retry_of=str(source),changes=overrides)
            (directory/'spe_inputs.json').write_text(json.dumps(record,indent=2)+'\n')
        job=dict(role=role,system=system,surface=surface,functional=functional,directory=str(directory),
                 source_directory=str(source),target_systems=targets,staged_sha256=hashes,changes=overrides,
                 natoms=len(atoms),cell=atoms.cell.tolist(),status='prepared')
        (directory/'completion_inputs.json').write_text(json.dumps(job,indent=2)+'\n')
        jobs.append(job)
    manifest=dict(created_at=datetime.now(timezone.utc).isoformat(),project_root=str(project),jobs=jobs,
        policy='Two isolated clean-slab relaxations and two finished-SPE electronic retries. Original jobs remain unchanged. Matching cell, atom count, constraints and metal potential; no energy scaling. NELM=300; original electronic tolerance retained.')
    (campaign/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (campaign/'joblist.txt').write_text(''.join(j['directory']+'\n' for j in jobs))
    (campaign/'logs').mkdir()
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--campaign',type=Path,required=True)
    args=p.parse_args()
    m=prepare(args.project_root.resolve(),args.campaign.resolve())
    print(json.dumps([{k:j[k] for k in ['role','system','functional','natoms','directory']} for j in m['jobs']],indent=2))
