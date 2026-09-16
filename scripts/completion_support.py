"""Discover only jobs recorded in local, isolated completion campaign manifests."""
import json
import hashlib
import re
from collections import Counter
from pathlib import Path


def completion_jobs(project):
    rows=[]
    for manifest in sorted((Path(project)/'completion_runs').glob('*/manifest.json')):
        record=json.loads(manifest.read_text())
        for job in record.get('jobs',[]):
            directory=Path(job['directory']).resolve()
            if manifest.parent.resolve() not in directory.parents:
                raise ValueError('Completion job escapes its campaign: '+str(directory))
            rows.append(job)
    return rows


def matching_slab(job, system, functional, complex_directory):
    """Reuse a completed clean slab only for the same surface, composition and cell."""
    match=re.fullmatch(r'.+_([A-Z][a-z]?\d+)(?:_.+)?',system)
    if job['role'] != 'slab' or job['functional'] != functional or not match or match.group(1) != job.get('surface'):
        return False
    import numpy as np
    from ase.io import read
    directory=Path(job['directory'])
    try:
        for name,digest in job['staged_sha256'].items():
            if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
                return False
        slab=read(directory/'POSCAR',format='vasp')
        complex_atoms=read(Path(complex_directory)/'POSCAR',format='vasp')
        metal=re.match(r'[A-Z][a-z]?',job['surface']).group()
        if Counter(slab.get_chemical_symbols()) != {metal:complex_atoms.get_chemical_symbols().count(metal)}:
            return False
        return bool(np.allclose(slab.cell,complex_atoms.cell,atol=1e-5,rtol=0))
    except (OSError,ValueError):
        return False
