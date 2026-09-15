"""Discover only jobs recorded in local, isolated completion campaign manifests."""
import json
import hashlib
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
    """Restrict a supplemental reference to its intended system and unchanged cell."""
    if job['role'] != 'slab' or job['functional'] != functional or system not in job['target_systems']:
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
        return bool(np.allclose(slab.cell,complex_atoms.cell,atol=1e-5,rtol=0))
    except (OSError,ValueError):
        return False
