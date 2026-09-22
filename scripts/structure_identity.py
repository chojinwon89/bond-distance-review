"""Conservative structure identities for duplicate selection, independent of energy."""
import hashlib
import json
from pathlib import Path
import numpy as np
from ase.io import read


def geometry_record(path):
    path=Path(path);atoms=read(path,format='vasp')
    if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell).all():
        raise ValueError('Nonfinite geometry')
    # Cartesian coordinates at 0.0001 Å precision, fixed atom order and constraints.
    # This intentionally does not infer equivalence from a formula or cell alone.
    constraints=[c.todict() for c in atoms.constraints]
    def serial(value):
        if hasattr(value,'tolist'):return value.tolist()
        raise TypeError(type(value).__name__)
    data=dict(symbols=atoms.get_chemical_symbols(),cell=np.round(atoms.cell.array,6).tolist(),
              positions=np.round(atoms.positions,4).tolist(),constraints=constraints)
    payload=json.dumps(data,sort_keys=True,separators=(',',':'),default=serial)
    return dict(file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                geometry_sha256=hashlib.sha256(payload.encode()).hexdigest(),
                convention='ordered symbols + cell rounded 1e-6 A + Cartesian positions rounded 1e-4 A + constraints',
                path=str(path),natoms=len(atoms))


def attach_geometry(result,directory):
    directory=Path(directory); result=dict(result); geometry={}
    for role,name in [('input','POSCAR'),('relaxed','CONTCAR')]:
        path=directory/name
        if path.is_file() and (role=='input' or result.get('status')=='converged' and (result.get('nsw') or 0)>0):
            try:geometry[role]=geometry_record(path)
            except (OSError,ValueError,IndexError):pass
    if geometry:result['structure_identity']=geometry
    return result


def same_structure(left,right,mode):
    a,b=[r.get('calculation',r) for r in [left,right]]
    role='input' if mode=='SPE' else 'relaxed'
    ga,gb=[r.get('structure_identity',{}).get(role,{}) for r in [a,b]]
    for key in ['file_sha256','geometry_sha256']:
        if ga.get(key) and ga.get(key)==gb.get(key):return True,key
    if mode=='SPE':
        def fingerprint(r):return r.get('structure_identity',{}).get('input',{}).get('file_sha256') or r.get('geometry_provenance',{}).get('poscar_sha256')
        if fingerprint(a) and fingerprint(a)==fingerprint(b):return True,'input_POSCAR_sha256'
    return False,'structure fingerprint missing or different'
