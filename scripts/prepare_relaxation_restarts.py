#!/usr/bin/env python3
"""Stage audited unfinished relaxations in a new, immutable-input campaign.

This command does not submit jobs. Source files are only read. The audit must
contain full OUTCAR hashes and results of the convergence/error-marker scan.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read

from add_missing_dft_images import incar_values
from component_store import identity, latest, read_store, relaxed
from extract_kestrel_singlepoint import component
from prepare_completion_jobs import sha


def restart_settings(settings, zbrent=False, scf_warning=False):
    values = dict(settings)
    values.update(ISTART='0', ICHARG='2', LWAVE='.FALSE.', LCHARG='.FALSE.')
    if zbrent:
        values.update(EDIFF=str(min(float(values['EDIFF']), 1e-6)),
                      IBRION='1', POTIM='0.2')
    if scf_warning:
        values.pop('IALGO', None)
        values.update(ALGO='Normal', NELM=str(max(int(values['NELM']), 300)))
    return values


def validate_coordinates(source):
    original = read(source/'POSCAR', format='vasp')
    final = read(source/'CONTCAR', format='vasp')
    if final.get_chemical_symbols() != original.get_chemical_symbols():
        raise ValueError('Atom identity/order changed: '+str(source))
    if not np.isfinite(final.positions).all() or not np.allclose(final.cell, original.cell, atol=1e-6, rtol=0):
        raise ValueError('Non-finite coordinates or changed cell: '+str(source))
    # VASP selective dynamics must survive copying CONTCAR to POSCAR.
    def constraints(atoms):
        return json.dumps([c.todict() for c in atoms.constraints], sort_keys=True,
                          default=lambda x:x.tolist())
    if constraints(final) != constraints(original):
        raise ValueError('Constraints changed: '+str(source))
    distances = final.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    if distances.min() < 0.5:
        raise ValueError('Near-coincident atoms in restart coordinates: '+str(source))
    return final


def completed_alternatives(rows, store):
    targets = {(r['surface'], r['molecule'], r['functional']) for r in rows}
    # Conservatively stop for any same-system completed relaxation; the caller
    # must compare input fingerprints before deciding whether to reuse or rerun.
    return [r for r in latest(read_store(store)) if r['role']=='complex'
            and relaxed(r['calculation'])
            and (r['surface'],r['molecule'],r['calculation']['functional']) in targets]


def prepare(project, campaign, audit, store):
    if campaign.parent != project/'completion_runs' or campaign.exists():
        raise ValueError('Use a new immediate subdirectory of completion_runs')
    rows = json.loads(audit.read_text())
    if not rows or len({r['directory'] for r in rows}) != len(rows):
        raise ValueError('Empty or duplicate source directory list')
    alternatives = completed_alternatives(rows, store)
    if alternatives:
        raise ValueError('Completed alternatives need identity review: '+str([
            r['calculation']['directory'] for r in alternatives]))
    prepared = []
    # Validate the entire campaign before writing any job inputs.
    for row in rows:
        source = Path(row['directory']).resolve()
        if project/'dft_jobs' not in source.parents:
            raise ValueError('Unexpected source tree: '+str(source))
        result = component(source)
        if result['status'] != 'unfinished' or row['status'] != 'unfinished':
            raise ValueError('Source is no longer unfinished: '+str(source))
        if sha(source/'OUTCAR') != row['full_output_sha256']:
            raise ValueError('Source OUTCAR changed since audit: '+str(source))
        if row['normal_completion_anywhere'] or row['ionic_convergence_anywhere']:
            raise ValueError('Completion marker needs manual review: '+str(source))
        ident = identity(source)
        if (ident['surface'],ident['molecule'],result['functional']) != (row['surface'],row['molecule'],row['functional']):
            raise ValueError('Audit/source identity mismatch: '+str(source))
        if (result.get('nsw') or 0) <= 0 or result['settings'].get('IBRION') not in ('1','2','3'):
            raise ValueError('Source is not a relaxation: '+str(source))
        atoms = validate_coordinates(source)
        names = ['POSCAR','CONTCAR','INCAR','POTCAR','KPOINTS']
        if result['functional']=='beef_vdw': names.append('vdw_kernel.bindat')
        hashes = {name:sha(source/name) for name in names}
        original = incar_values((source/'INCAR').read_text())
        settings = restart_settings(original, row['zbrent_fatal'], row['scf_limit_warning'])
        prepared.append((row, source, atoms, hashes, original, settings))
    campaign.mkdir(parents=True)
    jobs = []
    # Put a representative of each restart class first, for early inspection.
    prepared.sort(key=lambda p:(not p[0]['zbrent_fatal'], not p[0]['scf_limit_warning'], str(p[1])))
    representatives = []
    for predicate in [lambda p:p[0]['zbrent_fatal'],
                      lambda p:p[0]['scf_limit_warning'] and p[0]['functional']=='beef_vdw',
                      lambda p:p[0]['scf_limit_warning'] and p[0]['functional']=='r2scan',
                      lambda p:not p[0]['zbrent_fatal'] and not p[0]['scf_limit_warning']]:
        candidate = next((p for p in prepared if predicate(p)), None)
        if candidate is not None:
            prepared.remove(candidate)
            representatives.append(candidate)
    for index, (row, source, atoms, hashes, original, settings) in enumerate(representatives+prepared):
        directory = campaign/'relaxed'/source.parent.name/source.name
        directory.mkdir(parents=True)
        for name in hashes:
            if sha(source/name) != hashes[name]:
                raise ValueError('Source input changed during staging: '+str(source/name))
        shutil.copyfile(source/'CONTCAR', directory/'POSCAR')
        for name in ['POTCAR','KPOINTS','vdw_kernel.bindat']:
            if name in hashes: shutil.copyfile(source/name, directory/name)
        (directory/'INCAR').write_text(''.join(k+' = '+v+'\n' for k,v in settings.items()))
        staged = {p.name:sha(p) for p in directory.iterdir()}
        assert staged['POSCAR']==hashes['CONTCAR']
        assert staged['POTCAR']==hashes['POTCAR']
        changes = {k:dict(before=original.get(k),after=settings.get(k))
                   for k in sorted(set(original)|set(settings)) if original.get(k)!=settings.get(k)}
        job = dict(array_index=index, role='relaxed', system=source.parent.name,
                   surface=row['surface'], molecule=row['molecule'], functional=row['functional'],
                   directory=str(directory), source_directory=str(source), target_systems=[source.parent.name],
                   source_sha256=hashes, source_outcar_sha256=row['full_output_sha256'],
                   staged_sha256=staged, changes=changes, cause=row['cause'],
                   scf_limit_warning=row['scf_limit_warning'], natoms=len(atoms), cell=atoms.cell.tolist(), status='prepared')
        (directory/'completion_inputs.json').write_text(json.dumps(job,indent=2)+'\n')
        jobs.append(job)
    manifest = dict(schema_version=1, created_at=datetime.now(timezone.utc).isoformat(),
        project_root=str(project), audit_sha256=sha(audit), jobs=jobs,
        policy='Restart from saved CONTCAR with identical potentials, cell, constraints, functional and force threshold. '
               'ZBRENT: EDIFF<=1e-6, IBRION=1, POTIM=0.2. SCF warning: ALGO=Normal, NELM>=300. '
               'All inputs fingerprinted; source files retained. No automatic resubmissions.',
        duplicate_check='No completed same-system relaxation in merged Perlmutter/Kestrel component archive. '
                        'Separately reviewed Kestrel relaxation export; live Kestrel authentication unavailable.')
    (campaign/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (campaign/'joblist.txt').write_text(''.join(j['directory']+'\n' for j in jobs))
    (campaign/'logs').mkdir()
    shutil.copyfile(Path(__file__).with_name('completion_array.slurm'), campaign/'run_array.slurm')
    return manifest


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['project','campaign','audit','store']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.project.resolve(), args.campaign.resolve(), args.audit.resolve(), args.store.resolve())
    print('Prepared',len(result['jobs']),'isolated restart jobs in',args.campaign)
