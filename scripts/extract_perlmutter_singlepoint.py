#!/usr/bin/env python3
"""Audit Perlmutter SPE on original ML POSCARs using relaxed references.

Reuses the existing Kestrel completion, settings, composition and energy-screen
rules. Calculation files are read only; matching PAW potential identities are required.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from extract_kestrel_singlepoint import component, calculation_dirs, is_relaxed, assess, write_csv, FUNC_DIRS

from completion_support import completion_jobs, matching_slab
from molecule_names import canonical, ALIASES

ROOT = Path(__file__).resolve().parents[1]
def original_geometry(directory, expected_source=None):
    """Require the staged POSCAR to match the recorded original, never CONTCAR."""
    record = json.loads((directory / 'spe_inputs.json').read_text())
    expected_source = expected_source or directory.parent.parent / directory.name
    if Path(record['source']).resolve() != expected_source or Path(record['destination']).resolve() != directory:
        raise ValueError('Staging record source/destination disagrees with calculation directory')
    if record['functional'] != directory.name:
        raise ValueError('Staging record functional mismatch')
    digest = hashlib.sha256((directory / 'POSCAR').read_bytes()).hexdigest()
    if record['geometry_source'] != 'POSCAR':
        raise ValueError('SPE geometry is not the original POSCAR')
    if digest != record['source_sha256']['POSCAR'] or digest != record['staged_sha256']['POSCAR']:
        raise ValueError('SPE POSCAR differs from recorded original POSCAR')
    for name in ('INCAR', 'KPOINTS'):
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != record['staged_sha256'][name]:
            raise ValueError('Staged input changed: ' + name)
    return {'geometry_source': 'original ML POSCAR', 'poscar_sha256': digest,
            'source': record['source'], 'staged_sha256': record['staged_sha256']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, default=ROOT)
    args = parser.parse_args()
    project, output = args.project_root.resolve(), args.output_root.resolve()
    references = defaultdict(list)
    components = {}
    for role, root in [('slab', project/'vasp_slab'), ('molecule', project/'vasp_mol')]:
        for system in sorted(root.iterdir()):
            if not system.is_dir():
                continue
            label = canonical(system.name) if role == 'molecule' else re.sub(r'_n\d+$', '', system.name)
            for directory, expected in calculation_dirs(system):
                result = component(directory)
                components[str(directory)] = result
                if expected is not None and result['functional'] != expected:
                    continue
                if canonical(result.get('settings', {}).get('SYSTEM', '')) not in (canonical(system.name), label):
                    continue
                if is_relaxed(result, result['functional']):
                    references[(role, label, result['functional'])].append(result)
    completion = completion_jobs(project)
    for job in completion:
        if job['role']=='molecule':
            directory=Path(job['directory']);result=component(directory)
            result['completion_job']=job;components[str(directory)]=result
            if is_relaxed(result,job['functional']):references[('molecule',canonical(job['molecule']),job['functional'])].append(result)
            continue
        if job['role'] != 'slab':
            continue
        directory = Path(job['directory'])
        result = component(directory)
        result['completion_job'] = job
        components[str(directory)] = result
        expected = FUNC_DIRS.get(job['functional'],job['functional'])
        if result.get('settings', {}).get('SYSTEM') == job['surface'] and is_relaxed(result, expected):
            references[('slab', job['surface'], expected)].append(result)
    print('Audited references:', len(components), flush=True)
    rows = []
    for system in sorted((project/'dft_jobs').iterdir()):
        if not system.is_dir():
            continue
        parsed = re.fullmatch(r'(.+)_([A-Z][a-z]?\d+)(?:_(.+))?', system.name)
        if not parsed:
            continue
        raw_molecule, surface, site = parsed.groups()
        molecule = canonical(raw_molecule)
        metal = re.match(r'[A-Z][a-z]?', surface).group()
        candidates = [(dirname, functional, system/'singlepoint'/dirname, False) for dirname, functional in FUNC_DIRS.items()]
        candidates += [(j['functional'], FUNC_DIRS[j['functional']], Path(j['directory']), True)
                       for j in completion if j['role']=='spe' and j['system']==system.name]
        for dirname, functional, directory, retry in candidates:
            if not directory.is_dir():
                continue
            comp = component(directory)
            comp['role'] = 'complex'
            components[str(directory)] = comp
            slabs = [r for r in references[('slab', surface, functional)]
                     if r['composition'] == {metal: comp['composition'].get(metal, -1)}
                     and ('completion_job' not in r or matching_slab(r['completion_job'], system.name, dirname, directory))]
            adsorbate = dict(comp['composition']); adsorbate.pop(metal, None)
            molecules = [r for r in references[('molecule', molecule, functional)] if r['composition'] == adsorbate]
            slabs.sort(key=lambda r: (not np.allclose(np.array(r.get('cell', np.zeros((3,3))))[:2],
                                  np.array(comp.get('cell', np.ones((3,3))))[:2], atol=1e-5, rtol=0),
                                  Path(r['directory']).name != dirname, r['directory']))
            molecules.sort(key=lambda r: (Path(r['directory']).parent.name != raw_molecule,
                                          Path(r['directory']).name != dirname, r['directory']))
            slab, mol = slabs[0] if slabs else None, molecules[0] if molecules else None
            status, note, energy = assess(comp, slab, mol, functional)
            if not (directory/'OUTCAR').exists() and not (directory/'OUTCAR.gz').exists():
                status, note = 'output_missing', 'Prepared SPE input; no OUTCAR available at extraction'
            if comp.get('settings', {}).get('SYSTEM') != system.name:
                status, note = 'system_mismatch', 'INCAR SYSTEM differs from exact system directory'
            geometry = {}
            try:
                geometry = original_geometry(directory, system/dirname)
            except (OSError, ValueError, KeyError) as error:
                status, note = 'geometry_provenance_mismatch', str(error)
            comp['geometry_provenance'] = geometry
            rows.append(dict(surface=surface, molecule=molecule, functional=functional,
                system=system.name, site=site or '', retry=str(retry).lower(), status=status, note=note, E_ads_SPE=energy,
                complex_directory=str(directory), complex_NSW=comp['nsw'],
                slab_directory=slab['directory'] if slab else '', molecule_directory=mol['directory'] if mol else '',
                E_complex=comp['energy'], E_slab_relaxed=slab['energy'] if slab else None,
                E_molecule_relaxed=mol['energy'] if mol else None, slab_NSW=slab['nsw'] if slab else None,
                molecule_NSW=mol['nsw'] if mol else None, complex_status=comp['status'],
                poscar_sha256=geometry.get('poscar_sha256', ''), selected='false'))
    selected = {}
    for row in rows:
        if row['status'] != 'ok':
            continue
        key = row['surface'], row['molecule'], row['functional']
        # Keep separate adsorption-site candidates in the audit; no energy minimization.
        if key not in selected or (row['retry']=='true', bool(row['site']), row['complex_directory']) < (selected[key]['retry']=='true', bool(selected[key]['site']), selected[key]['complex_directory']):
            selected[key] = row
    for row in selected.values():
        row['selected'] = 'true'
    if not rows:
        raise ValueError('No SPE candidates; refusing to replace existing outputs')
    write_csv(output/'dft_perlmutter_singlepoint_audit.csv', rows, list(rows[0]))
    write_csv(output/'dft_perlmutter_singlepoint.csv', [selected[k] for k in sorted(selected)], list(rows[0]))
    metadata = dict(extracted_at=datetime.now(timezone.utc).isoformat(), project_root=str(project),
        formula='E_ads_SPE = E_complex_NSW0 - E_slab_relaxed - E_molecule_relaxed',
        policy='Require matching PAW potential identities as well as completion, electronic convergence, relaxed reference and functional/composition checks. Mismatched energies are diagnostic only. The final cluster policy revalidates published values. No scaling or offsets.',
        geometry='POSCAR hash must match both source and staged POSCAR hashes from spe_inputs.json. Check staged INCAR and KPOINTS hashes.',
        selection='Complex: prefer original over completion retry, then no site suffix, then lexical path; never by energy. References: matching functional and composition, prefer matching slab in-plane cell and original molecule spelling, then functional directory and lexical path.',
        counts=dict(Counter(row['status'] for row in rows)), screened_results=len(selected), components=components)
    (output/'dft_perlmutter_singlepoint_sources.json').write_text(json.dumps(metadata, indent=1)+'\n')
    print('SPE audit:', len(rows), metadata['counts'], 'Selected:', len(selected))


if __name__ == '__main__':
    main()
