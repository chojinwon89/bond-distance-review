#!/usr/bin/env python3
"""Audit relaxed VASP energies, convergence, composition and POTCAR identities.

Read calculation outputs only. Positive absolute energies are not a failure.
The +/-5 eV adsorption-energy screen is a review flag, not a convergence test.
"""
import argparse
import collections
import csv
import functools
import gzip
import importlib.util
import json
import math
import re
from pathlib import Path

FUNCTIONALS = {'PBE': 'pbe', 'PBE_D3': 'pbe_d3', 'r2scan': 'r2scan', 'beef_vdw': 'beef_vdw'}
ENERGY = re.compile(r'free  energy   TOTEN\s*=\s*([-+\d.Ee]+)')


def parse_output(head, tail):
    titles = re.findall(r'^\s*TITEL\s*=\s*(.+?)\s*$', head, re.M)
    species = re.findall(r'VRHFIN\s*=\s*([A-Z][a-z]?)\s*:', head)
    count = re.search(r'ions per type\s*=\s*([\d ]+)', head)
    counts = [int(x) for x in count.group(1).split()] if count else []
    metadata_ok = bool(counts) and len(species) == len(counts) == len(titles)
    composition = dict(zip(species, counts)) if metadata_ok else {}
    potentials = dict(zip(species, titles)) if metadata_ok else {}
    energies = list(ENERGY.finditer(tail))
    energy = float(energies[-1].group(1)) if energies else None
    nsw = re.search(r'\bNSW\s*=\s*(\d+)', head)
    if 'General timing and accounting' not in tail:
        convergence = 'unfinished'
    elif not nsw:
        convergence = 'unknown_NSW'
    elif energy is None:
        convergence = 'no_final_energy'
    else:
        start = energies[-2].end() if len(energies) > 1 else 0
        electronic = 'aborting loop because EDIFF is reached' in tail[start:energies[-1].start()]
        if not electronic:
            convergence = 'electronic_unconverged'
        elif int(nsw.group(1)) > 0 and 'reached required accuracy' not in tail:
            convergence = 'ionic_unconverged'
        else:
            convergence = 'converged'
    return dict(energy=energy, convergence=convergence, composition=composition,
                potentials=potentials, metadata_ok=metadata_ok)


@functools.lru_cache(None)
def read_output(path):
    actual = path if path.exists() else path.with_name(path.name + '.gz')
    if not actual.exists():
        result = parse_output('', '')
        result.update(convergence='missing', path=str(path))
        return result
    opener = gzip.open if actual.suffix == '.gz' else open
    with opener(actual, 'rb') as f:
        head = f.read(100000).decode(errors='replace')
        # Gzip seek needs the uncompressed size; retain a bounded tail while streaming.
        if actual.suffix == '.gz':
            tail_bytes = head.encode()
            while True:
                chunk = f.read(1000000)
                if not chunk:
                    break
                tail_bytes = (tail_bytes + chunk)[-1000000:]
        else:
            f.seek(max(0, actual.stat().st_size - 1000000))
            tail_bytes = f.read()
    result = parse_output(head, tail_bytes.decode(errors='replace'))
    result['path'] = str(actual)
    return result


def assess(complex_result, slab, molecule, metal):
    parts = [complex_result, slab, molecule]
    notes = []
    values = [p['energy'] for p in parts]
    raw = values[0] - values[1] - values[2] if all(v is not None and math.isfinite(v) for v in values) else None
    if raw is None:
        return 'error', 'Missing or non-finite component energy', raw
    if not all(p['metadata_ok'] for p in parts):
        return 'metadata_missing', 'Cannot verify OUTCAR species, atom counts and POTCAR TITEL identities', raw
    mismatches = []
    for label, ref in [('slab', slab), ('molecule', molecule)]:
        for element, title in ref['potentials'].items():
            job_title = complex_result['potentials'].get(element)
            if job_title != title:
                mismatches.append('{} {}: complex [{}], reference [{}]'.format(label, element, job_title or 'missing', title))
    if mismatches:
        notes.append('; '.join(mismatches))
    comp = complex_result['composition']
    if comp.get(metal, 0) != slab['composition'].get(metal, 0):
        notes.append('slab size mismatch: complex has {} {}, reference has {}'.format(comp.get(metal, 0), metal, slab['composition'].get(metal, 0)))
        status = 'slab_mismatch'
    else:
        expected = collections.Counter(slab['composition']) + collections.Counter(molecule['composition'])
        status = 'composition_mismatch' if dict(expected) != comp else 'ok'
        if status != 'ok':
            notes.append('Complex composition does not equal slab plus gas molecule')
    if mismatches:
        status = 'pseudopotential_mismatch'
    if status != 'ok':
        return status, '; '.join(notes), raw
    if not all(p['convergence'] == 'converged' for p in parts):
        return 'unconverged', 'See component convergence statuses', raw
    if abs(raw) > 5:
        return 'energy_review', '|E_ads| > 5 eV: review required; not evidence of failed convergence', raw
    return 'ok', '', raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vasp-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'dft_perlmutter_energy_audit.csv')
    args = parser.parse_args()
    base = args.vasp_root.resolve()
    spec = importlib.util.spec_from_file_location('collector', base / 'calc_binding_energy.py')
    collector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(collector)
    jobs = collector.discover_system_dirs(base / 'dft_jobs')
    rows = []
    for directory, functional in FUNCTIONALS.items():
        for job in jobs:
            surface, molecule = collector.parse_surface_molecule(job.name, True)
            gas_name = collector.MOLECULE_ALIASES.get(molecule, molecule)
            paths = [job / directory / 'OUTCAR', base / 'vasp_slab' / surface / directory / 'OUTCAR',
                     base / 'vasp_mol' / gas_name / directory / 'OUTCAR']
            parts = [read_output(p) for p in paths]
            status, note, raw = assess(*parts, metal=collector._surface_metal(surface))
            row = dict(functional=functional, system=job.name, surface=surface, molecule=molecule,
                       source_dir='dft_jobs', E_slab_mol=parts[0]['energy'], E_slab=parts[1]['energy'],
                       E_mol=parts[2]['energy'], E_ads=raw if status == 'ok' else None, status=status, note=note,
                       convergence='; '.join('{}: {}'.format(k, p['convergence']) for k, p in zip(['complex', 'slab', 'molecule'], parts)),
                       publishable=str(status == 'ok').lower())
            for label, part in zip(['complex', 'slab', 'molecule'], parts):
                row[label + '_outcar'] = str(Path(part['path']).relative_to(base))
            row['E_ads_raw'] = raw
            for label, part in zip(['complex', 'slab', 'molecule'], parts):
                row[label + '_potentials'] = json.dumps(part['potentials'], sort_keys=True)
                row[label + '_composition'] = json.dumps(part['composition'], sort_keys=True)
            rows.append(row)
        print('Audited', directory, flush=True)
    if not rows:
        raise ValueError('No jobs discovered; refusing to replace audit')
    with args.output.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('Audit status', collections.Counter(r['status'] for r in rows))
    print('Publishable by functional', collections.Counter(r['functional'] for r in rows if r['publishable'] == 'true'))


if __name__ == '__main__':
    main()
