#!/usr/bin/env python3
"""Portable, append-preserving VASP component snapshots and missing-energy requests.

Collect on either cluster; merge the resulting gzip JSONL on the website host.
This stores total energies independently of whether adsorption can be evaluated.
It never selects cross-cluster references or modifies published adsorption values.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess

from molecule_names import ALIASES, canonical, equivalent_names

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = 1


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def identity(directory):
    """Resolve explicit directory conventions, preserving molecule/isomer identity."""
    parts = Path(directory).parts
    if 'completion_runs' in parts:
        index=parts.index('completion_runs')
        if len(parts)>index+3 and parts[index+2]=='molecule':
            return dict(role='molecule',surface='',molecule=canonical(parts[index+3]),source_name=parts[index+3])
        if len(parts)>index+3 and parts[index+2] in ('spe','relaxed'):
            label=parts[index+3]
            match=re.fullmatch(r'(.+)_([A-Z][a-z]?\d+)(?:_.+)?',label)
            if match:
                return dict(role='complex',surface=match[2],molecule=canonical(match[1]),source_name=label)
        if len(parts)>index+3 and parts[index+2]=='slab':
            label=parts[index+3];match=re.search(r'([A-Z][a-z]?\d+)',label)
            return dict(role='slab',surface=match[1] if match else label,molecule='',source_name=label)
    for tree, role in [('vasp_mol', 'molecule'), ('vasp_slab_kestrel', 'slab'), ('vasp_slab', 'slab')]:
        if tree in parts:
            label = parts[parts.index(tree) + 1]
            return dict(role=role, surface=re.sub(r'_n\d+$', '', label) if role == 'slab' else '',
                        molecule=canonical(label) if role == 'molecule' else '', source_name=label)
    if 'dft_jobs' in parts:
        label = parts[parts.index('dft_jobs') + 1]
        match = re.fullmatch(r'(.+)_([A-Z][a-z]?\d+)(?:_.+)?', label)
        if match:
            return dict(role='complex', surface=match[2], molecule=canonical(match[1]), source_name=label)
    if 'best' in parts:
        rest = list(parts[parts.index('best') + 1:])
        if rest and re.fullmatch(r'C\d+', rest[0]):
            rest.pop(0)
        if rest:
            match = re.fullmatch(r'([A-Z][a-z]?\d+)_(.+)', rest[0])
            if match:
                return dict(role='complex', surface=match[1], molecule=canonical(match[2]), source_name=rest[0])
    return dict(role='unclassified', surface='', molecule='', source_name=Path(directory).name)


def record(cluster, calculation, provenance):
    if not cluster or not Path(calculation['directory']).is_absolute():
        raise ValueError('Cluster and absolute source directory are required')
    result = dict(schema_version=SCHEMA, cluster=cluster, **identity(calculation['directory']),
                  units='eV', energy_kind='VASP free energy TOTEN', calculation=calculation)
    # Provenance import dates do not duplicate identical calculation snapshots.
    result['snapshot_id'] = digest(result)
    result['provenance'] = provenance
    return result


def read_store(path):
    if not path.exists():
        return []
    with gzip.open(path, 'rt') as source:
        records = [json.loads(line) for line in source if line.strip()]
    for row in records:
        validate(row)
    return records


def validate(row):
    if row['schema_version'] != SCHEMA:
        raise ValueError('Unsupported component schema')
    core = {k: v for k, v in row.items() if k not in ('snapshot_id', 'provenance')}
    if digest(core) != row['snapshot_id']:
        raise ValueError('Component snapshot checksum mismatch')


def merge(store, incoming):
    records = {r['snapshot_id']: r for r in read_store(store)}
    for row in incoming:
        validate(row)
        records.setdefault(row['snapshot_id'], row)
    store.parent.mkdir(parents=True, exist_ok=True)
    temporary = store.with_name(store.name + '.tmp')
    # Fixed gzip metadata and deterministic ordering make repeat imports idempotent.
    with temporary.open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as output:
        for key in sorted(records):
            output.write((json.dumps(records[key], sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode())
    temporary.replace(store)
    return len(records)


def latest(records):
    chosen = {}
    for row in records:
        c = row['calculation']; key = (row['cluster'], c['directory'])
        # Observations of running/replaced outputs remain in the store. A genuinely
        # newer failure supersedes an older success in the *current* projection.
        rank = (c.get('outcar_mtime_ns', 0), row['provenance'].get('observed_at', ''), row['snapshot_id'])
        if key not in chosen or rank > chosen[key][0]:
            chosen[key] = rank, row
    return [chosen[key][1] for key in sorted(chosen)]


def collect(cluster, project):
    from extract_kestrel_singlepoint import component
    project = project.resolve(); records = []; seen = set()
    observed = datetime.now(timezone.utc).isoformat()
    for name in ['vasp_mol', 'vasp_slab', 'vasp_slab_kestrel', 'dft_jobs', 'poscar/best']:
        root = project / name
        if not root.is_dir():
            continue
        result = subprocess.run(['rg', '--files', '--hidden', '--no-ignore', str(root), '-g', 'INCAR'],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode not in (0, 1):
            raise RuntimeError(result.stderr)
        paths = sorted(result.stdout.splitlines())
        for path in paths:
            directory = Path(path).parent
            if root not in directory.resolve().parents or str(directory) in seen:
                continue
            seen.add(str(directory))
            records.append(record(cluster, component(directory), dict(kind='live-files', observed_at=observed,
                project_root=str(project), fingerprint='SHA256 of first 100000 bytes + NUL + last 1000000 bytes; not a full-file hash')))
        print('Collected', name, len(paths), flush=True)
    from completion_support import completion_jobs
    for job in completion_jobs(project):
        directory=Path(job['directory'])
        if job['role'] in ('molecule','slab','spe','relaxed'):
            records.append(record(cluster,component(directory),dict(kind='live-files',observed_at=observed,project_root=str(project),completion_manifest=str(directory/'completion_inputs.json'))))
    if not records:
        raise ValueError('No calculations found; store preserved')
    return records


def import_sources(cluster, source):
    data = json.loads(source.read_text())
    provenance = dict(kind='archived-source-export', source_file=source.name,
                      source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      observed_at=data.get('extracted_at', ''),
                      freshness='Saved output metadata; not a fresh cluster read')
    return [record(cluster, c, provenance) for c in data['components'].values()]


def cell_matches(left, right):
    return bool(left and right and len(left) == len(right) == 3
                and all(len(a) == len(b) == 3 for a, b in zip(left, right))
                and all(abs(a-b) <= 1e-5 for x, y in zip(left, right) for a, b in zip(x, y)))


def relaxed(c):
    return c['status'] == 'converged' and (c.get('nsw') or 0) > 0 and c.get('settings', {}).get('IBRION') in ('1', '2', '3')


def reference_candidates(complex_row, components, allow_unverified=False):
    from potential_matching import potential_match
    c = complex_row['calculation']; metal = re.match(r'[A-Z][a-z]?', complex_row['surface']).group()
    slab_composition = {metal: c['composition'].get(metal, 0)}
    molecule_composition = {k: v for k, v in c['composition'].items() if k != metal}
    candidates = {'slab': [], 'molecule': []}
    for row in components:
        r = row['calculation']; role = row['role']
        if role not in candidates or r['functional'] != c['functional'] or not relaxed(r):
            continue
        if role == 'slab':
            matches = row['surface'] == complex_row['surface'] and r['composition'] == slab_composition and (allow_unverified or cell_matches(c.get('cell'), r.get('cell')))
        else:
            matches = row['molecule'] == complex_row['molecule'] and r['composition'] == molecule_composition
        if matches and (allow_unverified or potential_match(c,r)):
            candidates[role].append(row['snapshot_id'])
    return dict(complex_snapshot_id=complex_row['snapshot_id'], required_slab_composition=slab_composition,
                required_molecule_composition=molecule_composition, required_cell_A=c.get('cell'),
                candidate_reference_ids=candidates,
                unresolved_roles=[role for role, ids in candidates.items() if not ids],
                required_potentials=c.get('potentials',{}),
                policy='Require matching PAW TITEL identities including variant/version; verify geometry, constraints and full settings before deriving energy')


def report(store, site):
    records = read_store(store); current = latest(records)
    fields = ['cluster', 'role', 'surface', 'molecule', 'functional', 'mode', 'status', 'E_TOTEN_eV',
              'composition', 'directory', 'snapshot_id', 'source_kind', 'note']
    with (site / 'dft_component_energies.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=fields); writer.writeheader()
        for row in current:
            c = row['calculation']; nsw = c.get('nsw')
            writer.writerow(dict(cluster=row['cluster'], role=row['role'], surface=row['surface'], molecule=row['molecule'],
                functional=c['functional'], mode='SPE' if nsw == 0 else 'relaxed' if nsw and nsw > 0 else 'unknown',
                status=c['status'], E_TOTEN_eV=c.get('energy'), composition=json.dumps(c['composition'], sort_keys=True),
                directory=c['directory'], snapshot_id=row['snapshot_id'], source_kind=row['provenance']['kind'], note=c.get('note', '')))
    coverage = list(csv.DictReader((site / 'dft_completion_coverage.csv').open()))
    relaxed_audit = list(csv.DictReader((site / 'dft_perlmutter_energy_audit.csv').open()))
    completed_relaxed_keys = {(r['surface'], canonical(r['molecule']), r['functional']) for r in relaxed_audit
                              if 'complex: converged;' in r['convergence']}
    index = {}
    for row in current:
        if row['role'] == 'complex':
            c = row['calculation']
            mode = 'SPE' if c.get('nsw') == 0 else 'relaxed' if (c.get('nsw') or 0) > 0 else 'unknown'
            index.setdefault((row['surface'], row['molecule'], c['functional'], mode), []).append(row)
    requests = []
    for entry in coverage:
        for mode in ['relaxed', 'SPE']:
            if entry['E_ads_' + mode]:
                continue
            molecule = canonical(entry['molecule'])
            candidates = index.get((entry['surface'], molecule, entry['functional'], mode), [])
            converged = [r for r in candidates if r['calculation']['status'] == 'converged']
            requests.append(dict(surface=entry['surface'], molecule=molecule, aliases=sorted(set([molecule] + list(equivalent_names(molecule, ALIASES)))),
                functional=entry['functional'], mode=mode, website_snapshot_status=entry[mode+'_status'],
                complex_candidate_ids=[r['snapshot_id'] for r in candidates],
                completed_complexes=[reference_candidates(r, current) for r in converged]))
    (site / 'dft_kestrel_search_requests.json').write_text(json.dumps(dict(schema_version=SCHEMA,
        instructions='Search both fully relaxed and NSW=0 jobs using all explicit aliases. Match the listed reference compositions/cells. These are unresolved website cells, not instructions to submit duplicate jobs. Import results with component_store.py; audit derived energies before updating the page.',
        requests=requests), separators=(',', ':')) + '\n')
    dependencies={}
    for request in requests:
        cell=(request['surface'],request['molecule'],request['functional'],request['mode'])
        for comp in request['completed_complexes']:
            for role in comp['unresolved_roles']:
                name=request['molecule'] if role=='molecule' else request['surface']
                composition=comp['required_molecule_composition'] if role=='molecule' else comp['required_slab_composition']
                key=(role,name,request['functional'],json.dumps(composition,sort_keys=True),
                     json.dumps(comp['required_cell_A']) if role=='slab' else '',
                     json.dumps({el:comp['required_potentials'].get(el) for el in composition},sort_keys=True))
                dependencies.setdefault(key,set()).add(cell)
    with (site/'dft_missing_reference_priorities.csv').open('w',newline='') as out:
        w=csv.writer(out);w.writerow(['role','name','functional','required_composition','required_slab_cell_A','required_potentials','affected_missing_cells','cells'])
        for key,cells in sorted(dependencies.items(),key=lambda x:(-len(x[1]),x[0])):
            w.writerow([*key,len(cells),json.dumps(sorted(cells))])
    summary = dict(schema_version=SCHEMA, component_snapshots=len(records), current_calculations=len(current),
        by_cluster=dict(Counter(r['cluster'] for r in current)),
        converged_by_cluster=dict(Counter(r['cluster'] for r in current if r['calculation']['status']=='converged')),
        missing_cells=len(requests), missing_cells_with_saved_converged_complex=sum(bool(r['completed_complexes']) for r in requests),
        relaxed_blank_reasons=dict(Counter(r['relaxed_status'] for r in coverage if not r['E_ads_relaxed'])),
        SPE_blank_reasons=dict(Counter(r['SPE_status'] for r in coverage if not r['E_ads_SPE'])),
        unconverged_relaxed_cells_with_converged_complex=sum(not r['E_ads_relaxed'] and r['relaxed_status']=='unconverged'
            and (r['surface'], canonical(r['molecule']), r['functional']) in completed_relaxed_keys for r in coverage),
        caveat='Website gap reasons use the dated completion snapshot. Live Perlmutter components and archived Kestrel components retain separate provenance. Completed candidate jobs may differ in geometry from the website target.')
    (site / 'dft_component_store_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    write_explanation(site, current, coverage, summary)
    print(json.dumps(summary, indent=2))


def write_explanation(site, current, coverage, summary):
    snapshot = json.loads((site / 'dft_completion_summary.json').read_text())['snapshot']
    escape = html.escape
    recovery_file=site/'dft_gas_reference_recovery.json'
    recovered=[]
    if recovery_file.exists():
        for job in json.loads(recovery_file.read_text())['jobs']:
            energy=job.get('energy_TOTEN_eV')
            value=f'; gas TOTEN {energy:.8f} eV' if job['calculation_status']=='converged' and energy is not None else ''
            recovered.append(escape(job['functional']+': '+job['scheduler_state']+', '+job['calculation_status']+value))
    recovery_status='<p><strong>Current recovery snapshot:</strong> '+'; '.join(recovered)+'. Converged retries now supply the matching table references; the failed originals below are retained for diagnosis.</p>' if recovered else ''
    example_paths = [
        ('perlmutter', '/dft_jobs/CH3OCH3_Ag111/PBE', 'DME + Ag111, relaxed'),
        ('perlmutter', '/dft_jobs/CH3OCH3_Ag111/singlepoint/PBE', 'DME + Ag111, SPE'),
        ('perlmutter', '/vasp_slab/Ag111/PBE', 'Standard clean Ag111 slab'),
        ('perlmutter', '/vasp_mol/CH3OCH3/PBE', 'DME gas reference (CH3OCH3)'),
        ('perlmutter', '/vasp_mol/DME/PBE', 'DME gas reference (DME alias)'),
        ('kestrel', '/vasp_slab/Ag111_n36/PBE', 'Saved Kestrel size-specific slab'),
    ]
    example_rows = []
    for cluster, suffix, label in example_paths:
        found = next((r for r in current if r['cluster'] == cluster and r['calculation']['directory'].endswith(suffix)), None)
        if not found:
            continue
        c = found['calculation']; energy = 'unavailable' if c.get('energy') is None else f"{c['energy']:.8f}"
        example_rows.append(f'<tr><td>{escape(label)}<br><small>{cluster}</small></td><td>{energy}</td>'
            f'<td>{escape(json.dumps(c["composition"], sort_keys=True))}</td><td>{escape(c["status"])}</td>'
            f'<td><code>{escape(c["directory"])}</code><br>{escape(c.get("note", ""))}</td></tr>')
    reason_labels = {
        'potential mismatch': 'The former component subtraction used different PAW potentials; a fully matched replacement is not yet available.',
        'refs unverified': 'The historical number has no fully validated matching-potential component set in the current store.',
        'unconverged': 'At least one component is unfinished or unconverged; this can be the gas reference, not the adsorbed complex.',
        'slab_mismatch': 'The audited slab atom count differs from the complex. Other reference failures may coexist.',
        'no matching Perlmutter calculation': 'No candidate matched the current Perlmutter audit layout and system key. Check Kestrel and other saved layouts.',
        'error': 'The reference/output is missing, unreadable or fails a metadata check; inspect the audit note.',
        'reference_unavailable': 'The SPE complex passed its completion checks, but a suitable converged slab or molecule reference was unavailable.',
        'complex_unconverged': 'The selected complex output did not pass completion/electronic-convergence checks.',
        'SPE pending': 'The separate NSW=0 calculation was queued at this snapshot. A finished relaxation does not finish its SPE job.',
        'SPE running': 'The separate NSW=0 calculation was running at this snapshot.',
        'output_missing': 'No output was available in the saved energy audit. The scheduler snapshot may be later; rerun extraction to check recently completed jobs.',
        'no matching Perlmutter calculation; check Kestrel sources': 'No matching Perlmutter SPE candidate; search Kestrel sources using the shared alias list.',
    }
    reason_rows = ''.join(f'<tr><td>{mode}</td><td>{escape(reason)}</td><td>{count}</td><td>{escape(reason_labels.get(reason, reason))}</td></tr>'
        for mode, counts in [('Relaxed', summary['relaxed_blank_reasons']), ('SPE', summary['SPE_blank_reasons'])]
        for reason, count in sorted(counts.items(), key=lambda x: -x[1]))
    functional_rows = ''.join(f'<tr><td>{func}</td><td>{sum(not r["E_ads_relaxed"] for r in coverage if r["functional"]==func)}</td>'
        f'<td>{sum(not r["E_ads_SPE"] for r in coverage if r["functional"]==func)}</td></tr>'
        for func in ['pbe', 'pbe_d3', 'r2scan', 'beef_vdw'])
    page = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Why completed DFT jobs can leave binding-energy gaps</title>
<style>body{{margin:0;background:#0f1116;color:#e6e6e6;font:16px/1.6 system-ui,sans-serif}}main{{max-width:1100px;margin:auto;padding:30px 22px 70px}}a{{color:#8abaff}}h1{{font-size:29px}}h2{{margin-top:32px;font-size:22px}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{border:1px solid #444;padding:10px;text-align:left}}th{{background:#20242e}}code{{overflow-wrap:anywhere}}.scroll{{overflow:auto}}.note{{padding:15px;background:#20242e;border-left:3px solid #e0a800}}small{{color:#abb4c5}}</style></head>
<body><main><a href="dft_comparison.html">← DFT comparison</a>
<h1>Completed calculations and missing binding energies</h1>
<div class="note"><strong>Matching-potential policy:</strong> the latest audit found 1,120 historical candidate subtractions
mixing metal PAW potentials, across Ag, Cu, Pd, Pt and Rh surfaces. A converged job alone cannot validate that subtraction.
The active table now requires a validated matching-potential component set. <code>potential mismatch</code> identifies known mixed-potential references;
<code>refs unverified</code> means the historical number cannot currently be validated, not that a mismatch has been proved.
51 deduplicated clean-slab calculations were submitted as array <code>58411454</code>, at most four concurrent, using the exact metal POTCAR block from the complexes.
<a href="dft_potential_mismatches.csv">Historical mismatch audit</a> · <a href="dft_potential_slab_recovery.json">Recovery inputs and submission</a>.</div>
<p>A blank binding-energy cell does <strong>not</strong> mean the adsorbed-complex calculation failed.
The page reports adsorption energies, while an individual VASP job reports a total energy.
Three suitable component results are needed:</p>
<p><code>E_ads = E(complex) − E(clean slab, relaxed) − E(gas molecule, relaxed)</code></p>
<p>For the SPE column, the complex must have NSW=0. The current convention uses relaxed slab and gas references for both columns.
A converged relaxation and a converged SPE are separate results, and neither supplies a missing gas or clean-slab calculation.</p>
<h2>Ag111 · DME (CH3OCH3), PBE</h2>
{recovery_status}
<p><strong>Both complex jobs are converged.</strong> Their total energies are available below and in the shared component download.</p>
<div class="scroll"><table><thead><tr><th>Component</th><th>TOTEN (eV)</th><th>Atoms</th><th>Output status</th><th>Source and diagnostic</th></tr></thead>
<tbody>{''.join(example_rows)}</tbody></table></div>
<p>The Perlmutter complex contains 36 Ag atoms, whereas the standard clean slab contains 64 and uses a different supercell.
Both spellings of the PBE gas reference end with <code>ZBRENT: fatal error in bracketing</code>, without successful ionic completion.
These are two independent blockers. The last gas TOTEN is retained as a diagnostic, not a converged reference.</p>
<p>Blindly subtracting those Perlmutter values produces <code>+65.17064489 eV</code>. That number removes 64 Ag atoms from a 36-Ag complex and uses an unfinished molecule.
It is stored in the raw audit, but is not a valid binding energy. Multiplying the slab energy by 36/64 would not repair the inconsistent surface reference.</p>
<p>The archived Kestrel <code>Ag111_n36/PBE</code> result has 36 Ag atoms and the same full cell as the complex.
It is a useful reference candidate; the archived Kestrel PBE DME gas outputs are also unfinished.
Fresh Kestrel files must be checked for a completed replacement. Every cross-cluster reference choice is recorded in the source-selection download.</p>
<p>Two isolated Perlmutter DME gas-reference retries (array <code>58376018</code>) now target PBE and PBE+D3.
Their <a href="dft_gas_reference_recovery.json">dated job-status snapshot</a> distinguishes scheduler state from verified convergence.
The completed BEEF-vdW gas geometry supplies starting coordinates only; each retry computes its own functional-specific energy.
The first retry array exhausted its electronic iteration limit and was stopped with outputs preserved; the replacement uses the Davidson electronic solver.</p>
<p>For BEEF-vdW, the saved Perlmutter and Kestrel complex SPE totals differ by about 0.000001 eV.
Switching that complex does not resolve the unusual binding energy: the current matched references give approximately
−25.227 eV relaxed and −25.218 eV SPE. The potential audit subsequently identified Ag/Ag_pv mixing in this subtraction; these historical numbers are now withheld from the active table until a matching-potential reference is available.
The <a href="dft_cluster_duplicates.csv">duplicate comparison</a> records both energies and the available structure evidence.</p>
<p>DME and CH3OCH3 already resolve to the same molecule. Ethanol is kept distinct despite the shared gross formula C2H6O.
Potential variants remain recorded; the current policy requires matching PAW identities for complex and references.
Converged, composition-matched unusual binding energies are shown with an amber review marker; magnitude alone does not cause these blanks.</p>
<h2>What the blanks mean across the page</h2>
<p>Website/scheduler snapshot: <code>{escape(snapshot)}</code>. There are 415 systems × 4 functionals = 1,660 possible entries in each energy column.
{sum(not r['E_ads_relaxed'] for r in coverage)} relaxed binding energies and {sum(not r['E_ads_SPE'] for r in coverage)} SPE binding energies are missing.
Counts below are <strong>table cells</strong>, not individual job directories; aliases, adsorption sites and retry layouts can generate multiple job candidates.</p>
<div class="scroll"><table><thead><tr><th>Column</th><th>Primary audit status</th><th>Missing cells</th><th>Meaning</th></tr></thead><tbody>{reason_rows}</tbody></table></div>
<p>Statuses report the first or grouped audit blocker, not every possible failure. For example, the Ag111–DME relaxed row is classified as slab mismatch even though its gas reference also failed.
An “unconverged” relaxed row can therefore have a converged complex: this is true for {summary['unconverged_relaxed_cells_with_converged_complex']} of the {summary['relaxed_blank_reasons'].get('unconverged',0)} rows with that label.
Historical published energies are retained if a newer audit cannot replace them.</p>
<table><thead><tr><th>Functional</th><th>Relaxed blanks</th><th>SPE blanks</th></tr></thead><tbody>{functional_rows}</tbody></table>
<p>The shared store contains a converged complex candidate for {summary['missing_cells_with_saved_converged_complex']} of the {summary['missing_cells']} missing cells.
Candidates can have different geometry or site provenance, so this is a search lead, not a count of immediately publishable binding energies.</p>
<h2>Reusable Perlmutter and Kestrel data</h2>
<p><strong>Size-specific Kestrel slabs are available.</strong> The saved vasp_slab records include Ag111_n36,
Au111_n36, Pd111_n36 and Pt111_n36 under all four functionals, along with other slab sizes.
The active Kestrel location supplied by the user is <code>/scratch/jcho5/goad-global-optimization/vasp_slab</code>.
Archived records retain their original <code>/kfs3/scratch/...</code> paths.</p>
<p>The shared-reference derivation now fills missing table cells using these converged, composition- and cell-matched slabs,
with matching functional, cutoff and stored core settings. The subsequent source-selection pass prefers usable Perlmutter candidates;
an unusual result can use a Kestrel replacement only with verified structure equivalence. Historical values lacking a validated matching-potential component set are withheld.
Unusual values remain visible with review markers. See the
<a href="dft_shared_reference_energies.csv" download>derived energies with all three component totals and source paths</a>
and <a href="dft_shared_reference_audit.json">reference-matching policy and audit</a>.
This uses saved Kestrel records; it is not a fresh read of the Kestrel filesystem.</p>
<p>The store currently represents {summary['current_calculations']} distinct cluster/path calculations:
{summary['by_cluster'].get('perlmutter',0)} Perlmutter and {summary['by_cluster'].get('kestrel',0)} archived Kestrel calculations.
It keeps total energies even when adsorption is unavailable, together with cluster, absolute source path, canonical molecule,
functional, NSW, composition, cell, potentials, convergence, settings and available output fingerprints.
Repeated imports preserve distinct output snapshots and do not erase the other cluster.</p>
<p>Fresh Perlmutter reads and archived Kestrel exports are explicitly distinguished. Archived results are not represented as a live Kestrel search.
The portable collector searches INCAR files under vasp_mol, vasp_slab, vasp_slab_kestrel, dft_jobs and poscar/best,
including nested fully_relaxed and singlepoint layouts. NSW determines calculation type.</p>
<ul><li><a href="dft_component_energies.csv" download>All component total energies and status (CSV)</a></li>
<li><a href="dft_component_store.jsonl.gz" download>Reusable component snapshots with full metadata (gzip JSONL)</a></li>
<li><a href="dft_kestrel_search_requests.json" download>Missing-cell Kestrel search requests, aliases, required atoms/cells and candidate IDs</a></li>
<li><a href="dft_missing_reference_priorities.csv" download>Missing references grouped by the cells they can complete</a></li>
<li><a href="dft_cluster_selection.csv" download>Perlmutter-first selection decisions</a> · <a href="dft_cluster_duplicates.csv" download>Cross-cluster duplicate comparisons</a></li>
<li><a href="dft_component_store_summary.json">Store summary</a> · <a href="dft_completion_coverage.csv" download>Every page cell and its audit status</a></li>
<li><a href="https://github.com/chojinwon89/bond-distance-review/blob/main/scripts/COMPONENT_STORE.md">Portable collection and merge instructions</a></li></ul>
<p>Missing binding energies also leave dependent ML differences and comparison points empty. Review-marked energies remain in tables but are excluded from paired figure statistics.
The 15 missing DFT contact measurements are a separate geometry-source issue requiring their Kestrel CONTCAR files.</p>
</main></body></html>'''
    (site / 'dft_data_gaps.html').write_text(page)
    comparison = (site / 'dft_comparison.html').read_text()
    comparison = re.sub(r'<!-- component-availability -->.*?<!-- /component-availability -->\n?', '', comparison, flags=re.S)
    note = '''<!-- component-availability -->
<div class="note info"><b>A blank binding energy can still have a converged complex.</b>
Binding energies require compatible complex, clean-slab and gas-molecule results. Ag111–DME PBE has completed relaxed and SPE complexes; its original gas reference failed and its standard Perlmutter slab has the wrong atom count. Converged recovery references are applied when available.
Size-specific Kestrel slabs are used only when their PAW potential identities also match the complex. The potential audit supersedes the earlier variant-permissive policy.
<br><a href="dft_data_gaps.html">Detailed explanation, Ag111–DME component energies, and counts of every gap reason</a> &middot;
<a href="dft_component_energies.csv" download>Component total energies (Perlmutter + archived Kestrel)</a> &middot;
<a href="dft_kestrel_search_requests.json" download>Kestrel search requests</a>.</div>
<!-- /component-availability -->
'''
    comparison = comparison.replace('<h2>Per-system structure', note + '<h2>Per-system structure')
    (site / 'dft_comparison.html').write_text(comparison)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--store', type=Path, default=ROOT / 'dft_component_store.jsonl.gz')
    commands = p.add_subparsers(dest='command', required=True)
    collect_parser = commands.add_parser('collect')
    collect_parser.add_argument('--cluster', required=True)
    collect_parser.add_argument('--project-root', type=Path, required=True)
    import_parser = commands.add_parser('import-sources')
    import_parser.add_argument('--cluster', required=True)
    import_parser.add_argument('--source', type=Path, required=True)
    merge_parser = commands.add_parser('merge')
    merge_parser.add_argument('--input', type=Path, required=True)
    report_parser = commands.add_parser('report')
    report_parser.add_argument('--site-root', type=Path, default=ROOT)
    args = p.parse_args()
    if args.command == 'report':
        report(args.store, args.site_root)
    else:
        if args.command == 'collect':
            incoming = collect(args.cluster, args.project_root)
        elif args.command == 'import-sources':
            incoming = import_sources(args.cluster, args.source)
        else:
            incoming = read_store(args.input)
            if not incoming:
                raise ValueError('No imported observations; store preserved')
        print('Stored snapshots:', merge(args.store, incoming))


if __name__ == '__main__':
    main()
