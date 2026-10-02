"""Publish an explicit, dated allowlist of research artifacts; never copy POTCARs.

Usage: python scripts/publish_research_snapshot_20261002.py SOURCE_PROJECT BEEF_STATUS_JSON
Refresh the source CPU report before running this script. Scheduler status is a
timestamped input, not inferred from the existence of prepared calculation files.
"""
from pathlib import Path
import csv
import datetime
import hashlib
import html
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(sys.argv[1]).resolve()
OUT = ROOT / 'research_updates/2026-10-02'
OUT.mkdir(parents=True, exist_ok=True)
inventory = []


def copy(source, destination):
    source = Path(source)
    target = OUT / destination
    assert source.is_file(), source
    assert 'POTCAR' not in source.name.upper()
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    inventory.append(dict(path=str(target.relative_to(OUT)), source=str(source),
                          bytes=target.stat().st_size,
                          sha256=hashlib.sha256(target.read_bytes()).hexdigest()))


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


spe = SOURCE / 'figures/website_spe_parity_20261002'
for name in ['mlip_vs_dft_spe_zoom.png', 'mlip_vs_dft_spe_zoom.pdf',
             'mlip_vs_dft_spe_full_range.png', 'mlip_vs_dft_spe_full_range.pdf',
             'mlip_vs_dft_spe_interactive.html', 'plotted_pairs.csv',
             'missing_ml_pairs.csv', 'validation.json', 'README.txt', 'plot_spe_comparison.py']:
    copy(spe / name, 'spe/' + name)
for name in ['dft_comparison.html', 'dft_comparison_singlepoint.csv',
             'dft_singlepoint_vs_sevennet.csv', 'download_manifest.json']:
    copy(spe / 'sources' / name, 'spe/sources/' + name)

pym = SOURCE / 'benchmarks/pymatgen_best_20261002'
for name in ['best_structures.csv', 'all_candidates.csv', 'selection.json']:
    copy(pym / name, 'pymatgen/' + name)
best = rows(pym / 'best_structures.csv')
for row in best:
    for ext in ['cif', 'vasp']:
        name = row['name'] + '.' + ext
        copy(pym / name, 'pymatgen/' + name)

cpu = SOURCE / 'benchmarks/ag_glycerol_cpu_contact_20261002'
for name in ['status.json', 'all_results.csv', 'paired_comparison.csv', 'comparison.png',
             'comparison.pdf', 'PROTOCOL.txt', 'cpu_references.json',
             'initialization_checks.json', 'preflight.json', 'best_endpoints.json']:
    copy(cpu / name, 'ag_glycerol_cpu/' + name)
for name in ['run_ag_glycerol_cpu_contact.py', 'report_ag_glycerol_cpu_contact.py',
             'check_ag_glycerol_initialization.py']:
    copy(SOURCE / 'benchmarks' / name, 'ag_glycerol_cpu/source/' + name)
for row in rows(cpu / 'all_results.csv'):
    run = cpu / 'runs' / (row['mode'] + '_seed' + row['seed'])
    for name in ['summary.json', 'final.extxyz']:
        copy(run / name, 'ag_glycerol_cpu/runs/' + run.name + '/' + name)

kestrel = SOURCE / 'benchmarks/kestrel_pymatgen_20261002'
for name in ['PROTOCOL.txt', 'manifest.json', 'case_index.csv', 'bundle_audit.json']:
    copy(kestrel / name, 'kestrel/' + name)
beef = SOURCE / 'benchmarks/pymatgen_beef_20261002'
for name in ['PROTOCOL.md', 'job_index.csv', 'submission.json', 'preflight_checks.json']:
    copy(beef / name, 'beef/' + name)
copy(Path(sys.argv[2]).resolve(), 'beef/queue_snapshot.json')

status = json.loads((cpu / 'status.json').read_text())
pairs = rows(cpu / 'paired_comparison.csv')
valid_pairs = [r for r in pairs if r['both_valid'] == 'True']
improved = sum(float(r['contact25_minus_baseline_eV']) < 0 for r in valid_pairs)
queue = json.loads((OUT / 'beef/queue_snapshot.json').read_text())
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
local = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=-6))).strftime('%B %d, %Y, %I:%M %p MDT')
cpu_label = 'Complete' if status['finished'] else 'In progress'
cpu_best = json.loads((cpu / 'best_endpoints.json').read_text())
endpoint_links = ' · '.join('<a href="ag_glycerol_cpu/runs/{}_seed{}/final.extxyz">{}: {:.4f} eV (seed {})</a>'.format(
    mode, record['seed'], 'Original GOAD' if mode == 'baseline' else 'Initial O–Ag = 2.5 Å',
    record['final_E_ads_eV'], record['seed']) for mode, record in cpu_best.items() if record)
pym_rows = ''.join('<tr><td>{}</td><td>{:.4f}</td><td>{:.3f}</td><td>{}/{}</td><td><a href="pymatgen/{}.cif">CIF</a> · <a href="pymatgen/{}.vasp">VASP</a></td></tr>'.format(
    html.escape(r['case'].replace('_', ' · ')), float(r['E_ads_eV']), float(r['nearest_O_metal_A']),
    r['converged_starts'], r['total_starts'], r['name'], r['name']) for r in best)
queue_text = ', '.join('{} {}'.format(n, s.lower()) for s, n in queue['counts'].items())
page = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Research update · 2 October 2026 | Bond distance review</title>
<style>
:root{color-scheme:dark;--bg:#0f1116;--card:#191c24;--fg:#e6e6e6;--muted:#adb3bf;--line:#343a48}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.6 system-ui,sans-serif}
main{max-width:1150px;margin:auto;padding:30px 20px 70px}a{color:#8fbcff}h1{font-size:32px;line-height:1.2}h2{font-size:24px}h3{font-size:18px}
nav{display:flex;gap:18px;flex-wrap:wrap}.muted,figcaption{color:var(--muted);font-size:14px}
section{margin-top:30px;padding:22px;background:var(--card);border:1px solid var(--line);border-radius:12px}
.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin-top:24px}.summary a{display:block;padding:18px;border:1px solid var(--line);border-radius:10px;text-decoration:none}.summary strong{display:block;font-size:25px;color:var(--fg)}
figure{margin:20px 0}img{display:block;width:100%;height:auto;background:white;border-radius:6px}figcaption{margin-top:10px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:11px 10px;border-bottom:1px solid var(--line)}th{color:#afcfff}td{white-space:nowrap}.badge{color:#ffdb9b;font-size:14px}code{overflow-wrap:anywhere}details{margin:15px 0}summary{cursor:pointer;font-weight:600}
</style></head><body><main>
<nav><a href="../../index.html">Structure gallery</a><a href="../../dft_comparison.html">DFT comparison</a><a href="https://github.com/chojinwon89/bond-distance-review">Repository</a></nav>
<h1>Adsorption search and DFT validation</h1>
<p>Research update · 2 October 2026</p>
<p class="muted">Published snapshot: __LOCAL__. Calculation and queue records retain their own UTC check times. This page does not poll the clusters.</p>
<div class="summary"><a href="#spe"><strong>1,450 pairs</strong>Four-functional SPE comparison</a><a href="#pymatgen"><strong>5 systems</strong>Independent pymatgen DFT starts</a><a href="#cpu"><strong>__DONE__/20 searches</strong>Ag–O initialization experiment</a><a href="#kestrel"><strong>50 systems</strong>Kestrel inputs prepared</a></div>
<section id="spe"><h2>GOAD + SevenNet-Omni versus DFT single points</h2>
<p>The horizontal axis is the DFT single-point adsorption energy at the MLIP geometry; the vertical axis is the MLIP adsorption energy. Each panel uses one functional: PBE, PBE+D3, r²SCAN, or BEEF-vdW. These are the values in the published comparison table, rather than DFT-relaxed energies.</p>
<p><a href="spe/mlip_vs_dft_spe_interactive.html">Explore all molecule and surface labels interactively</a> · <a href="spe/mlip_vs_dft_spe_zoom.pdf">PDF</a> · <a href="spe/plotted_pairs.csv">All plotted values (CSV)</a></p>
<figure><img src="spe/mlip_vs_dft_spe_zoom.png" alt="Four panels comparing MLIP adsorption energies against DFT single-point adsorption energies"><figcaption>370 PBE, 351 PBE+D3, 363 r²SCAN, and 366 BEEF-vdW pairs. Open circles retain the source table’s SPE review flags. Filled circles mean no such flag was present, not independent verification.</figcaption></figure>
<p>The zoomed view covers −5 to +1 eV; 108 points lie outside it (1 PBE+D3 and 107 BEEF-vdW). All finite pairs remain in the full-range view and downloads. Another 74 entries lack MLIP energies and are listed separately. Large differences and flagged reference or potential issues require checking before interpreting model accuracy.</p>
<details><summary>Show full range, including outliers</summary><figure><img loading="lazy" src="spe/mlip_vs_dft_spe_full_range.png" alt="Full-range four-functional comparison including all finite outliers"></figure><a href="spe/mlip_vs_dft_spe_full_range.pdf">Full-range PDF</a></details>
<p class="muted"><a href="spe/missing_ml_pairs.csv">Missing MLIP entries</a> · <a href="spe/validation.json">Numerical checks</a> · <a href="spe/plot_spe_comparison.py">Plotting script</a> · <a href="spe/sources/download_manifest.json">Source retrieval record</a></p></section>
<section id="pymatgen"><h2>Lowest converged structures from independent pymatgen starts</h2>
<p>These five selections are all <strong>PBE+D3</strong>. Each is the lowest-energy converged candidate within its own system, using matched slab and molecular references. They were not selected by comparing total energies across functionals. The audit found 54 converged candidates among 60 starting structures.</p>
<div class="scroll"><table><thead><tr><th>System</th><th>DFT E<sub>ads</sub> (eV)</th><th>Nearest O–metal (Å)</th><th>Converged / starts</th><th>Final structure</th></tr></thead><tbody>__PYM_ROWS__</tbody></table></div>
<p>DME denotes dimethyl ether, CH₃OCH₃. Structure filenames retain the initial site, anchor, and tilt for provenance; these labels do not necessarily describe the final relaxed site. Lowest among tested starts does not establish a global minimum.</p>
<p><a href="pymatgen/best_structures.csv">Selected structures and references</a> · <a href="pymatgen/all_candidates.csv">Every candidate</a> · <a href="pymatgen/selection.json">Selection audit and source paths</a></p></section>
<section id="cpu"><h2>Does an initial O–Ag distance of 2.5 Å change GOAD results?</h2>
<p class="badge">__CPU_LABEL__ · __DONE__/20 completed · __VALID__ valid endpoints · job 59224069</p>
<p>Ten paired seeds compare the original GOAD initialization against a rigid vertical shift of each initial candidate to a nearest O–Ag distance of 2.5 Å. The original rule sets a 2.3 Å vertical height, which differs from an interatomic distance. Both groups use GOAD without charge guidance and the same SevenNet-MF model (7net-mf-ompa, omat24), CPU references, and evaluation budget.</p>
<p>The distance is applied only to the initial population. Subsequent search and final relaxation have no O–Ag distance constraint. Orientations, internal conformations, and lateral positions are unchanged by the initial shift. This separates a placement effect from a change in search algorithm.</p>
<figure><img loading="lazy" src="ag_glycerol_cpu/comparison.png" alt="Paired GOAD seed comparisons of final adsorption energy and nearest O–Ag distance"><figcaption>__PAIR_COUNT__ valid matched pairs are available; the 2.5 Å initialization gives a lower final MLIP energy in __IMPROVED__ of these pairs. Unmatched completed runs are shown without connecting lines. These endpoints have not received new DFT validation.</figcaption></figure>
<p>Lower MLIP energy alone does not demonstrate improved DFT accuracy. Promising distinct endpoints need comparison under the same DFT functional and reference convention. This CPU model is separate from the SevenNet-Omni model in the SPE panels above.</p>
<p>Lowest MLIP-energy endpoint in each group (download final extended XYZ): __ENDPOINT_LINKS__. Compare these best-of-group results alongside the paired seeds, rather than treating success in some seeds as a consistent improvement.</p>
<p><a href="ag_glycerol_cpu/all_results.csv">All completed searches</a> · <a href="ag_glycerol_cpu/paired_comparison.csv">Paired differences</a> · <a href="ag_glycerol_cpu/comparison.pdf">PDF</a> · <a href="ag_glycerol_cpu/PROTOCOL.txt">Protocol</a> · <a href="ag_glycerol_cpu/status.json">Snapshot time</a></p></section>
<section id="beef"><h2>New BEEF-vdW validation on Perlmutter</h2>
<p class="badge">Scheduler snapshot: __QUEUE__ · __QUEUE_TIME__</p>
<p>71 individual regular CPU jobs were submitted: 60 independent pymatgen starts, five GOAD endpoints, and six slab or gas references. They cover Ag(111)–DME, acetaldehyde, and glycerol, plus Cu(111)–DME and ethanol.</p>
<p>This campaign uses fresh component references and BEEF-vdW with ZAB_VDW = −1.8867. Its energies must be kept separate from archived calculations with different settings. A scheduler completion by itself will not qualify a structure as converged: electronic and ionic convergence and the required references must also pass.</p>
<p><a href="beef/job_index.csv">Submitted job IDs and paths</a> · <a href="beef/queue_snapshot.json">Queue snapshot</a> · <a href="beef/PROTOCOL.md">Calculation and restart protocol</a></p></section>
<section id="kestrel"><h2>Kestrel expansion: five metals × ten molecules</h2>
<p class="badge">Inputs prepared and audited; submission not confirmed from this session.</p>
<p>Ag(111), Pt(111), Pd(111), Ir(111), and Cu(111) are crossed with glycerol, formate, formic acid, methane, O₂, NH₃, acetaldehyde, methanol, 1-propanol, and ethanol. Independent pymatgen placements sample top, bridge, and hollow sites with different anchors and tilts.</p>
<p>The prepared bundle contains <strong>660 complex relaxations, 15 component references, and four scaling tests</strong>. It uses PBE+D3, ENCUT = 450 eV, a 2 × 2 × 1 slab mesh, and a force threshold of 0.05 eV/Å. Slab cells, fixed layers, and species potentials match within each adsorption-energy calculation. New gas references include triplet O₂ and a neutral CHO₂ radical for the formate convention.</p>
<p>Ag/Pt/Pd cells contain 36 metal atoms; Cu/Ir cells contain 64. Coverage therefore differs across metals. This initial campaign uses one isolated conformer and does not exhaust molecular conformations.</p>
<p>The proposed shared-node setup uses 26 MPI ranks per complex (13 for gas references), with four tests comparing 26 and 52 ranks before choosing production resources. The default submission limit is 20 concurrent dependency lanes, without arrays. These are planned settings, not measured optimal performance.</p>
<p>Target layout: <code>/projects/ccpc/jcho/pymatgen/&lt;surface&gt;_&lt;molecule&gt;/PBE_D3/&lt;start&gt;/</code></p>
<p><a href="kestrel/case_index.csv">All 50 systems and start counts</a> · <a href="kestrel/PROTOCOL.txt">Full protocol and submission procedure</a> · <a href="kestrel/bundle_audit.json">Input audit</a></p></section>
<section><h2>Data preservation and next checks</h2><p>This dated snapshot preserves numerical results, source paths, hashes, final structures, and calculation status separately. Existing website energy tables remain unchanged. Future cluster collections should append new records and retain their functional, model, geometry, cell, potential, and component-reference provenance.</p>
<p>Next checks are to validate distinct initialization-test endpoints with DFT, collect converged BEEF-vdW outputs and references, and submit the Kestrel campaign from an authenticated session after its scaling tests. DFT remains the reference for assessing prediction accuracy.</p>
<p><a href="snapshot_manifest.json">Artifact inventory and SHA-256 hashes</a> · <a href="../../scripts/publish_research_snapshot_20261002.py">Snapshot generation script</a></p></section>
</main></body></html>
'''
for key, value in dict(LOCAL=local, DONE=status['completed'], VALID=status['valid'],
                       CPU_LABEL=cpu_label, PYM_ROWS=pym_rows, PAIR_COUNT=len(valid_pairs),
                       IMPROVED=improved, QUEUE=queue_text, QUEUE_TIME=queue['checked_utc'],
                       ENDPOINT_LINKS=endpoint_links).items():
    page = page.replace('__' + key + '__', str(value))
(OUT / 'index.html').write_text(page)
(OUT / 'snapshot_manifest.json').write_text(json.dumps(dict(published_utc=now,
    cpu_status=status, beef_queue_checked_utc=queue['checked_utc'],
    kestrel_status='prepared; no submission confirmed by this session', artifacts=inventory), indent=2))
print('Published {} artifacts to {}'.format(len(inventory), OUT))
