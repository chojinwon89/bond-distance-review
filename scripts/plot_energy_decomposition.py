#!/usr/bin/env python3
"""Descriptive paired analysis; explicitly marked energy-review values stay in tables only."""
import csv
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = {'pbe': 'PBE', 'pbe_d3': 'PBE+D3', 'r2scan': 'r2SCAN', 'beef_vdw': 'BEEF-vdW'}
COLORS = ['#167b9a', '#b65172', '#398354', '#9c7023']


def decompose(ml, spe, relaxed):
    return ml - spe, spe - relaxed, ml - relaxed


def load_rows(path):
    rows = []
    with path.open(newline='') as source:
        for row in csv.DictReader(source):
            if row.get('relaxed_review') == 'true' or row.get('SPE_review') == 'true':
                continue
            try:
                ml, spe, relaxed = [float(row[k]) for k in (
                    'E_ads_ML_relaxed_displayed', 'E_ads_DFT_SP', 'E_ads_DFT_relaxed_displayed')]
            except ValueError:
                continue
            if not np.isfinite([ml, spe, relaxed]).all():
                continue
            a, b, c = decompose(ml, spe, relaxed)
            rows.append(dict(surface=row['surface'], molecule=row['molecule'], functional=row['functional'],
                             ML=ml, SPE=spe, relaxed=relaxed, residual=a, relaxation_gap=b, total=c,
                             provenance=row['provenance'], complex_directory=row['complex_directory']))
    return rows


def common_spe(path):
    systems = {}
    with path.open(newline='') as source:
        for row in csv.DictReader(source):
            if row.get('SPE_review') == 'true':
                continue
            try:
                energy = float(row['E_ads_DFT_SP'])
            except ValueError:
                continue
            if np.isfinite(energy):
                systems.setdefault((row['surface'], row['molecule']), {})[row['functional']] = energy
    return [dict(surface=s, molecule=m, **v) for (s, m), v in systems.items()
            if set(v) == set(FUNCTIONALS)]


def save(fig, name):
    fig.savefig(ROOT / name, dpi=180, facecolor='white', bbox_inches='tight')
    plt.close(fig)


def main():
    source = ROOT / 'dft_comparison_singlepoint.csv'
    with source.open(newline='') as f: source_rows=list(csv.DictReader(f))
    review_excluded=sum(r.get('relaxed_review')=='true' or r.get('SPE_review')=='true' for r in source_rows)
    rows = load_rows(source)
    grouped = {f: [r for r in rows if r['functional'] == f] for f in FUNCTIONALS}
    by_system = {}
    for r in rows:
        by_system.setdefault((r['surface'], r['molecule']), {})[r['functional']] = r
    common = common_spe(source)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(3, 1, figsize=(9, 12), layout='constrained')
    for ax, key, title in zip(axes, ['residual', 'relaxation_gap', 'total'],
                              ['ML - SPE: model / method / reference residual',
                               'SPE - relaxed DFT: observed relaxation gap',
                               'ML - relaxed DFT: total discrepancy']):
        for i, (f, label) in enumerate(FUNCTIONALS.items()):
            values = [r[key] for r in grouped[f]]
            box = ax.boxplot([values], positions=[i], widths=.48, patch_artist=True,
                             flierprops={'markersize': 3, 'alpha': .5})
            box['boxes'][0].set_facecolor(COLORS[i])
        ax.set_xticks(range(4), [f'{label}\nn={len(grouped[f])}' for f, label in FUNCTIONALS.items()])
        ax.axhline(0, color='#555555', linewidth=.8, linestyle='--')
        ax.set_ylabel('Energy difference (eV)')
        ax.set_title(title, loc='left', fontsize=13)
        ax.grid(axis='y', alpha=.2)
    save(fig, 'energy_decomposition_distributions.png')

    fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
    shifts = {}
    for i, f in enumerate(list(FUNCTIONALS)[1:]):
        values = [v[f] - v['pbe'] for v in common]
        shifts[f] = {'n': len(values), 'mean': float(np.mean(values)), 'median': float(np.median(values)),
                     'mean_absolute_shift': float(np.mean(np.abs(values)))}
        box = ax.boxplot([values], positions=[i], widths=.45, patch_artist=True,
                         flierprops={'markersize': 4, 'alpha': .6})
        box['boxes'][0].set_facecolor(COLORS[i + 1])
    ax.set_xticks(range(3), [FUNCTIONALS[f] for f in list(FUNCTIONALS)[1:]])
    ax.axhline(0, color='#555555', linestyle='--', linewidth=.8)
    ax.set_ylabel('SPE(functional) - SPE(PBE) (eV)')
    ax.set_title(f'Functional sensitivity: same {len(common)} systems for every comparison', loc='left', fontsize=13)
    ax.grid(axis='y', alpha=.2)
    save(fig, 'energy_functional_shifts.png')

    fig, axes = plt.subplots(4, 1, figsize=(8, 17), layout='constrained', sharex=True, sharey=True)
    bound = max(max(abs(r['residual']), abs(r['relaxation_gap'])) for r in rows) * 1.08
    for ax, (f, label), color in zip(axes, FUNCTIONALS.items(), COLORS):
        values = grouped[f]
        ax.scatter([r['residual'] for r in values], [r['relaxation_gap'] for r in values],
                   s=18, color=color, alpha=.65, edgecolors='none')
        ax.plot([-bound, bound], [bound, -bound], '--', color='#777777', linewidth=1,
                label='Total discrepancy = 0')
        ax.axhline(0, color='#cccccc', linewidth=.8)
        ax.axvline(0, color='#cccccc', linewidth=.8)
        ax.set(xlim=(-bound, bound), ylim=(-bound, bound), ylabel='SPE - relaxed DFT (eV)',
               title=f'{label} | n={len(values)}')
        ax.legend(loc='upper right', fontsize=9)
    axes[-1].set_xlabel('ML - SPE (eV)')
    save(fig, 'energy_decomposition_scatter.png')

    summaries = []
    for f, label in FUNCTIONALS.items():
        v = grouped[f]
        summaries.append(dict(functional=label, n=len(v),
            mean_residual=float(np.mean([r['residual'] for r in v])),
            mae_residual=float(np.mean([abs(r['residual']) for r in v])),
            mean_gap=float(np.mean([r['relaxation_gap'] for r in v])),
            mae_total=float(np.mean([abs(r['total']) for r in v])),
            negative_gap_count=sum(r['relaxation_gap'] < -.01 for r in v)))
    for name, data in [('energy_decomposition.csv', rows), ('energy_decomposition_summary.csv', summaries),
                       ('energy_functional_cohort.csv', common)]:
        with (ROOT / name).open('w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    metadata = dict(source=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    displayed_SPE_rows=len(source_rows),energy_review_rows_excluded=review_excluded,
                    paired_rows=len(rows), systems=len(by_system), common_four_functional_systems=len(common),
                    summary=summaries, functional_shifts=shifts,
                    max_identity_error=max(abs(r['residual'] + r['relaxation_gap'] - r['total']) for r in rows))
    (ROOT / 'energy_decomposition_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    table = ''.join(f'<tr><td>{s["functional"]}</td><td>{s["n"]}</td>'
                    f'<td>{s["mean_residual"]:+.3f}</td><td>{s["mae_residual"]:.3f}</td>'
                    f'<td>{s["mean_gap"]:+.3f}</td><td>{s["mae_total"]:.3f}</td></tr>' for s in summaries)
    shift_text = '; '.join(f'{FUNCTIONALS[f]}: {v["mean"]:+.3f} eV mean ({v["median"]:+.3f} median)'
                           for f, v in shifts.items())
    section = f'''<!-- energy-decomposition -->
<section id="energy-decomposition" aria-labelledby="energy-analysis-title">
<h2 id="energy-analysis-title">Energy differences: method residual and relaxation</h2>
<p><b>(ML &minus; relaxed DFT) = (ML &minus; SPE) + (SPE &minus; relaxed DFT).</b>
All terms here are adsorption energies in eV. Subtracting the same-geometry residual from the total leaves
<b>(ML &minus; relaxed DFT) &minus; (ML &minus; SPE) = SPE &minus; relaxed DFT</b>, not ML &minus; relaxed DFT &minus; SPE.</p>
<p><b>Interpretation:</b> ML &minus; SPE measures a same-geometry model discrepancy only when the ML and SPE structures
and adsorption-reference conventions match. It includes ML approximation error, functional differences and reference offsets;
it is not a pure functional correction. SPE &minus; relaxed DFT is a relaxation contribution only with compatible references,
potentials, numerical settings and a linked initial/final geometry. It is not a functional-independent "true error".</p>
<p><b>Current evidence is descriptive:</b> {len(rows)} complete ML/SPE/relaxed triples across {len(by_system)} exact
surface/molecule pairs out of {len(source_rows)} displayed SPE entries. Missing inputs are omitted, never treated as zero.
{review_excluded} rows with explicitly marked energy-review values are excluded from these paired statistics; the tables retain their values.
All panels use the same triples within each functional. Historical and audited results are mixed; geometry identity and
cross-run reference compatibility have not been established. Kestrel SPE values currently subtract relaxed molecule/slab
references; the newly submitted reference SPE jobs are not included. POTCAR variants were allowed by the existing publication
policy, which does not establish physical reference equivalence. The existing energy screen can bias these distributions.</p>
<p>One ML value is used consistently per row: the displayed relaxed-comparison ML energy. Historical SPE differences sometimes
used another ML reference; those published differences remain unchanged and are not used in these plots.
Displayed ML/relaxed precision limits the derived differences. Formate and formic acid remain separate.</p>
<h3>Paired discrepancy distributions</h3>
<img src="energy_decomposition_distributions.png" alt="Three box plots of ML minus SPE, SPE minus relaxed DFT, and total discrepancy by functional" loading="lazy" style="width:100%;height:auto;max-width:850px">
<p class="legend">Median line; box = 25th&ndash;75th percentiles; whiskers = 1.5 interquartile ranges; all outliers shown. Different functionals have different coverage.</p>
<div class="energy-scroll"><table class="big"><thead><tr><th>Functional</th><th>n</th><th>Mean ML&minus;SPE</th><th>MAE ML&minus;SPE</th><th>Mean SPE&minus;relaxed</th><th>MAE ML&minus;relaxed</th></tr></thead><tbody>{table}</tbody></table></div>
<h3>Functional sensitivity on a common cohort</h3>
<img src="energy_functional_shifts.png" alt="Paired single-point adsorption energy shifts relative to PBE on the common four-functional cohort" loading="lazy" style="width:100%;height:auto;max-width:850px">
<p>{len(common)} systems have SPE values for all four functionals; ML and relaxed energies are not required for this panel. {shift_text}.
Negative shifts mean stronger adsorption than PBE. Pairing removes system-composition imbalance, but does not prove identical
geometry or settings. These are observed method/reference shifts, not isolated exchange-correlation effects or an accuracy ranking.</p>
<h3>Cancellation versus reinforcement</h3>
<img src="energy_decomposition_scatter.png" alt="Scatter plots of residual against observed relaxation gap; the diagonal marks zero total discrepancy" loading="lazy" style="width:100%;height:auto;max-width:720px">
<p>Each point is an exact surface/molecule/functional triple. Its x + y value is the total ML&minus;relaxed discrepancy.
Opposite signs cancel; equal signs reinforce. A negative gap below &minus;0.01 eV occurs in
{sum(s['negative_gap_count'] for s in summaries)} triples: these require provenance/geometry review, not clipping.
For genuinely matched calculations relaxed to a lower energy with identical references, the gap should be nonnegative.</p>
<p><b>Legacy-column warning:</b> the retained ML &minus; SPE &minus; DFT (relaxed) column is literal arithmetic,
not an error decomposition or a corrected physical error. It is excluded from this analysis.</p>
<p><a href="energy_decomposition.csv" download>Paired data and provenance</a> &middot;
<a href="energy_decomposition_summary.csv" download>Summary statistics</a> &middot;
<a href="energy_functional_cohort.csv" download>Common functional cohort</a> &middot;
<a href="energy_decomposition_metadata.json">Cohort counts and source checksum</a> &middot;
<a href="https://vasp.at/wiki/index.php/NSW">VASP ionic-step definition</a> &middot;
<a href="scripts/ENERGY_DECOMPOSITION.md">Analysis methodology</a></p>
</section>
<!-- /energy-decomposition -->
'''
    path = ROOT / 'dft_comparison.html'
    page = path.read_text()
    page = re.sub(r'<!-- energy-decomposition -->.*?<!-- /energy-decomposition -->\n?', '', page, flags=re.S)
    page = page.replace('<h2>Per-system structure', section + '<h2>Per-system structure')
    path.write_text(page)
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
