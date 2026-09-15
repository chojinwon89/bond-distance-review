#!/usr/bin/env python3
"""Add the published single-point data to existing comparison cards.

Run from any directory: python3 scripts/update_comparison_energies.py
Exact surface/molecule/functional joins avoid conflating distinct structures.
Verified Perlmutter energies refresh relaxed cells; existing geometry is preserved.
"""
import csv
from datetime import datetime, timezone
import html
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = {'PBE': 'pbe', 'PBE+D3': 'pbe_d3', 'r²SCAN': 'r2scan', 'BEEF-vdW': 'beef_vdw'}


def read_csv(name):
    with (ROOT / name).open(newline="") as source:
        return list(csv.DictReader(source))


def update():
    kestrel_enabled = (ROOT / 'dft_kestrel_singlepoint.csv').exists()
    retained_spe_csv = (ROOT / 'dft_comparison_singlepoint.csv').read_bytes() if kestrel_enabled else None
    lookup = {}
    for row in read_csv('dft_singlepoint_vs_sevennet.csv'):
        key = row['surface'], row['molecule'], row['functional']
        if key in lookup:
            raise ValueError(f'Duplicate single-point key: {key}')
        values = tuple(float(row[k]) for k in ('E_ads_DFT', 'E_ads_ML', 'diff_eV'))
        assert all(math.isfinite(v) for v in values)
        assert abs(values[1] - values[0] - values[2]) < 2e-6
        lookup[key] = values
    aliases = {'C2H4': 'ethene', 'C2H6': 'ethane', 'CH3CH2OH': 'ethanol',
               'C2H5OH': 'ethanol', 'CH3CHO': 'acetaldehyde', 'CH3COOH': 'acetic_acid',
               'CH3OCH3': 'DME', 'CH3OH': 'methanol', 'CH4': 'methane',
               'H2CO': 'formaldehyde', 'HCOOH': 'formic_acid'}
    relaxed = {}
    review = {}
    audited = {}
    for row in read_csv('dft_perlmutter_energy_audit.csv'):
        key = row['surface'], aliases.get(row['molecule'], row['molecule']), row['functional']
        audited.setdefault(key, []).append(row)
        if row['status'] == 'energy_review' and row.get('E_ads_raw'):
            # Keep the explicit base system before site variants; never choose a
            # review value by its unusually low energy.
            if key not in review or (len(row['system'].split('_')), row['system']) < (len(review[key]['system'].split('_')), review[key]['system']):
                review[key] = row
        if row['publishable'] != 'true':
            continue
        # If multiple adsorption sites exist, use the lowest verified energy.
        if key not in relaxed or float(row['E_ads']) < float(relaxed[key]['E_ads']):
            relaxed[key] = row
    ml = {(r['surface'], r['molecule'], r['functional']): float(r['E_ads_ML'])
          for r in read_csv('dft_vs_mlip_pairs.csv')}
    previous = {(r['surface'], r['molecule'], r['functional']): float(r['E_ads_DFT'])
                for r in read_csv('dft_vs_mlip_pairs.csv')}
    published = {(r['surface'], r['molecule'], r['functional']): r
                 for r in read_csv('dft_comparison_published_relaxed.csv')}
    ml_by_system = {}
    for (surface, molecule, functional), value in ml.items():
        ml_by_system.setdefault((surface, molecule), set()).add(value)
    refreshed = []
    reviewed = []
    retained = []
    path = ROOT / 'dft_comparison.html' 
    page = path.read_text().replace('minmax(330px,1fr)', 'minmax(min(100%,480px),1fr)')
    # Make regeneration idempotent.
    page = re.sub(r'<!-- single-point-note -->.*?<!-- /single-point-note -->\n?', '', page, flags=re.S)
    page = re.sub(r'<td class="sp-energy(?: [^"]*)?"[^>]*>.*?</td>', '', page)
    page = page.replace('<div class="energy-scroll">', '').replace('</table><!-- energy-scroll --></div>', '</table>')
    page = re.sub(r'<tr class="energy-groups">.*?</tr>', '', page)
    page = re.sub(r'<th class="sp-energy(?: [^"]*)?"[^>]*>.*?</th>', '', page)
    matched, systems = [], set()

    def card(match):
        prefix, surface, molecule, body = match.groups()
        surface, molecule = html.unescape(surface), html.unescape(molecule)
        def table(tm):
            content = tm.group(1)
            def row(rm):
                cells = rm.group(1)
                label = re.search(r'<td class="l">([^<]+)</td>', cells)
                if label is None:
                    return rm.group(0)
                functional = FUNCTIONALS[html.unescape(label.group(1))]
                key = surface, molecule, functional
                current = relaxed.get(key)
                tds = re.findall(r'<td[^>]*>.*?</td>', cells)
                displayed = html.unescape(re.sub('<[^>]+>', '', tds[1]))
                try: has_displayed = math.isfinite(float(displayed))
                except ValueError: has_displayed = False
                is_review = False
                if not current and key in review and (not has_displayed or 'energy-review' in tds[1]) and key not in published:
                    current = dict(review[key], E_ads=review[key]['E_ads_raw'])
                    is_review = True
                if current:
                    energy = float(current['E_ads'])
                    annotation = ('Converged Perlmutter result; ENERGY REVIEW: |E_ads| > 5 eV; excluded from figure statistics. ' if is_review else 'Screened Perlmutter result: ') + current['system']
                    title = html.escape(annotation, quote=True)
                    marker = ' class="energy-review"' if is_review else ''
                    tds[1] = f'<td{marker} title="{title}">{energy:.3f}</td>'
                    ml_energy = ml.get(key)
                    candidates = ml_by_system.get((surface, molecule), set())
                    if ml_energy is None and len(candidates) == 1:
                        ml_energy = next(iter(candidates))
                    if ml_energy is None:
                        # The legacy table may have an ML value even when its DFT reference was absent.
                        raw = re.sub('<[^>]+>', '', tds[2])
                        try:
                            ml_energy = float(raw)
                        except ValueError:
                            pass
                    if ml_energy is not None:
                        delta = ml_energy - energy
                        delta_title = 'ML minus relaxed DFT' + ('; includes an energy-review value' if is_review else '')
                        tds[3] = f'<td{marker} title="{delta_title}">{delta:+.3f}</td>'
                    else:
                        tds[3] = '<td>&mdash;</td>'
                    result = {'surface': surface, 'molecule': molecule, 'functional': functional,
                                      'system': current['system'], 'E_ads_DFT': energy,
                                      'E_ads_ML': ml_energy, 'delta_eV': None if ml_energy is None else ml_energy-energy}
                    (reviewed if is_review else refreshed).append(result)
                else:
                    candidates = audited.get(key, [])
                    labels = {'pseudopotential_mismatch': 'POTCAR mismatch',
                              'slab_mismatch': 'Slab size mismatch',
                              'composition_mismatch': 'Composition mismatch',
                              'metadata_missing': 'Metadata missing',
                              'energy_review': 'Energy needs review',
                              'unconverged': 'Not converged', 'error': 'Output missing/incomplete'}
                    if candidates:
                        candidate = next((r for r in candidates if r['status'] == 'pseudopotential_mismatch'), candidates[0])
                        status = labels.get(candidate['status'], 'Needs review')
                        detail = ' | '.join(r['system'] + ': ' + r['note'] + '; ' + r['convergence'] for r in candidates)
                    else:
                        status, detail = 'Not audited', 'No matching result in the current Perlmutter audit'
                    # Retain a previously displayed number when the new audit cannot
                    # replace it. The committed snapshot also recovers values hidden
                    # by the earlier audit-only rendering change.
                    displayed = html.unescape(re.sub('<[^>]+>', '', tds[1]))
                    try:
                        has_displayed_energy = math.isfinite(float(displayed))
                    except ValueError:
                        has_displayed_energy = False
                    saved = published.get(key)
                    if not has_displayed_energy and saved:
                        displayed = saved['E_ads_DFT']
                        tds[2] = '<td>{}</td>'.format(html.escape(saved['E_ads_ML']))
                        tds[3] = '<td>{}</td>'.format(html.escape(saved['delta_eV']))
                    if has_displayed_energy or saved:
                        annotation = 'Previously published energy retained; current audit: ' + status + '. ' + detail
                        marker = ' class="energy-review"' if 'energy-review' in tds[1] else ''
                        tds[1] = '<td{} title="{}">{}</td>'.format(marker, html.escape(annotation, quote=True), html.escape(displayed))
                        retained.append(key)
                    else:
                        annotation = 'No previously published energy; current audit: ' + status + '. ' + detail
                        tds[1] = '<td title="{}">&mdash;</td>'.format(html.escape(annotation, quote=True))
                        tds[3] = '<td>&mdash;</td>'
                cells = ''.join(tds)
                values = lookup.get(key)
                if values is None:
                    extra = '<td class="sp-energy" title="No published single-point result for this exact system and functional">&mdash;</td>' * 3
                else:
                    matched.append({'surface':surface, 'molecule':molecule, 'functional':functional,
                                    'E_ads_DFT_SP':values[0], 'E_ads_ML_SP':values[1], 'delta_SP':values[2]})
                    systems.add((surface, molecule))
                    extra = ''.join(f'<td class="sp-energy">{v:+.3f}</td>' for v in values)
                return '<tr>' + cells + extra + '</tr>'
            content = re.sub(r'<tr>(.*?)</tr>', row, content, flags=re.S)
            content = content.replace('<th>&Delta;</th>', '<th>&Delta;</th><th class="sp-energy">DFT</th><th class="sp-energy">ML</th><th class="sp-energy">&Delta;</th>')
            groups = '<tr class="energy-groups"><th></th><th colspan="3">Relaxed comparison</th><th colspan="3">Single-point comparison</th></tr>'
            return '<div class="energy-scroll"><table class="mini">' + groups + content + '</table><!-- energy-scroll --></div>'
        return prefix + re.sub(r'<table class="mini">(.*?)</table>', table, body, flags=re.S)
    page = re.sub(r'(<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">)(.*?)(?=<div class="g" data-surf=|<script>)', card, page, flags=re.S)
    added = sum((r['surface'], r['molecule'], r['functional']) not in previous for r in refreshed)
    changed = sum(abs(r['E_ads_DFT'] - previous.get((r['surface'], r['molecule'], r['functional']), r['E_ads_DFT'])) > 2e-6 for r in refreshed)
    image_note = ('Existing structure images are from the earlier geometry extraction; added Kestrel images carry their own functional and source audit. '
                  'Geometry metrics retain the earlier extraction, and an image does not imply a new adsorption-energy result;'
                  if (ROOT / 'dft_structure_sources.json').exists() else 'Structure images are from the earlier geometry extraction;')
    note = f'''<!-- single-point-note -->
  <div class="note info"><b>Single-point adsorption energies added:</b> {len(matched)} functional results across {len(systems)} of the 415 systems below, from the site's published single-point dataset.
  <br><b>Perlmutter energy refresh ({datetime.now(timezone.utc).date().isoformat()}):</b> {len(refreshed)} screened relaxed results are included;
  {added} fill previously unavailable functional entries and {changed} update earlier values.
  These results pass completion, electronic/ionic convergence, atom-count and composition checks,
  and the |E<sub>ads</sub>| &le; 5 eV review screen.
  <br><b>Completed results outside the review range:</b> {len(reviewed)} additional relaxed values are shown in amber with a dagger (†).
  They passed convergence and composition checks but exceed |E<sub>ads</sub>| = 5 eV; they are not screened benchmark values.
  They and any derived differences are excluded from the energy-figure statistics.
  <a href="dft_comparison_review.csv" download>Displayed relaxed results requiring energy review</a>.
  Positive absolute total energies are allowed; their sign does not determine convergence or reference compatibility.
  <br><b>Project reference policy:</b> POTCAR variants such as Ag and Ag_pv are treated as equivalent for filtering.
  Potential identity differences do not exclude a result; the original potential names remain in the audit.
  No energy offset or correction is applied.
  <br><b>Previously published values preserved:</b> {len(retained)} additional relaxed energies and their published ML differences are retained
  when the current audit cannot supply a replacement. These are historical results, not newly screened results.
  Hover over an energy or dash for its provenance and current audit findings; audit findings do not replace published energies.
  A dash in a relaxed DFT cell means no retained published value or completed, composition-matched result is available.
  <a href="dft_comparison_published_relaxed.csv" download>Published relaxed-energy snapshot</a>.
  {image_note} the single-point dataset has not been re-audited against these OUTCARs.
  <a href="dft_comparison_perlmutter.csv" download>Screened relaxed energies</a> &middot;
  <a href="dft_perlmutter_energy_audit.csv" download>Full extraction and status report</a>.
  <br>All table energies are in <b>eV</b>; &Delta; = E<sub>ads</sub>(ML) &minus; E<sub>ads</sub>(DFT).
  E<sub>ads</sub> = E(slab + molecule) &minus; E(slab) &minus; E(molecule).
  The single-point columns compare energies evaluated without further geometry relaxation; the images show relaxed structures.
  Each comparison retains its own ML value. A dash means no single-point record for that exact system and functional.
  <a href="dft_comparison_singlepoint.csv" download>Download the displayed single-point energies</a> &middot;
  <a href="dft_singlepoint_vs_sevennet.csv" download>Full single-point dataset</a>.</div>
<!-- /single-point-note -->
'''
    page = page.replace('<h2>Per-system structure', note + '<h2>Per-system structure')
    if '.energy-scroll{' not in page:
        page = page.replace('</style>', '.energy-scroll{overflow-x:auto;}\n table.mini th,table.mini td{white-space:nowrap;}\n table.mini .sp-energy{background:rgba(76,120,168,.10);}\n</style>')
    if '.energy-review::after' not in page:
        page = page.replace('</style>', 'table.mini td.energy-review{color:#ffd479;background:rgba(224,168,0,.12);}\n.energy-review::after{content:" †";}\n</style>')
    path.write_text(page)
    with (ROOT / 'dft_comparison_singlepoint.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(matched[0]))
        writer.writeheader()
        writer.writerows(matched)
    with (ROOT / 'dft_comparison_perlmutter.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['surface', 'molecule', 'functional', 'system', 'E_ads_DFT', 'E_ads_ML', 'delta_eV'])
        writer.writeheader()
        writer.writerows(refreshed)
    with (ROOT / 'dft_comparison_review.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['surface', 'molecule', 'functional', 'system', 'E_ads_DFT', 'E_ads_ML', 'delta_eV'])
        writer.writeheader(); writer.writerows(reviewed)
    print(f'Added {len(matched)} single-point results across {len(systems)} systems; refreshed {len(refreshed)} relaxed results')
    if kestrel_enabled:
        from update_kestrel_singlepoint import update as update_singlepoint
        (ROOT / 'dft_comparison_singlepoint.csv').write_bytes(retained_spe_csv)
        update_singlepoint(ROOT)

if __name__ == '__main__':
    update()
