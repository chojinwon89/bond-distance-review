#!/usr/bin/env python3
"""Add the published single-point data to existing comparison cards.

Run from any directory: python3 scripts/update_comparison_energies.py
Exact surface/molecule/functional joins avoid conflating distinct structures.
Verified Perlmutter energies refresh relaxed cells; existing geometry is preserved.
"""
import csv
import html
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = {'PBE': 'pbe', 'PBE+D3': 'pbe_d3', 'r²SCAN': 'r2scan', 'BEEF-vdW': 'beef_vdw'}


def update():
    lookup = {}
    for row in csv.DictReader((ROOT / 'dft_singlepoint_vs_sevennet.csv').open()):
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
    for row in csv.DictReader((ROOT / 'dft_perlmutter_energy_audit.csv').open()):
        if row['publishable'] != 'true':
            continue
        key = row['surface'], aliases.get(row['molecule'], row['molecule']), row['functional']
        # If multiple adsorption sites exist, use the lowest verified energy.
        if key not in relaxed or float(row['E_ads']) < float(relaxed[key]['E_ads']):
            relaxed[key] = row
    ml = {(r['surface'], r['molecule'], r['functional']): float(r['E_ads_ML'])
          for r in csv.DictReader((ROOT / 'dft_vs_mlip_pairs.csv').open())}
    previous = {(r['surface'], r['molecule'], r['functional']): float(r['E_ads_DFT'])
                for r in csv.DictReader((ROOT / 'dft_vs_mlip_pairs.csv').open())}
    ml_by_system = {}
    for (surface, molecule, functional), value in ml.items():
        ml_by_system.setdefault((surface, molecule), set()).add(value)
    refreshed = []
    path = ROOT / 'dft_comparison.html' 
    page = path.read_text().replace('minmax(330px,1fr)', 'minmax(min(100%,480px),1fr)')
    # Make regeneration idempotent.
    page = re.sub(r'<!-- single-point-note -->.*?<!-- /single-point-note -->\n?', '', page, flags=re.S)
    page = re.sub(r'<td class="sp-energy"[^>]*>.*?</td>', '', page)
    page = page.replace('<div class="energy-scroll">', '').replace('</table><!-- energy-scroll --></div>', '</table>')
    page = re.sub(r'<tr class="energy-groups">.*?</tr>', '', page)
    page = re.sub(r'<th class="sp-energy"[^>]*>.*?</th>', '', page)
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
                if current:
                    energy = float(current['E_ads'])
                    title = html.escape('Verified Perlmutter result: ' + current['system'], quote=True)
                    tds[1] = f'<td title="{title}">{energy:.3f}</td>'
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
                        tds[3] = f'<td title="ML minus relaxed DFT">{delta:+.3f}</td>'
                    else:
                        tds[3] = '<td>&mdash;</td>'
                    refreshed.append({'surface': surface, 'molecule': molecule, 'functional': functional,
                                      'system': current['system'], 'E_ads_DFT': energy,
                                      'E_ads_ML': ml_energy, 'delta_eV': None if ml_energy is None else ml_energy-energy})
                elif 'title=' not in tds[1]:
                    tds[1] = tds[1].replace('<td', '<td title="Previously published value; not verified by the current Perlmutter extraction"', 1)
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
    note = f'''<!-- single-point-note -->
  <div class="note info"><b>Single-point adsorption energies added:</b> {len(matched)} functional results across {len(systems)} of the 415 systems below, from the site's published single-point dataset.
  <br><b>Perlmutter energy refresh (2026-09-09):</b> {len(refreshed)} verified relaxed results are included;
  {added} fill previously unavailable functional entries and {changed} update earlier values.
  Only completed, converged calculations with compatible slab references and passing the collector's energy checks are used for this refresh.
  Other relaxed values remain from the previous publication and are identified by their tooltips. Structure images are from the earlier geometry extraction.
  <a href="dft_comparison_perlmutter.csv" download>Verified relaxed energies</a> &middot;
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
    path.write_text(page)
    with (ROOT / 'dft_comparison_singlepoint.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(matched[0]))
        writer.writeheader()
        writer.writerows(matched)
    with (ROOT / 'dft_comparison_perlmutter.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(refreshed[0]))
        writer.writeheader()
        writer.writerows(refreshed)
    print(f'Added {len(matched)} single-point results across {len(systems)} systems; refreshed {len(refreshed)} relaxed results')

if __name__ == '__main__':
    update()
