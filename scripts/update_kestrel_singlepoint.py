#!/usr/bin/env python3
"""Fill missing SPE cells and replace their ML column with ML - SPE - relaxed DFT."""
import csv
from decimal import Decimal, InvalidOperation
import html
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = {"PBE": "pbe", "PBE+D3": "pbe_d3", "r\u00b2SCAN": "r2scan", "BEEF-vdW": "beef_vdw"}
CARD = re.compile(r'(<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">)(.*?)(?=<div class="g" data-surf=|<script>)', re.S)


def number(value):
    try:
        parsed = Decimal(html.unescape(re.sub(r"<[^>]+>", "", value)).strip())
        return parsed if parsed.is_finite() else None
    except InvalidOperation:
        return None


def combined(ml, spe, relaxed):
    return ml - spe - relaxed if all(value is not None for value in (ml, spe, relaxed)) else None


def read_rows(path):
    if not path.exists():
        return {}
    with path.open(newline="") as source:
        return {(row["surface"], row["molecule"], row["functional"]): row for row in csv.DictReader(source)}


def relaxed_cells(page):
    return [re.findall(r'<td\b[^>]*>.*?</td>', row, re.S)[:4]
            for row in re.findall(r'<tr>(.*?)</tr>', page, re.S) if '<td class="l">' in row]


def render(page, screened, previous):
    displayed = []
    additions = []
    old_relaxed = relaxed_cells(page)

    def card(match):
        prefix, surface, molecule, body = match.groups()
        surface, molecule = html.unescape(surface), html.unescape(molecule)

        def table(match):
            table_body = match.group(1)

            def row(match):
                cells = re.findall(r'<td\b[^>]*>.*?</td>', match.group(1), re.S)
                if len(cells) != 7:
                    return match.group(0)
                functional = FUNCTIONALS[html.unescape(re.sub(r'<[^>]+>', '', cells[0]))]
                key = (surface, molecule, functional)
                saved = previous.get(key, {})
                candidate = screened.get(key)
                ml, relaxed = number(cells[2]), number(cells[1])
                spe = number(cells[4])
                existing_spe = spe is not None
                delta = number(cells[6])
                display_spe = html.unescape(re.sub(r'<[^>]+>', '', cells[4]))
                display_delta = html.unescape(re.sub(r'<[^>]+>', '', cells[6]))
                provenance = saved.get("provenance", "Previously published single-point result")
                source = saved.get("complex_directory", "")
                legacy_ml = saved.get("E_ads_ML_SP", "")
                if existing_spe:
                    precise = number(saved.get("E_ads_DFT_SP", ""))
                    if precise is not None and abs(precise - spe) <= Decimal("0.0005"):
                        spe = precise
                    precise_delta = number(saved.get("delta_SP", ""))
                    if delta is not None and precise_delta is not None and abs(precise_delta - delta) <= Decimal("0.0005"):
                        delta = precise_delta
                elif number(saved.get("E_ads_DFT_SP", "")) is not None:
                    spe = number(saved["E_ads_DFT_SP"])
                    display_spe = f"{spe:+.3f}"
                    delta = number(saved.get("delta_SP", ""))
                    display_delta = f"{delta:+.3f}" if delta is not None else ""
                elif candidate:
                    spe = number(candidate["E_ads_SPE"])
                    display_spe = f"{spe:+.3f}"
                    delta = ml - spe if ml is not None else None
                    display_delta = f"{delta:+.3f}" if delta is not None else ""
                    provenance = candidate.get("provenance", "Kestrel NSW=0 complex with converged relaxed slab and molecule references")
                    source = candidate["complex_directory"]
                    legacy_ml = str(ml) if ml is not None else ""
                    additions.append(key)
                value = combined(ml, spe, relaxed)
                if spe is None:
                    title = "No published SPE value or screened cluster result; see the source audits"
                    spe_cell = f'<td class="sp-energy" title="{title}">&mdash;</td>'
                else:
                    title = html.escape(provenance + ("; source: " + source if source else ""), quote=True)
                    spe_cell = f'<td class="sp-energy" title="{title}">{html.escape(display_spe)}</td>'
                    displayed.append({"surface": surface, "molecule": molecule, "functional": functional,
                        "E_ads_DFT_SP": str(spe), "E_ads_ML_SP": legacy_ml,
                        "delta_SP": str(delta) if delta is not None else "",
                        "ML_minus_SPE_minus_relaxed_DFT": str(value) if value is not None else "",
                        "E_ads_ML_relaxed_displayed": str(ml) if ml is not None else "",
                        "E_ads_DFT_relaxed_displayed": str(relaxed) if relaxed is not None else "",
                        "provenance": provenance, "complex_directory": source})
                if value is None:
                    combined_cell = '<td class="sp-energy sp-combined" title="Requires ML and DFT in the relaxed comparison and an SPE value">&mdash;</td>'
                else:
                    title = html.escape(f"ML - SPE - relaxed DFT = {ml} - ({spe}) - ({relaxed}) eV", quote=True)
                    combined_cell = f'<td class="sp-energy sp-combined" title="{title}">{value:+.3f}</td>'
                delta_title = "ML minus SPE; published differences retain their original ML reference in the download"
                delta_cell = f'<td class="sp-energy" title="{delta_title}">{html.escape(display_delta) if delta is not None else "&mdash;"}</td>'
                return '<tr>' + ''.join(cells[:4]) + spe_cell + combined_cell + delta_cell + '</tr>'

            table_body = re.sub(r'<tr>(.*?)</tr>', row, table_body, flags=re.S)
            table_body = re.sub(r'<th class="sp-energy[^\"]*"[^>]*>.*?</th>', '', table_body, flags=re.S)
            header = ('<th class="sp-energy" title="Single-point adsorption energy with relaxed references">SPE</th>'
                      '<th class="sp-energy sp-combined" title="ML minus SPE minus relaxed DFT adsorption energy">ML &minus; SPE<br>&minus; DFT (relaxed)</th>'
                      '<th class="sp-energy" title="ML minus single-point adsorption energy">&Delta;</th>')
            table_body = table_body.replace('<th>&Delta;</th>', '<th>&Delta;</th>' + header)
            return '<table class="mini">' + table_body + '</table>'

        return prefix + re.sub(r'<table class="mini">(.*?)</table>', table, body, flags=re.S)

    page = CARD.sub(card, page)
    assert relaxed_cells(page) == old_relaxed, "Relaxed values, differences or provenance changed"
    return page, displayed, additions


def update(root=ROOT):
    root = Path(root)
    page_path = root / "dft_comparison.html"
    page = page_path.read_text()
    snapshot = root / "dft_comparison_singlepoint_published.csv"
    if not snapshot.exists():
        snapshot.write_bytes((root / "dft_comparison_singlepoint.csv").read_bytes())
    previous = read_rows(snapshot)
    previous.update(read_rows(root / "dft_comparison_singlepoint.csv"))
    screened = read_rows(root / "dft_kestrel_singlepoint.csv")
    perlmutter = read_rows(root / "dft_perlmutter_singlepoint.csv")
    for key, row in perlmutter.items():
        row = dict(row, provenance="Perlmutter NSW=0 original ML POSCAR with converged relaxed slab and molecule references")
        screened.setdefault(key, row)
    page, displayed, additions = render(page, screened, previous)
    legacy_count = sum(not row["complex_directory"] for row in displayed)
    perlmutter_count = sum(row['provenance'].startswith('Perlmutter') for row in displayed)
    kestrel_count = len(displayed) - legacy_count - perlmutter_count
    systems = len({(row["surface"], row["molecule"]) for row in displayed})
    summary = (f'<b>Single-point adsorption energies:</b> {len(displayed)} functional results across {systems} of the 415 systems below; '
               f'{legacy_count} original published results, {kestrel_count} Kestrel results and {perlmutter_count} Perlmutter results; cluster additions use relaxed references.')
    page = re.sub(r'<b>Single-point adsorption energies(?: added)?:</b>.*?(?=\n  <br><b>Perlmutter)', summary, page, flags=re.S)
    page = page.replace('the single-point dataset has not been re-audited against these OUTCARs.',
                        'historical single-point values are retained; new Kestrel results have a separate source audit.')
    page = page.replace('Each comparison retains its own ML value. A dash means no single-point record for that exact system and functional.',
                        'The combined column uses the ML and DFT values displayed in the relaxed comparison. A dash means a required input is unavailable.')
    page = re.sub(r'<!-- kestrel-singlepoint-note -->.*?<!-- /kestrel-singlepoint-note -->\n?', '', page, flags=re.S)
    note = '''<!-- kestrel-singlepoint-note -->
  <div class="note info"><b>Kestrel single-point convention:</b>
  E<sub>ads</sub>(SPE) = E(complex, NSW=0) &minus; E(relaxed slab) &minus; E(relaxed molecule).
  New entries require completed, electronically converged single-point complexes and converged relaxed references
  with matching functional and composition, using the existing |E<sub>ads</sub>| &le; 5 eV review screen.
  <br><b>Combined column:</b> E<sub>ads</sub>(ML) &minus; E<sub>ads</sub>(SPE) &minus; E<sub>ads</sub>(DFT, relaxed).
  It uses the displayed relaxed-comparison inputs and the full-precision SPE value; all three must be available.
  Published SPE energies and single-point &Delta; values are preserved; their original ML references remain in the download.
  New single-point &Delta; values use the displayed relaxed-comparison ML value.
  <br><a href="dft_kestrel_singlepoint.csv" download>Screened Kestrel SPE results and relaxed reference paths</a> &middot;
  <a href="dft_kestrel_singlepoint_audit.csv" download>Candidate audit</a> &middot;
  <a href="dft_kestrel_singlepoint_sources.json">Calculation settings and source records</a> &middot;
  <a href="dft_comparison_singlepoint_published.csv" download>Original published single-point snapshot</a>.</div>
<!-- /kestrel-singlepoint-note -->
'''
    page = page.replace('<h2>Per-system structure', note + '<h2>Per-system structure')
    page = re.sub(r'<!-- perlmutter-singlepoint-note -->.*?<!-- /perlmutter-singlepoint-note -->\n?', '', page, flags=re.S)
    source_path = root / 'dft_perlmutter_singlepoint_sources.json'
    if source_path.exists():
        metadata = json.loads(source_path.read_text())
        counts = metadata['counts']
        completed = sum(c.get('status') == 'converged' for path, c in metadata['components'].items()
                        if c.get('role') == 'complex' or ('/dft_jobs/' in path and '/singlepoint/' in path))
        stamp = metadata['extracted_at'][:10]
        note = f'''<!-- perlmutter-singlepoint-note -->
  <div class="note info"><b>Perlmutter SPE refresh ({stamp}):</b> {sum(counts.values())} prepared calculations audited;
  {completed} completed with electronic convergence. {len(perlmutter)} adsorption energies pass the existing reference and energy checks;
  {perlmutter_count} fill previously missing table entries. All earlier published SPE values and differences are retained.
  <br>These calculations use the <b>original ML POSCAR</b>, verified against the recorded staging hashes, with NSW=0.
  E<sub>ads</sub>(SPE) = E(complex) &minus; E(relaxed slab) &minus; E(relaxed molecule), using the same functional and matching composition.
  POTCAR variants follow the existing project policy; no scaling, offsets or sign-based rejection is applied.
  Missing, incomplete and reference-limited results remain in the audit; the SPE batch is still in progress at this snapshot.
  <br><a href="dft_perlmutter_singlepoint.csv" download>Screened Perlmutter SPE energies and reference paths</a> &middot;
  <a href="dft_perlmutter_singlepoint_audit.csv" download>All SPE candidates and component total energies</a> &middot;
  <a href="dft_perlmutter_singlepoint_sources.json">Settings, convergence and geometry provenance</a>.</div>
<!-- /perlmutter-singlepoint-note -->
'''
        page = page.replace('<h2>Per-system structure', note + '<h2>Per-system structure')
    page_path.write_text(page)
    with (root / "dft_comparison_singlepoint.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(displayed[0]))
        writer.writeheader()
        writer.writerows(displayed)
    print(f"Displayed {len(displayed)} SPE results ({kestrel_count} Kestrel, {perlmutter_count} Perlmutter, {legacy_count} original published); filled {len(additions)} missing cells this run.")


if __name__ == "__main__":
    update()
