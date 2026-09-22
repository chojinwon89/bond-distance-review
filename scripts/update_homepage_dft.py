"""Refresh gallery DFT references from validated selections, without changing ML energies."""
import csv
import json
import re
from pathlib import Path
from component_store import read_store
from molecule_names import canonical

ROOT = Path(__file__).resolve().parents[1]
FUNCTIONALS = ['r2scan', 'pbe_d3', 'pbe', 'beef_vdw']


def choose_reference(row, candidates, snapshots):
    compatible = []
    for candidate in candidates:
        if candidate.get('status') not in ('ok', 'energy_review') or not candidate.get('E_ads'):
            continue
        record = snapshots[candidate['complex_snapshot_id']]
        comp = record['calculation']['composition']
        metal = re.match(r'[A-Z][a-z]?', row['surface']).group()
        metal_count = comp.get(metal, 0)
        if metal_count != row['n_metal'] or sum(comp.values()) - metal_count != row['n_mol']:
            continue
        compatible.append(candidate)
    # A relaxed reference is preferred. Preserve the card's functional when
    # available in that mode; otherwise use a fixed order, never energy ranking.
    compatible.sort(key=lambda c: (c['mode'] != 'relaxed',
        c['functional'] != row.get('eads_dft_func'), FUNCTIONALS.index(c['functional'])))
    return compatible[0] if compatible else None


def main(root=ROOT):
    page = (root/'index.html').read_text()
    marker = 'const ROWS = '
    start = page.index(marker) + len(marker)
    rows, length = json.JSONDecoder().raw_decode(page[start:])
    snapshots = {r['snapshot_id']: r for r in read_store(root/'dft_component_store.jsonl.gz')}
    grouped = {}
    with (root/'dft_cluster_selection.csv').open() as f:
        for r in csv.DictReader(f):
            if r.get('application_status') == 'selected candidate applied':
                grouped.setdefault((r['surface'], canonical(r['molecule'])), []).append(r)
    audit = []
    for r in rows:
        selected = choose_reference(r, grouped.get((r['surface'], canonical(r['molecule'])), []), snapshots)
        old, old_func = r.get('eads_dft'), r.get('eads_dft_func', '')
        r['eads_dft'] = float(selected['E_ads']) if selected else None
        r['eads_dft_func'] = selected['functional'] if selected else ''
        r['eads_dft_mode'] = selected['mode'] if selected else ''
        r['eads_dft_review'] = bool(selected and selected['status'] == 'energy_review')
        r['eads_dft_source'] = selected['complex_cluster'] if selected else ''
        if selected or old is not None:
            audit.append(dict(name=r['name'],surface=r['surface'],molecule=r['molecule'],
                n_metal=r['n_metal'],n_mol=r['n_mol'],old_energy=old,old_functional=old_func,
                new_energy=r['eads_dft'],functional=r['eads_dft_func'],mode=r['eads_dft_mode'],
                status=selected['status'] if selected else 'no validated size-matched selection',
                complex_snapshot_id=selected['complex_snapshot_id'] if selected else ''))
    page = page[:start] + json.dumps(rows,ensure_ascii=False) + page[start+length:]
    page = page.replace("if(v===null) return '<div class=\"eads none\">E<sub>ads</sub> &mdash;</div>';", '')
    old = 'let s=`<div class="eads" style="color:${col}" title="${tip}">E<sub>ads</sub> ${v.toFixed(2)} <small>eV</small></div>`;'
    new = 'let s=v===null?\'<div class="eads none">E<sub>ads</sub> &mdash;</div>\':`<div class="eads" style="color:${col}" title="${tip}">E<sub>ads</sub> ${v.toFixed(2)} <small>eV</small></div>`;'
    assert old in page or new in page
    page = page.replace(old,new)
    old = 's+=`<div class="eadsdft">DFT ${r.eads_dft_func} ${Number(r.eads_dft).toFixed(2)} eV</div>`;'
    new = 's+=`<div class="eadsdft" title="DFT reference for the same system and atom counts; this is not a claim of identical gallery coordinates. Source: ${r.eads_dft_source}">DFT ref. ${r.eads_dft_func} (${r.eads_dft_mode}) ${Number(r.eads_dft).toFixed(2)} eV${r.eads_dft_review?" · review":""} <a href="dft_comparison.html">details</a></div>`;'
    assert old in page or new in page
    page = page.replace(old,new)
    page = re.sub(r'(?m)^[ \t]+$', '', page)
    (root/'index.html').write_text(page)
    with (root/'homepage_dft_refresh.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    print('Gallery DFT references:',sum(r['eads_dft'] is not None for r in rows),'of',len(rows))


if __name__ == '__main__':
    main()
