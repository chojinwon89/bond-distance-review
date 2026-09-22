"""Append-only, checksum-verified publication observations, independent of validation."""
import csv
import hashlib
import json
from pathlib import Path
from bs4 import BeautifulSoup
from molecule_names import canonical

REL = Path('data/adsorption_results/publications.jsonl')


def key(record):
    return tuple(record[k] for k in ('surface', 'molecule', 'functional', 'mode'))


def checksum(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


def read(root):
    path = root / REL
    records = []
    if path.exists():
        for line in path.read_text().splitlines():
            record = json.loads(line)
            identity = record.pop('record_id')
            if checksum(record) != identity:
                raise ValueError('Publication ledger checksum mismatch')
            records.append(dict(record, record_id=identity))
    return records


def capture(root, page, source, spe_rows=()):
    from update_kestrel_singlepoint import FUNCTIONALS, number
    records = read(root)
    seen = {r['record_id'] for r in records}
    spe = {(r['surface'], canonical(r['molecule']), r['functional']): r for r in spe_rows}
    for card in BeautifulSoup(page, 'html.parser').select('.g'):
        for row in card.select('table tr'):
            cells = row.select('td')
            if len(cells) != 7:
                continue
            f = FUNCTIONALS[cells[0].get_text()]
            base = dict(surface=card['data-surf'], molecule=canonical(card['data-mol']), functional=f)
            for mode, index in [('relaxed', 1), ('SPE', 4)]:
                cell = cells[index]
                value = number(str(cell))
                if value is None or 'historical-energy' in cell.get('class', []):
                    continue
                saved = spe.get((base['surface'], base['molecule'], f), {}) if mode == 'SPE' else {}
                record = dict(base, schema_version=1, mode=mode, energy_eV=str(value),
                              displayed_text=cell.get_text(), source=source,
                              title=cell.get('title', ''), original_cell_html=str(cell),
                              spe_metadata=saved, ml_displayed=cells[2].get_text())
                identity = checksum(record)
                if identity not in seen:
                    records.append(dict(record, record_id=identity)); seen.add(identity)
    path = root / REL
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(''.join(json.dumps(r, sort_keys=True) + '\n' for r in records))
    temporary.replace(path)
    selected = {}
    for r in records:
        if r['source'] != 'local-before-selection':
            selected[key(r)] = r
    return selected


if __name__ == '__main__':
    import argparse
    import subprocess
    import io
    parser = argparse.ArgumentParser(description='Import an immutable published Git snapshot into the local ledger')
    parser.add_argument('--commit', required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(['git', 'rev-parse', args.commit], cwd=root, text=True).strip()
    def content(path):
        return subprocess.check_output(['git', 'show', commit + ':' + path], cwd=root, text=True)
    records = capture(root, content('dft_comparison.html'), 'git:' + commit,
                      list(csv.DictReader(io.StringIO(content('dft_comparison_singlepoint.csv')))))
    print('Stored keys:', len(records), 'observations:', len(read(root)))
