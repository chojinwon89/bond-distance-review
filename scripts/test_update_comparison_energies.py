"""Regression checks that audits annotate, rather than erase, published energies."""
import contextlib
import csv
import html
import io
import math
import re
import shutil
import tempfile
import unittest
from pathlib import Path
import update_comparison_energies as updater

ROOT = Path(__file__).resolve().parents[1]


def table_values(page):
    result = {}
    for surface, molecule, body in re.findall(r'<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">(.*?)(?=<div class="g" data-surf=|<script>)', page, re.S):
        for row in re.findall(r'<tr>(.*?)</tr>', body, re.S):
            cells = re.findall(r'<td[^>]*>(.*?)</td>', row, re.S)
            if cells:
                cells = [html.unescape(re.sub('<[^>]+>', '', c)) for c in cells]
                result[(html.unescape(surface), html.unescape(molecule), updater.FUNCTIONALS[cells[0]])] = cells[1:4]
    return result


class PublishedEnergyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        for name in ['dft_comparison.html', 'dft_perlmutter_energy_audit.csv',
                     'dft_singlepoint_vs_sevennet.csv', 'dft_vs_mlip_pairs.csv',
                     'dft_comparison_published_relaxed.csv']:
            shutil.copyfile(str(ROOT / name), str(self.root / name))
        self.original_root = updater.ROOT
        updater.ROOT = self.root

    def tearDown(self):
        updater.ROOT = self.original_root
        self.temp.cleanup()

    def update(self):
        with contextlib.redirect_stdout(io.StringIO()):
            updater.update()

    def test_restore_all_published_numbers_and_idempotence(self):
        before = table_values((self.root / 'dft_comparison.html').read_text())
        self.update()
        page = (self.root / 'dft_comparison.html').read_text()
        values = table_values(page)
        with (self.root / 'dft_comparison_published_relaxed.csv').open() as f:
            published = list(csv.DictReader(f))
        with (self.root / 'dft_comparison_perlmutter.csv').open() as f:
            refreshed = {(r['surface'], r['molecule'], r['functional']) for r in csv.DictReader(f)}
        for row in published:
            key = (row['surface'], row['molecule'], row['functional'])
            self.assertTrue(math.isfinite(float(values[key][0])))
            if key not in refreshed:
                self.assertEqual(values[key], before[key])
        self.assertEqual(len(published), 1087)
        self.assertNotIn('class="audit-status"', page)
        self.assertIn('Previously published energy retained; current audit:', page)
        self.update()
        self.assertEqual((self.root / 'dft_comparison.html').read_text(), page)

    def test_later_audit_failure_keeps_displayed_number_and_delta(self):
        self.update()
        before = table_values((self.root / 'dft_comparison.html').read_text())
        path = self.root / 'dft_perlmutter_energy_audit.csv'
        with path.open() as f:
            rows = list(csv.DictReader(f))
        for target in rows:
            target.update(publishable='false', status='unconverged', note='Reference rerun in progress')
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        self.update()
        after = table_values((self.root / 'dft_comparison.html').read_text())
        self.assertEqual(before, after)
        self.assertIn('Reference rerun in progress', (self.root / 'dft_comparison.html').read_text())

    def test_review_result_fills_blank_with_flag_but_not_screened_export(self):
        from bs4 import BeautifulSoup
        path=self.root/'dft_comparison.html'
        page=path.read_text()
        pattern=r'(<div class="g" data-surf="Rh100" data-mol="CH3">)(.*?)(?=<div class="g" data-surf=|<script>)'
        def blank(match):
            body=re.sub(r'(<tr><td class="l">BEEF-vdW</td>)<td[^>]*>.*?</td>',r'\1<td>&mdash;</td>',match.group(2),count=1)
            return match.group(1)+body
        path.write_text(re.sub(pattern,blank,page,flags=re.S))
        self.update()
        page=BeautifulSoup(path.read_text(),'html.parser')
        card=page.select_one('.g[data-surf="Rh100"][data-mol="CH3"]')
        row=next(r for r in card.select('tr') if r.select('td') and r.select('td')[0].get_text()=='BEEF-vdW')
        self.assertIn('energy-review',row.select('td')[1].get('class',[]))
        self.assertEqual(row.select('td')[1].get_text(),'-28.388')
        with (self.root/'dft_comparison_perlmutter.csv').open() as f:
            self.assertFalse(any(r['system']=='CH3_Rh100' and r['functional']=='beef_vdw' for r in csv.DictReader(f)))


if __name__ == '__main__':
    unittest.main()
