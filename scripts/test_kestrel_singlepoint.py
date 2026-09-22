from decimal import Decimal
import contextlib
import io
from pathlib import Path
import re
import shutil
import tempfile
import unittest

from extract_kestrel_singlepoint import assess, infer_functional
from update_kestrel_singlepoint import combined, relaxed_cells, render


class SPEAuditTests(unittest.TestCase):
    def setUp(self):
        self.comp = dict(nsw=0, functional="pbe", status="converged", energy=100.0,
                         composition={"Ag": 36, "C": 1, "O": 2, "H": 1})
        self.slab = dict(nsw=1000, functional="pbe", status="converged", energy=80.0,
                         composition={"Ag": 36}, settings={"IBRION": "2"})
        self.mol = dict(nsw=500, functional="pbe", status="converged", energy=21.0,
                        composition={"C": 1, "O": 2, "H": 1}, settings={"IBRION": "2"})
        for c in [self.comp,self.slab,self.mol]:c['potentials']={el:'PAW_PBE '+el+' fixture' for el in c['composition']}

    def test_positive_absolute_energies_are_valid(self):
        status, _, energy = assess(self.comp, self.slab, self.mol, "pbe")
        self.assertEqual((status, energy), ("ok", -1.0))

    def test_relaxed_complex_is_not_singlepoint(self):
        self.comp["nsw"] = 1000
        self.assertEqual(assess(self.comp, self.slab, self.mol, "pbe")[0], "not_singlepoint")

    def test_reference_requires_converged_relaxation(self):
        for key, value in (("nsw", 0), ("status", "ionic_unconverged"), ("functional", "beef_vdw")):
            with self.subTest(key=key):
                reference = dict(self.slab, **{key: value})
                self.assertEqual(assess(self.comp, reference, self.mol, "pbe")[0], "reference_unavailable")

    def test_formic_acid_cannot_replace_formate_reference(self):
        self.mol["composition"]["H"] = 2
        self.assertEqual(assess(self.comp, self.slab, self.mol, "pbe")[0], "composition_mismatch")

    def test_slab_size_must_match(self):
        self.slab["composition"]["Ag"] = 64
        self.assertEqual(assess(self.comp, self.slab, self.mol, "pbe")[0], "composition_mismatch")

    def test_existing_energy_review_screen(self):
        self.comp["energy"] = 90.0
        self.assertEqual(assess(self.comp, self.slab, self.mol, "pbe")[0], "energy_review")

    def test_d3_damping_variants_and_r2scan_labels(self):
        for ivdw in ("11", "12"):
            self.assertEqual(infer_functional({"GGA": "PE", "IVDW": ivdw}), "pbe_d3")
        self.assertEqual(infer_functional({"METAGGA": "R2SCAN"}), "r2scan")
        self.assertEqual(infer_functional({"GGA": "PE", "LUSE_VDW": ".TRUE."}), "unknown")


def fixture(spe="&mdash;", delta="&mdash;"):
    return f'''<div class="g" data-surf="Ag111" data-mol="formate">
<img src="dftcmp/png/Ag111_formate_dft.png" alt="unchanged">
<table class="mini"><tr><th class="l">func</th><th>E_ads DFT</th><th>E_ads ML</th><th>&Delta;</th><th class="sp-energy">DFT</th><th class="sp-energy">ML</th><th class="sp-energy">&Delta;</th></tr><tr><td class="l">PBE</td><td title="published">-3.21</td><td>-1.98</td><td>+1.23</td><td class="sp-energy">{spe}</td><td class="sp-energy">-1.98</td><td class="sp-energy">{delta}</td></tr></table></div><script></script>'''


class SPETableTests(unittest.TestCase):
    def setUp(self):
        self.key = ("Ag111", "formate", "pbe")
        self.screened = {self.key: {"E_ads_SPE": "-3.2039268", "complex_directory": "/verified/singlepoint/PBE"}}

    def test_exact_requested_expression(self):
        self.assertEqual(combined(Decimal("-1.98"), Decimal("-3.2039268"), Decimal("-3.21")), Decimal("4.4339268"))
        self.assertIsNone(combined(None, Decimal("-3"), Decimal("-3.2")))

    def test_fill_missing_and_preserve_relaxed_cells_and_image(self):
        page = fixture()
        result, rows, additions = render(page, self.screened, {})
        self.assertEqual(additions, [self.key])
        self.assertEqual(rows[0]["ML_minus_SPE_minus_relaxed_DFT"], "4.4339268")
        self.assertEqual(relaxed_cells(result), relaxed_cells(page))
        self.assertIn('<img src="dftcmp/png/Ag111_formate_dft.png" alt="unchanged">', result)
        self.assertNotIn('<th class="sp-energy">ML</th>', result)
        previous = {self.key: rows[0]}
        again, _, additions = render(result, self.screened, previous)
        self.assertEqual(result, again)
        self.assertEqual(additions, [])

    def test_published_spe_and_delta_are_not_replaced(self):
        page = fixture("-0.123", "+0.456")
        result, rows, additions = render(page, self.screened, {})
        self.assertEqual(additions, [])
        self.assertEqual(rows[0]["E_ads_DFT_SP"], "-0.123")
        self.assertEqual(rows[0]["delta_SP"], "0.456")
        self.assertIn('>+0.456</td>', result)

    def test_review_value_is_visible_flagged_and_excluded_until_screened_replacement(self):
        reviewed={self.key:dict(E_ads_SPE='-20.123456',complex_directory='/review/SPE',status='energy_review',provenance='Perlmutter ENERGY REVIEW')}
        page,rows,_=render(fixture(),reviewed,{})
        self.assertIn('energy-review',page)
        self.assertEqual(rows[0]['SPE_review'],'true')
        self.assertEqual(rows[0]['analysis_eligible'],'false')
        self.assertEqual(rows[0]['E_ads_DFT_SP'],'-20.123456')
        again,_,_=render(page,reviewed,{self.key:rows[0]})
        self.assertEqual(page,again)
        screened={self.key:dict(self.screened[self.key],status='ok')}
        _,updated,_=render(page,screened,{self.key:rows[0]})
        self.assertEqual(updated[0]['SPE_review'],'false')
        self.assertEqual(updated[0]['E_ads_DFT_SP'],'-3.2039268')

    def test_queue_status_does_not_become_an_energy(self):
        page,rows,_=render(fixture(),{},{},{self.key:dict(label='queued',title='SPE queued; task 123')})
        self.assertEqual(rows,[])
        self.assertIn('>queued</td>',page)
        self.assertIn('task 123',page)


class RegenerationTests(unittest.TestCase):
    def test_relaxed_updater_retains_spe_additions_and_layout(self):
        import update_comparison_energies as legacy
        source = Path(__file__).resolve().parents[1]
        files = ['dft_comparison.html', 'dft_perlmutter_energy_audit.csv',
                 'dft_singlepoint_vs_sevennet.csv', 'dft_vs_mlip_pairs.csv',
                 'dft_comparison_published_relaxed.csv', 'dft_comparison_singlepoint.csv',
                 'dft_comparison_singlepoint_published.csv', 'dft_kestrel_singlepoint.csv',
                 'dft_structure_sources.json', 'dft_perlmutter_singlepoint.csv',
                 'dft_perlmutter_singlepoint_sources.json']
        files += [name for name in ['dft_cluster_selection.csv','dft_cluster_selection_summary.json','dft_gas_reference_recovery.json','dft_potential_mismatches.csv','dft_perlmutter_singlepoint_audit.csv','dft_shared_reference_energies.csv']
                  if (source/name).exists()]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in files:
                shutil.copyfile(source / name, root / name)
            original = (root / 'dft_comparison.html').read_text()
            previous_root = legacy.ROOT
            legacy.ROOT = root
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    legacy.update()
                    first = (root / 'dft_comparison.html').read_text()
                    legacy.update()
                    second = (root / 'dft_comparison.html').read_text()
            finally:
                legacy.ROOT = previous_root
            self.assertEqual(first, second)
            self.assertEqual(relaxed_cells(original), relaxed_cells(first))
            self.assertEqual(re.findall(r'<img\b[^>]*>', original), re.findall(r'<img\b[^>]*>', first))
            import csv
            with (source / 'dft_comparison_singlepoint.csv').open() as f:
                expected_count = len(list(csv.DictReader(f)))
            self.assertIn(f'{expected_count} functional results', first)
            self.assertIn('sp-combined', first)
            self.assertNotIn('<th class="sp-energy">ML</th>', first)


if __name__ == "__main__":
    unittest.main()
