import unittest
import contextlib
import io
from pathlib import Path
import re
import shutil
import tempfile
from molecule_names import canonical,equivalent_names


class MoleculeNameTests(unittest.TestCase):
    def test_formula_and_name_equivalences(self):
        for name,expected in [('CH3OCH3','DME'),('CH₃OCH₃','DME'),('dimethyl ether','DME'),
                              ('C2H5OH','ethanol'),('CH3CH2OH','ethanol'),('CH3OH','methanol'),
                              ('CH3CHO','acetaldehyde'),('CH3COOH','acetic_acid'),('H2CO','formaldehyde'),
                              ('HCOOH','formic_acid'),('ethylene','ethene'),('CH3O','methoxy')]:
            with self.subTest(name=name):self.assertEqual(canonical(name),expected)

    def test_isomers_radicals_and_distinct_stoichiometries_stay_separate(self):
        for a,b in [('CH3OCH3','CH3CH2OH'),('formate','HCOOH'),('CH3O','CH2OH'),
                    ('propanol','isopropanol'),('CO','Co')]:
            self.assertNotEqual(canonical(a),canonical(b))
        self.assertEqual(canonical('C2H6O'),'C2H6O')
        self.assertEqual(canonical('C3H8O'),'C3H8O')

    def test_reference_alias_order_preserves_spelling_and_excludes_isomer(self):
        names=['ethanol','DME','CH3OCH3','CH3CH2OH']
        self.assertEqual(equivalent_names('CH3OCH3',names),['CH3OCH3','DME'])
        self.assertEqual(equivalent_names('C2H5OH',names,preferred='CH3CH2OH'),['CH3CH2OH','ethanol'])

    def test_notation_regeneration_preserves_all_energy_tables(self):
        from update_molecule_notation import update
        source=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in ['dft_comparison.html','dft_perlmutter_energy_audit.csv','dft_perlmutter_singlepoint_audit.csv']:
                shutil.copyfile(source/name,root/name)
            (root/'dft_jobs'/'CH3OCH3_Ag100').mkdir(parents=True)
            (root/'vasp_mol'/'DME').mkdir(parents=True)
            before=(root/'dft_comparison.html').read_text()
            with contextlib.redirect_stdout(io.StringIO()):
                update(root,root)
                first=(root/'dft_comparison.html').read_text()
                update(root,root)
            self.assertEqual(re.findall(r'<table\b.*?</table>',before,re.S),re.findall(r'<table\b.*?</table>',first,re.S))
            self.assertEqual((root/'dft_comparison.html').read_text(),first)
            self.assertIn('DME (CH3OCH3)',first)


if __name__=='__main__':unittest.main()
