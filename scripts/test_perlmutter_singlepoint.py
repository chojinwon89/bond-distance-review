import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from decimal import Decimal

from extract_perlmutter_singlepoint import original_geometry, canonical
from update_kestrel_singlepoint import render, relaxed_cells
from test_kestrel_singlepoint import fixture


class GeometryProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)/'CO_Ag100'/'singlepoint'/'PBE'
        self.directory.mkdir(parents=True)
        hashes = {}
        for name in ('POSCAR', 'INCAR', 'KPOINTS'):
            (self.directory/name).write_text(name + ' original\n')
            hashes[name] = hashlib.sha256((self.directory/name).read_bytes()).hexdigest()
        self.record = dict(source=str(self.directory.parent.parent/'PBE'), destination=str(self.directory),
                           functional='PBE', geometry_source='POSCAR', source_sha256=hashes, staged_sha256=hashes)
        self.write_record()

    def write_record(self):
        (self.directory/'spe_inputs.json').write_text(json.dumps(self.record))

    def test_original_poscar_accepted(self):
        self.assertEqual(original_geometry(self.directory)['geometry_source'], 'original ML POSCAR')

    def test_relaxed_geometry_not_accepted(self):
        self.record['geometry_source'] = 'CONTCAR'; self.write_record()
        with self.assertRaises(ValueError): original_geometry(self.directory)

    def test_changed_geometry_or_settings_rejected(self):
        for name in ('POSCAR', 'INCAR', 'KPOINTS'):
            with self.subTest(name=name):
                p = self.directory/name; data = p.read_bytes(); p.write_text('changed')
                with self.assertRaises(ValueError): original_geometry(self.directory)
                p.write_bytes(data)

    def test_wrong_system_record_rejected(self):
        self.record['source'] = str(self.directory.parent.parent/'PBE_D3'); self.write_record()
        with self.assertRaises(ValueError): original_geometry(self.directory)

    def test_formula_aliases_preserve_chemical_identity(self):
        self.assertEqual(canonical('C2H6'), 'ethane')
        self.assertEqual(canonical('HCOOH'), 'formic_acid')
        self.assertNotEqual(canonical('HCOOH'), canonical('formate'))


class PerlmutterRenderingTests(unittest.TestCase):
    def test_source_and_expression_and_preservation(self):
        key = ('Ag111', 'formate', 'pbe')
        candidate = {key: dict(E_ads_SPE='-3.2039268', complex_directory='/perlmutter/singlepoint/PBE',
                               provenance='Perlmutter original ML POSCAR')}
        before = fixture()
        after, rows, additions = render(before, candidate, {})
        self.assertEqual(additions, [key])
        self.assertEqual(rows[0]['provenance'], 'Perlmutter original ML POSCAR')
        self.assertEqual(Decimal(rows[0]['ML_minus_SPE_minus_relaxed_DFT']), Decimal('4.4339268'))
        self.assertEqual(relaxed_cells(before), relaxed_cells(after))
        existing = fixture('-0.123', '+0.456')
        after, rows, additions = render(existing, candidate, {})
        self.assertEqual(additions, [])
        self.assertEqual(rows[0]['E_ads_DFT_SP'], '-0.123')
        self.assertEqual(rows[0]['delta_SP'], '0.456')


if __name__ == '__main__':
    unittest.main()
