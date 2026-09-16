"""Regression checks for the reference-mismatch audit (standard library only)."""
import copy
import gzip
import tempfile
import unittest
from pathlib import Path
from extract_perlmutter_energies import assess, parse_output, read_output


def component(energy, composition, potentials, convergence='converged'):
    return dict(energy=energy, composition=composition, potentials=potentials,
                convergence=convergence, metadata_ok=True)


class ReferenceAuditTests(unittest.TestCase):
    def setUp(self):
        self.parts = [component(7, {'Ag': 36, 'C': 1, 'H': 3}, {'Ag': 'Ag', 'C': 'C', 'H': 'H'}),
                      component(20, {'Ag': 36}, {'Ag': 'Ag'}),
                      component(-12, {'C': 1, 'H': 3}, {'C': 'C', 'H': 'H'})]

    def test_positive_absolute_energies_can_be_valid(self):
        self.assertEqual(assess(*self.parts, metal='Ag'), ('ok', '', -1))

    def test_metal_potential_variants_are_rejected(self):
        self.parts[1]['potentials']['Ag'] = 'Ag_pv'
        status, note, raw = assess(*self.parts, metal='Ag')
        self.assertEqual(status, 'potential_mismatch')
        self.assertIn('potential',note)
        self.assertEqual(raw, -1)

    def test_gas_potential_variants_are_rejected(self):
        self.parts[2]['potentials']['C'] = 'C_h'
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'potential_mismatch')

    def test_atom_counts_and_composition(self):
        parts = copy.deepcopy(self.parts)
        parts[1]['composition']['Ag'] = 64
        self.assertEqual(assess(*parts, metal='Ag')[0], 'slab_mismatch')
        self.parts[2]['composition']['H'] = 4
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'composition_mismatch')

    def test_incomplete_reference_and_missing_metadata(self):
        self.parts[2]['convergence'] = 'unfinished'
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'unconverged')
        self.parts[1]['metadata_ok'] = False
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'metadata_missing')

    def test_outlier_is_review_not_convergence_failure(self):
        self.parts[0]['energy'] = 40
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'energy_review')
        self.parts[0]['energy'] = float('nan')
        self.assertEqual(assess(*self.parts, metal='Ag')[0], 'error')

    def test_final_electronic_step_must_converge(self):
        head = 'NSW = 10\n TITEL = PAW_PBE Ag 02Apr2005\n VRHFIN =Ag: d10\n ions per type = 36\n'
        tail = ('aborting loop because EDIFF is reached\nfree  energy   TOTEN = 3 eV\n'
                'free  energy   TOTEN = 2 eV\nreached required accuracy\nGeneral timing and accounting')
        self.assertEqual(parse_output(head, tail)['convergence'], 'electronic_unconverged')
        tail = tail.replace('free  energy   TOTEN = 2', 'aborting loop because EDIFF is reached\nfree  energy   TOTEN = 2')
        parsed = parse_output(head, tail)
        self.assertEqual(parsed['convergence'], 'converged')
        self.assertEqual(parsed['potentials']['Ag'], 'PAW_PBE Ag 02Apr2005')
        self.assertEqual(parsed['composition'], {'Ag': 36})
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'OUTCAR'
            with gzip.open(str(path) + '.gz', 'wt') as f:
                f.write(head + tail)
            self.assertEqual(read_output(path)['energy'], 2)
            self.assertTrue(read_output(path)['path'].endswith('.gz'))


if __name__ == '__main__':
    unittest.main()
