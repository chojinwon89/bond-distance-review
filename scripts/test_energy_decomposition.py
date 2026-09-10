import unittest
from pathlib import Path
from plot_energy_decomposition import decompose, load_rows, common_spe, FUNCTIONALS


class DecompositionTests(unittest.TestCase):
    def test_identity(self):
        for values in [(-1, -1, -1), (-1.61, -1.8390008, -1.857), (0.2, -2, -1)]:
            a, b, total = decompose(*values)
            self.assertAlmostEqual(a + b, total)
            self.assertAlmostEqual(total - a, b)

    def test_perfect_match(self):
        self.assertEqual(decompose(-1, -1, -1), (0, 0, 0))

    def test_published_pairs(self):
        rows = load_rows(Path(__file__).resolve().parents[1] / 'dft_comparison_singlepoint.csv')
        keys = [(r['surface'], r['molecule'], r['functional']) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertGreater(len(rows), 0)
        for r in rows:
            self.assertAlmostEqual(r['residual'] + r['relaxation_gap'], r['total'])

    def test_functional_cohort(self):
        cohort = common_spe(Path(__file__).resolve().parents[1] / 'dft_comparison_singlepoint.csv')
        self.assertGreater(len(cohort), 0)
        for row in cohort:
            self.assertTrue(set(FUNCTIONALS).issubset(row))


if __name__ == '__main__':
    unittest.main()
