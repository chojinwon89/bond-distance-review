import unittest

import numpy as np

from add_missing_dft_images import final_positions_forces, incar_values, normalize_setting


class StructureAuditTests(unittest.TestCase):
    def test_negative_coordinates_do_not_end_force_table(self):
        text = """ POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      2.04500 1.44603 15.00000 0.000374 -0.006203 0.170165
     -0.00311 0.00184 19.37421 -0.003350 -0.007901 0.027463
 -----------------------------------------------------------------------------------
 total drift: 0 0 0
"""
        data = final_positions_forces(text + text.replace("-0.00311", "-0.00999"))
        self.assertEqual(data.shape, (2, 6))
        self.assertAlmostEqual(data[1, 0], -0.00999)
        self.assertTrue(np.isfinite(data).all())

    def test_incar_comments_do_not_disable_actual_relaxation(self):
        settings = incar_values("! single-point\nNSW = 1000; IBRION = 2 # NSW = 0\nEDIFFG = -5E-02\n")
        self.assertEqual(settings, {"NSW": "1000", "IBRION": "2", "EDIFFG": "-5E-02"})

    def test_vasp_truncated_metagga_and_boolean_settings(self):
        self.assertEqual(normalize_setting("METAGGA", "R2SCA"), "R2SCAN")
        self.assertNotEqual(normalize_setting("METAGGA", "SCAN"), "R2SCAN")
        self.assertEqual(normalize_setting("LUSE_VDW", ".TRUE."), "T")

if __name__ == "__main__":
    unittest.main()
