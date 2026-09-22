import unittest
from update_homepage_dft import choose_reference


class HomepageReferenceTests(unittest.TestCase):
    def test_size_and_validation_precede_functional_preference(self):
        row=dict(surface='Ag111',n_metal=36,n_mol=9,eads_dft_func='r2scan')
        records={'small':{'calculation':{'composition':{'Ag':36,'C':2,'H':6,'O':1}}},
                 'large':{'calculation':{'composition':{'Ag':64,'C':2,'H':6,'O':1}}}}
        candidates=[dict(status='ok',E_ads='-8',mode='relaxed',functional='r2scan',complex_snapshot_id='large'),
                    dict(status='ok',E_ads='-.3',mode='relaxed',functional='pbe',complex_snapshot_id='small')]
        self.assertEqual(choose_reference(row,candidates,records)['functional'],'pbe')
        candidates[1]['status']='potential_mismatch'
        self.assertIsNone(choose_reference(row,candidates,records))

    def test_review_value_retained_without_energy_ranking(self):
        row=dict(surface='Ag111',n_metal=36,n_mol=9,eads_dft_func='beef_vdw')
        records={'s':{'calculation':{'composition':{'Ag':36,'C':2,'H':6,'O':1}}}}
        candidates=[dict(status='energy_review',E_ads='-9',mode='relaxed',functional='beef_vdw',complex_snapshot_id='s'),
                    dict(status='ok',E_ads='-.5',mode='SPE',functional='pbe',complex_snapshot_id='s')]
        self.assertEqual(choose_reference(row,candidates,records)['E_ads'],'-9')


class RecoveryIdentityTests(unittest.TestCase):
    def test_spe_retry_is_a_complex_with_canonical_molecule(self):
        from component_store import identity
        result=identity('/pscratch/example/completion_runs/campaign/spe/C2H6_Ir100/r2scan')
        self.assertEqual(result['role'],'complex')
        self.assertEqual(result['surface'],'Ir100')
        self.assertEqual(result['molecule'],'ethane')
