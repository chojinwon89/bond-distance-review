import tempfile
import unittest
from pathlib import Path

from component_store import identity, latest, merge, read_store, record, reference_candidates


def fixture(cluster='perlmutter', path='/data/dft_jobs/CH3OCH3_Ag111/PBE', energy=-138., mtime=1):
    return record(cluster, dict(directory=path, energy=energy, status='converged', functional='pbe',
        nsw=100, outcar_mtime_ns=mtime, composition={'Ag':36,'C':2,'O':1,'H':6},
        settings={'IBRION':'2'}, cell=[[8.,0,0],[4.,7.,0],[0,0,37.]]), {'kind':'live-files','observed_at':'2026-09-15'})


class StoreTests(unittest.TestCase):
    def test_aliases_and_nested_layouts(self):
        self.assertEqual(identity('/data/poscar/best/C2/Ag111_CH3OCH3/beef_vdw/fully_relaxed')['molecule'], 'DME')
        self.assertEqual(identity('/data/vasp_slab/Ag111_n36/PBE')['surface'], 'Ag111')
        self.assertNotEqual(identity('/data/vasp_mol/C2H5OH/PBE')['molecule'], 'DME')
        self.assertNotEqual(identity('/data/vasp_mol/C2H6O/PBE')['molecule'], 'DME')

    def test_reimport_preserves_both_clusters_and_versions(self):
        with tempfile.TemporaryDirectory() as temp:
            store = Path(temp)/'store.gz'
            old = fixture(); other = fixture(cluster='kestrel'); new = fixture(energy=-139.,mtime=2)
            merge(store,[old,other]); merge(store,[new,old]); before=store.read_bytes()
            merge(store,[old,new]); self.assertEqual(before,store.read_bytes())
            self.assertEqual(len(read_store(store)),3)
            current=latest(read_store(store)); self.assertEqual(len(current),2)
            self.assertEqual(next(r for r in current if r['cluster']=='perlmutter')['calculation']['energy'],-139.)
            merge(store,[]); self.assertEqual(before,store.read_bytes())

    def test_reference_candidates_require_identity_composition_cell_and_convergence(self):
        comp=fixture()
        slab=fixture('kestrel','/data/vasp_slab/Ag111_n36/PBE')
        slab['calculation']['composition']={'Ag':36}
        gas=fixture('kestrel','/data/vasp_mol/DME/PBE')
        gas['calculation']['composition']={'C':2,'O':1,'H':6}
        result=reference_candidates(comp,[slab,gas])
        self.assertEqual(result['unresolved_roles'],[])
        gas['calculation']['status']='unfinished'
        self.assertEqual(reference_candidates(comp,[slab,gas])['unresolved_roles'],['molecule'])
        slab['calculation']['composition']={'Ag':64}
        self.assertEqual(reference_candidates(comp,[slab,gas])['unresolved_roles'],['slab','molecule'])
        slab['calculation']['composition']={'Ag':36}; slab['calculation']['cell'][0][0]=9.
        self.assertEqual(reference_candidates(comp,[slab])['unresolved_roles'],['slab','molecule'])

    def test_new_failure_is_not_hidden_by_old_success(self):
        old=fixture(); new=fixture(mtime=2); new['calculation']['status']='unfinished'
        self.assertEqual(latest([new,old])[0]['calculation']['status'],'unfinished')

    def test_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'store.gz'; row=fixture(); merge(path,[row]); before=path.read_bytes()
            row['calculation']['energy']=99.
            with self.assertRaisesRegex(ValueError,'checksum'):
                merge(path,[row])
            self.assertEqual(before,path.read_bytes())


if __name__=='__main__':unittest.main()
