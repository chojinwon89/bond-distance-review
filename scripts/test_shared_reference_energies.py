import copy
from decimal import Decimal
import unittest

from component_store import record
from shared_reference_energies import derive, fill_relaxed


def components():
    settings=dict(ENCUT='450',ISPIN='2',ISMEAR='1',SIGMA='0.05',IBRION='2')
    cell=[[8.,0,0],[4.,7.,0],[0,0,37.]]
    def make(cluster,path,energy,composition):
        return record(cluster,dict(directory=path,energy=energy,composition=composition,cell=copy.deepcopy(cell),
            status='converged',functional='beef_vdw',nsw=100,settings=dict(settings),potentials={}),{'kind':'archived-source-export'})
    return [make('perlmutter','/data/dft_jobs/CH3OCH3_Ag111/beef_vdw',-10.65320050,{'Ag':36,'C':2,'O':1,'H':6}),
            make('kestrel','/data/vasp_slab/Ag111_n36/beef_vdw',56.86657689,{'Ag':36}),
            make('perlmutter','/data/vasp_mol/DME/beef_vdw',-42.29312929,{'C':2,'O':1,'H':6})]


class SharedReferenceTests(unittest.TestCase):
    def test_unusual_value_and_positive_slab_total_are_preserved(self):
        result=derive(*components(),'relaxed')
        self.assertEqual(Decimal(result['E_ads']),Decimal('-25.22664810'))
        self.assertEqual(result['status'],'energy_review')
        self.assertEqual(result['slab_cluster'],'kestrel')

    def test_invalid_reference_cannot_be_promoted(self):
        for kind in ['atom_count','cell','unfinished','isomer','cutoff','mode']:
            c,s,g=copy.deepcopy(components())
            if kind=='atom_count':s['calculation']['composition']={'Ag':64}
            if kind=='cell':s['calculation']['cell'][0][0]=9.
            if kind=='unfinished':g['calculation']['status']='unfinished'
            if kind=='isomer':g['molecule']='ethanol'
            if kind=='cutoff':s['calculation']['settings']=dict(s['calculation']['settings'],ENCUT='500')
            with self.subTest(kind=kind),self.assertRaises(ValueError):
                derive(c,s,g,'SPE' if kind=='mode' else 'relaxed')

    def test_page_fill_marks_review_preserves_existing_and_is_idempotent(self):
        data=derive(*components(),'relaxed')
        page='<div class="g" data-surf="Ag111" data-mol="DME"><table><tr><td class="l">BEEF-vdW</td><td>—</td><td>-0.20</td><td>—</td></tr></table></div><script></script>'
        filled=fill_relaxed(page,[data]);self.assertIn('energy-review',filled);self.assertIn('>-25.227</td>',filled)
        self.assertEqual(fill_relaxed(filled,[data]),filled)
        published=page.replace('<td>—</td>','<td>-1.234</td>',1)
        self.assertEqual(fill_relaxed(published,[data]),published)


if __name__=='__main__':unittest.main()
