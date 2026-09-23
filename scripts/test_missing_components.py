import tempfile,unittest
from pathlib import Path
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write
from prepare_missing_components import bottom_constraints,settings_for,validate_inputs

class MissingComponentTests(unittest.TestCase):
    def test_constraints_follow_layers_not_reordered_indices(self):
        original=Atoms('Ag4C',positions=[[0,0,1],[2,0,1],[0,0,3],[2,0,3],[0,0,5]],cell=[8,8,20],pbc=True)
        original.set_constraint(FixAtoms(indices=[0,1]))
        target=original[[4,2,0,3,1]];target.set_constraint()
        self.assertEqual(bottom_constraints(target,original,'Ag'),2)
        self.assertEqual(set(target.constraints[0].get_indices()),{2,4})
        original.set_constraint(FixAtoms(indices=[2,3]))
        with self.assertRaises(ValueError):bottom_constraints(target,original,'Ag')

    def test_mismatched_surface_is_rejected(self):
        a=Atoms('Ag2',positions=[[0,0,1],[0,0,3]],cell=[8,8,20],pbc=True);a.set_constraint(FixAtoms(indices=[0]))
        b=a.copy();b.cell[0,0]=9
        with self.assertRaises(ValueError):bottom_constraints(b,a,'Ag')

    def test_spe_has_no_ionic_motion_and_functional_is_preserved(self):
        initial=dict(GGA='BF',LUSE_VDW='.TRUE.',AGGAC='0.0000',IALGO='48',NELM='150',EDIFF='1e-5',EDIFFG='-0.05',ENCUT='450')
        spe=settings_for(initial,'beef_vdw','spe')
        self.assertEqual((spe['NSW'],spe['IBRION']),('0','-1'))
        self.assertNotIn('IALGO',spe);self.assertNotIn('EDIFFG',spe)
        for key in ['GGA','LUSE_VDW','AGGAC','ENCUT']:self.assertEqual(spe[key],initial[key])
        relax=settings_for(initial,'beef_vdw','relaxed')
        self.assertEqual(relax['EDIFFG'],initial['EDIFFG'])
        self.assertGreater(int(relax['NSW']),0)
        with self.assertRaises(ValueError):settings_for(initial,'pbe','spe')

    def test_potential_order_is_checked(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp)
            write(d/'POSCAR',Atoms('AgO',positions=[[1,1,1],[1,1,3]],cell=[8,8,20],pbc=True),format='vasp',vasp5=True)
            (d/'INCAR').write_text('GGA=PE\nNSW=0\nIBRION=-1\n')
            (d/'KPOINTS').write_text('Gamma\n0\nGamma\n1 1 1\n0 0 0\n')
            def block(el):return 'TITEL = PAW_PBE '+el+' fixture\nVRHFIN ='+el+':\nEnd of Dataset\n'
            (d/'POTCAR').write_text(block('O')+block('Ag'))
            with self.assertRaisesRegex(ValueError,'ordering'):validate_inputs(d)
            (d/'POTCAR').write_text(block('Ag')+block('O'))
            self.assertEqual(len(validate_inputs(d)),2)

if __name__=='__main__':unittest.main()
