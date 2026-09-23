import unittest
import numpy as np
from ase import Atoms
from functional_geometry import measure, correspondence
from component_store import reference_candidates

class GeometryTests(unittest.TestCase):
    def test_periodic_surface_contact_uses_nearest_image(self):
        a=Atoms('AgCO',positions=[[0.1,1,0],[9.9,1,2],[9.9,1,3.2]],cell=[10,10,20],pbc=[True,True,False])
        result=measure(a,'Ag')
        self.assertAlmostEqual(result['nearest_contact_A'],np.sqrt(4.04),places=6)
        self.assertEqual(result['contacts'][0]['label'],'C1')

    def test_heavy_atom_mapping_handles_permutation(self):
        a=Atoms('AgCCO',positions=[[0,0,0],[2,2,2],[3.5,2,2],[4.8,2.3,2]],cell=[15,15,20],pbc=True)
        b=a[[0,3,1,2]]
        mapping,note=correspondence(a,b,'Ag')
        self.assertEqual(mapping,{1:2,2:3,3:1})

    def test_dissociation_is_not_paired_as_same_connectivity(self):
        a=Atoms('AgCO',positions=[[0,0,0],[2,2,2],[3.2,2,2]],cell=[15,15,20],pbc=True)
        b=a.copy();b.positions[2]=[7,7,7]
        self.assertEqual(correspondence(a,b,'Ag')[0],{})

    def test_linear_molecule_has_angle_but_no_torsion(self):
        a=Atoms('AgOCO',positions=[[0,0,0],[1,1,2],[2.2,1,2],[3.4,1,2]],cell=[15,15,20],pbc=True)
        m=measure(a,'Ag');self.assertEqual(len(m['angles']),1)
        self.assertAlmostEqual(m['angles'][0]['degrees'],180)
        self.assertFalse(m['torsions'])

    def test_single_atom_has_no_invented_angles(self):
        a=Atoms('AgH',positions=[[0,0,0],[0,0,2]],cell=[10,10,20],pbc=True)
        m=measure(a,'Ag');self.assertFalse(m['angles']);self.assertFalse(m['torsions'])

if __name__=='__main__':unittest.main()
