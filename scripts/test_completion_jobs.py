import json
import hashlib
from pathlib import Path
import tempfile
import unittest
from prepare_completion_jobs import incar_text, metal_potential
from completion_support import completion_jobs, matching_slab


class CompletionTests(unittest.TestCase):
    def test_supplemental_slab_rejects_wrong_target_cell_and_changed_input(self):
        from ase import Atoms
        from ase.io import write
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);slab=root/'slab';comp=root/'complex';slab.mkdir();comp.mkdir()
            atoms=Atoms('Cu',positions=[[0,0,0]],cell=[5,5,15],pbc=True)
            write(slab/'POSCAR',atoms,format='vasp');write(comp/'POSCAR',atoms,format='vasp')
            job=dict(role='slab',functional='beef_vdw',target_systems=['C2H4_Cu100'],directory=str(slab),
                     staged_sha256={'POSCAR':hashlib.sha256((slab/'POSCAR').read_bytes()).hexdigest()})
            self.assertTrue(matching_slab(job,'C2H4_Cu100','beef_vdw',comp))
            self.assertFalse(matching_slab(job,'C2H6_Cu100','beef_vdw',comp))
            self.assertFalse(matching_slab(job,'C2H4_Cu100','r2scan',comp))
            atoms.cell[0,0]=6;write(comp/'POSCAR',atoms,format='vasp')
            self.assertFalse(matching_slab(job,'C2H4_Cu100','beef_vdw',comp))
            write(slab/'POSCAR',atoms,format='vasp')
            self.assertFalse(matching_slab(job,'C2H4_Cu100','beef_vdw',comp))

    def test_potential_subset_exact_element(self):
        text='TITEL = PAW_PBE Cu 2001\nVRHFIN =Cu: d10\nEnd of Dataset\nTITEL = PAW_PBE C 2001\nVRHFIN =C: s2\nEnd of Dataset\n'
        metal=metal_potential(text,'Cu')
        self.assertIn('VRHFIN =Cu:',metal)
        self.assertNotIn('VRHFIN =C:',metal)
        with self.assertRaises(ValueError):metal_potential(text,'Ir')

    def test_incar_retry_retains_electronic_tolerance_and_functional(self):
        text='SYSTEM = C2H6_Ir100\nNSW = 0\nIBRION = -1\nEDIFF = 1E-05\nNELM = 150\nMETAGGA = R2SCAN\n'
        changed=incar_text(text,{'NELM':'300'})
        self.assertIn('EDIFF = 1E-05',changed)
        self.assertIn('METAGGA = R2SCAN',changed)
        self.assertIn('NSW = 0',changed)
        self.assertIn('NELM = 300',changed)
        self.assertNotIn('NELM = 150',changed)

    def test_campaign_path_must_stay_inside_campaign(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);campaign=root/'completion_runs'/'test';campaign.mkdir(parents=True)
            p=campaign/'manifest.json'
            p.write_text(json.dumps({'jobs':[{'directory':str(campaign/'job')}]}))
            self.assertEqual(len(completion_jobs(root)),1)
            p.write_text(json.dumps({'jobs':[{'directory':'/outside/campaign'}]}))
            with self.assertRaises(ValueError):completion_jobs(root)


if __name__=='__main__':unittest.main()
