import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write

from component_store import identity, record, merge
from prepare_completion_jobs import sha
from prepare_relaxation_restarts import prepare, restart_settings, validate_coordinates
from report_relaxation_restarts import report


class RestartTests(unittest.TestCase):
    def test_optimizer_and_scf_changes_preserve_physics(self):
        original = dict(EDIFF='1E-05',EDIFFG='-0.05',IBRION='2',IALGO='48',
                        NELM='150',GGA='BF',ENCUT='450',LUSE_VDW='.TRUE.')
        changed = restart_settings(original,True,True)
        self.assertNotIn('IALGO',changed)
        self.assertEqual(changed['ALGO'],'Normal')
        self.assertEqual(changed['IBRION'],'1')
        self.assertEqual(changed['NELM'],'300')
        self.assertLessEqual(float(changed['EDIFF']),1e-6)
        for key in ['EDIFFG','GGA','ENCUT','LUSE_VDW']:
            self.assertEqual(changed[key],original[key])
        untouched=restart_settings(original)
        for key,value in original.items():self.assertEqual(untouched[key],value)

    def fixture(self, root):
        source=root/'dft_jobs'/'CO_Ag100'/'PBE';source.mkdir(parents=True)
        atoms=Atoms('AgCO',positions=[[1,1,1],[1,1,3],[1,1,4.2]],cell=[8,8,20],pbc=True)
        atoms.set_constraint(FixAtoms(indices=[0]))
        write(source/'POSCAR',atoms,format='vasp')
        atoms.positions[1,2]+=0.1
        write(source/'CONTCAR',atoms,format='vasp')
        (source/'INCAR').write_text('NSW=1000\nIBRION=2\nEDIFF=1E-05\nEDIFFG=-0.05\nNELM=150\nIALGO=48\nGGA=PE\n')
        (source/'POTCAR').write_text('test fixture only; not a potential')
        (source/'KPOINTS').write_text('test k points')
        (source/'OUTCAR').write_text('Interrupted test output')
        row=dict(directory=str(source),surface='Ag100',molecule='CO',functional='pbe',
                 status='unfinished',full_output_sha256=sha(source/'OUTCAR'),
                 normal_completion_anywhere=False,ionic_convergence_anywhere=False,
                 zbrent_fatal=False,scf_limit_warning=False,cause='Slurm time limit')
        audit=root/'audit.json';audit.write_text(json.dumps([row]))
        result=dict(directory=str(source),status='unfinished',functional='pbe',nsw=1000,settings={'IBRION':'2'})
        return source,audit,result

    def test_staging_preserves_sources_constraints_and_potentials(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source,audit,result=self.fixture(root)
            before={p.name:sha(p) for p in source.iterdir()}
            campaign=root/'completion_runs'/'trial'
            with patch('prepare_relaxation_restarts.component',return_value=result):
                m=prepare(root,campaign,audit,root/'absent.gz')
                with self.assertRaises(ValueError):prepare(root,campaign,audit,root/'absent.gz')
            job=m['jobs'][0];dest=Path(job['directory'])
            self.assertEqual(before,{p.name:sha(p) for p in source.iterdir()})
            self.assertEqual(sha(dest/'POSCAR'),sha(source/'CONTCAR'))
            self.assertEqual(sha(dest/'POTCAR'),sha(source/'POTCAR'))
            self.assertEqual(identity(dest)['role'],'complex')
            self.assertEqual(identity(dest)['molecule'],'CO')
            self.assertNotIn('OUTCAR',[p.name for p in dest.iterdir()])

    def test_changed_output_or_completed_alternative_prevents_staging(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source,audit,result=self.fixture(root)
            campaign=root/'completion_runs'/'trial';store=root/'store.gz'
            completed=dict(result,status='converged',directory='/scratch/other/dft_jobs/CO_Ag100/PBE')
            merge(store,[record('kestrel',completed,{})])
            with self.assertRaisesRegex(ValueError,'Completed alternatives'):prepare(root,campaign,audit,store)
            self.assertFalse(campaign.exists())
            (source/'OUTCAR').write_text('A newer calculation')
            with patch('prepare_relaxation_restarts.component',return_value=result):
                with self.assertRaisesRegex(ValueError,'changed since audit'):prepare(root,campaign,audit,root/'absent.gz')
            self.assertFalse(campaign.exists())

    def test_changed_constraints_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            source,_,_=self.fixture(Path(temp))
            from ase.io import read
            atoms=read(source/'CONTCAR');atoms.set_constraint()
            write(source/'CONTCAR',atoms,format='vasp')
            with self.assertRaisesRegex(ValueError,'Constraints changed'):validate_coordinates(source)

    def test_scheduler_completion_does_not_imply_ionic_convergence(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);campaign=root/'campaign';campaign.mkdir()
            audit=campaign/'source_audit.json'
            audit.write_text(json.dumps([dict(cause='Slurm time limit',scf_limit_warning=False)]*3))
            jobs=[dict(array_index=i,directory='/restart/'+str(i),source_directory='/original/'+str(i),
                       source_sha256={'CONTCAR':'hash'},surface='Ag100',molecule='CO',functional='pbe',
                       cause='Slurm time limit',changes={}) for i in range(3)]
            (campaign/'manifest.json').write_text(json.dumps(dict(audit_sha256=sha(audit),jobs=jobs)))
            (campaign/'submission.json').write_text(json.dumps(dict(job_id='123',submitted_at='2026-09-23')))
            (campaign/'cross_cluster_review.json').write_text(json.dumps(dict(live_verified=False)))
            for name in ['run_array.slurm','joblist.txt']:(campaign/name).write_text('fixture')
            accounting='123_1|COMPLETED|0:0|00:01:00|12:00:00|start|end\n123_2|COMPLETED|0:0|00:01:00|12:00:00|start|end\n'
            calculations=[dict(status='invalid',nsw=1000),dict(status='unfinished',nsw=1000),dict(status='converged',nsw=1000)]
            with patch('report_relaxation_restarts.subprocess.check_output',side_effect=[accounting,'123_0|PENDING\n']),patch('report_relaxation_restarts.component',side_effect=calculations):
                result=report(campaign,root)
            self.assertEqual(result['summary']['statuses'],{'queued':1,'needs review':1,'converged':1})


if __name__=='__main__':unittest.main()
