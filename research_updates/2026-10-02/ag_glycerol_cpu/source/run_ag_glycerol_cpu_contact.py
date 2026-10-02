"""Matched GOAD CPU experiment: original versus nearest O--Ag 2.5 A initialization."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import sys
import time
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from run_oxygenated_pilot import Meter, MeteredGOAD, StopGA, bottom_two, relax, save, sha
from relaxed_endpoint_metrics import endpoint_metrics
from cgas.geometry import internal_distances

OUT = ROOT/'benchmarks/ag_glycerol_cpu_contact_20261002'
SEEDS = [101,202,303,404,505,606,707,808,909,1010]

def contacts(atoms, n):
    oi = np.flatnonzero(atoms.numbers[n:]==8)+n
    delta = atoms.positions[oi,None,:]-atoms.positions[None,:n,:]
    _,ds = find_mic(delta.reshape(-1,3),atoms.cell,[True,True,False])
    ds=ds.reshape(len(oi),n).min(axis=1)
    return dict(oxygen_indices=(oi-n).tolist(),O_Ag_distances_A=ds.tolist(),
                nearest_O_Ag_A=float(ds.min()),closest_oxygen_index=int(oi[np.argmin(ds)]-n))

def project(atoms,n,target=2.5):
    """Vertical translation giving min over ALL O--Ag pairs exactly target.

    Taking the highest upper-sphere intersection keeps every oxygen outside
    the Ag-centered target spheres and avoids burying other O atoms to bring
    a selected upper oxygen down. Molecular conformation/orientation unchanged.
    """
    assert np.all(atoms.numbers[:n]==47)
    assert np.allclose(atoms.cell.array[:2,2],0,atol=1e-8)
    oi=np.flatnonzero(atoms.numbers[n:]==8)+n
    delta=atoms.positions[oi,None,:]-atoms.positions[None,:n,:]
    delta[:,:,2]=0
    _,lateral=find_mic(delta.reshape(-1,3),atoms.cell,[True,True,False])
    lateral=lateral.reshape(len(oi),n)
    reachable=lateral<=target
    shifts=atoms.positions[None,:n,2]+np.sqrt(np.maximum(0,target**2-lateral**2))-atoms.positions[oi,None,2]
    shift=float(shifts[reachable].max())
    atoms.positions[n:,2]+=shift
    assert abs(contacts(atoms,n)['nearest_O_Ag_A']-target)<1e-7
    return shift

class InitialContactGOAD(MeteredGOAD):
    def __init__(self,*args,initial_mode,**kwargs):
        self.initial_mode=initial_mode
        self.initial_records=[]
        self.initial_structures=[]
        super().__init__(*args,**kwargs)

    def _initialize_population(self):
        super()._initialize_population()
        for i,individual in enumerate(self.population):
            a=super()._create_system(individual)
            before=contacts(a,len(self.surface))
            shift=project(a,len(self.surface)) if self.initial_mode=='contact25' else 0.0
            individual['position'][2]+=shift
            recreated=super()._create_system(individual)
            assert np.allclose(a.positions,recreated.positions,atol=1e-9)
            self.initial_structures.append(recreated)
            self.initial_records.append(dict(candidate=i,before=before,after=contacts(recreated,len(self.surface)),vertical_shift_A=shift))
        self.initial_population_hash=hashlib.sha256(b''.join(a.positions.tobytes() for a in self.initial_structures)).hexdigest()

def load_inputs():
    return read(OUT/'inputs/clean_slab.extxyz'),read(OUT/'inputs/gas_molecule.extxyz')

def calculator():
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID'), 'Compute-node srun required'
    import torch
    from sevenn.calculator import SevenNetCalculator
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    model=SevenNetCalculator('7net-mf-ompa',modal='omat24',device='cpu')
    digest=hashlib.sha256()
    for k,v in sorted(model.model.state_dict().items()):
        digest.update(k.encode());digest.update(str(v.dtype).encode());digest.update(str(tuple(v.shape)).encode())
        digest.update(v.detach().cpu().contiguous().numpy().tobytes())
    manifest=json.loads((OUT/'manifest.json').read_text())
    assert digest.hexdigest()==manifest['model_state_sha256']
    for p,h in manifest['input_sha256'].items():assert sha(p)==h,p
    return model

def preflight():
    model=calculator();slab,mol=load_inputs()
    old=json.loads((OUT/'inputs/references.json').read_text())
    refs={}
    for label,a in [('slab',slab),('molecule',mol)]:
        a.calc=model
        energy=float(a.get_potential_energy());delta=energy-old[label]['energy_eV']
        assert abs(delta)<.002,(label,delta)
        refs[label]=dict(energy_eV=energy,old_gpu_energy_eV=old[label]['energy_eV'],CPU_minus_GPU_eV=delta)
    refs['reference_energy_eV']=sum(refs[k]['energy_eV'] for k in ['slab','molecule'])
    save(OUT/'cpu_references.json',refs)
    endpoint=read(OUT/'inputs/old_goad_seed202.extxyz');endpoint.calc=model
    started=time.perf_counter();energy=float(endpoint.get_potential_energy());endpoint.get_forces()
    elapsed=time.perf_counter()-started
    old_summary=json.loads((OUT/'inputs/old_summary.json').read_text())
    ads=energy-refs['reference_energy_eV']
    assert abs(ads-old_summary['final_E_ads_eV'])<.002
    result=dict(job_id=os.environ['SLURM_JOB_ID'],hostname=os.uname().nodename,device='cpu',
                old_endpoint_E_ads_CPU_eV=ads,old_endpoint_E_ads_GPU_eV=old_summary['final_E_ads_eV'],
                single_evaluation_seconds=elapsed,model_state_verified=True)
    save(OUT/'preflight.json',result);print(json.dumps(result),flush=True)

def run(task):
    import torch
    seed=SEEDS[task//2];mode=['baseline','contact25'][task%2]
    destination=OUT/'runs'/f'{mode}_seed{seed}'
    destination.mkdir(parents=True,exist_ok=False)
    model=calculator();slab,mol=load_inputs();n=len(slab)
    refs=json.loads((OUT/'cpu_references.json').read_text())
    topology=json.loads((OUT/'inputs/topology.json').read_text())
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
    save(destination/'started.json',dict(mode=mode,seed=seed,job_id=os.environ['SLURM_JOB_ID'],
         hostname=os.uname().nodename,device='cpu',script_sha256=sha(__file__)))
    meter=Meter(model,1500)
    ga=InitialContactGOAD(surface=slab,molecule=mol,calculator=meter,
        surface_energy=refs['slab']['energy_eV'],molecule_energy=refs['molecule']['energy_eV'],
        calculator_type=None,n_workers=1,n_fixed_layers=2,population_size=60,generations=200,
        mutation_rate=.3,crossover_rate=.7,elite_size=5,verbose=False,initial_mode=mode)
    started=time.perf_counter();termination='generation_or_stagnation_limit'
    try:ga.run()
    except StopGA:termination='evaluation_budget'
    assert ga.best_individual is not None
    save(destination/'initial_contacts.json',ga.initial_records)
    write(destination/'initial_population.extxyz',ga.initial_structures)
    candidate=ga.best_individual['structure'].copy();candidate.calc=None
    write(destination/'best.extxyz',candidate)
    final,result=relax(candidate,model,200,bottom_two(slab),destination/'final.extxyz')
    result.update(endpoint_metrics(final,n,topology['bonds'],internal_distances(mol)))
    save(destination/'polish.json',result)
    summary=dict(mode=mode,seed=seed,device='cpu',search_evaluations=meter.evaluations,
        search_E_ads_eV=float(ga.best_energy),search_termination=termination,
        final_E_ads_eV=result['energy_eV']-refs['reference_energy_eV'],
        valid=bool(result['converged'] and result['intact'] and result['admissible']),
        elapsed_seconds=time.perf_counter()-started,initial_population_sha256=ga.initial_population_hash,
        **contacts(final,n),**{k:v for k,v in result.items() if k!='trace'})
    write(destination/'final.cif',final)
    save(destination/'search_trace.json',meter.trace)
    save(destination/'summary.json',summary)
    save(destination/'complete.json',dict(finished=True,valid=summary['valid']))
    print(json.dumps(summary),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preflight',action='store_true');parser.add_argument('--task',type=int)
    args=parser.parse_args()
    if args.preflight:preflight()
    else:run(args.task if args.task is not None else int(os.environ['SLURM_PROCID']))
