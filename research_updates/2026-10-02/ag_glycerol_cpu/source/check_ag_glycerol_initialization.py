"""Verify the contact intervention without MLIP energy evaluations."""
import json
import random
import numpy as np
from ase.calculators.emt import EMT
from run_ag_glycerol_cpu_contact import OUT,SEEDS,InitialContactGOAD,load_inputs,save,contacts
from run_oxygenated_pilot import MeteredGOAD

slab,mol=load_inputs()
kwargs=dict(surface=slab,molecule=mol,calculator=EMT(),surface_energy=0,molecule_energy=0,
            calculator_type=None,n_workers=1,n_fixed_layers=2,population_size=60,generations=200,
            mutation_rate=.3,crossover_rate=.7,elite_size=5,verbose=False)
records=[]
for seed in SEEDS:
    runs={}
    for mode in ['original','baseline','contact25']:
        np.random.seed(seed);random.seed(seed)
        ga=MeteredGOAD(**kwargs) if mode=='original' else InitialContactGOAD(**kwargs,initial_mode=mode)
        ga._initialize_population()
        states=(np.random.get_state(),random.getstate())
        runs[mode]=(ga,[ga._create_system(g) for g in ga.population],states)
    for original,baseline,contact in zip(runs['original'][1],runs['baseline'][1],runs['contact25'][1]):
        assert np.array_equal(original.positions,baseline.positions)
        assert np.allclose(baseline.positions[:len(slab)],contact.positions[:len(slab)])
        delta=contact.positions[len(slab):]-baseline.positions[len(slab):]
        assert np.allclose(delta[:,:2],0,atol=1e-9)
        assert np.allclose(delta[:,2],delta[0,2],atol=1e-9)
        assert abs(contacts(contact,len(slab))['nearest_O_Ag_A']-2.5)<1e-7
    for mode in ['baseline','contact25']:
        assert np.array_equal(runs['original'][2][0][1],runs[mode][2][0][1])
        assert runs['original'][2][0][2:]==runs[mode][2][0][2:]
        assert runs['original'][2][1]==runs[mode][2][1]
    counts={i:sum(contacts(a,len(slab))['closest_oxygen_index']==i for a in runs['contact25'][1]) for i in [0,3,5]}
    assert all(counts.values()),counts
    records.append(dict(seed=seed,initial_closest_oxygen_counts=counts,baseline_unchanged=True,random_state_unchanged=True))
save(OUT/'initialization_checks.json',dict(passed=True,seeds=records))
print(json.dumps(records),flush=True)
