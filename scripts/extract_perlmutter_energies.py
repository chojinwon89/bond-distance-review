#!/usr/bin/env python3
"""Uses the VASP project collector, then requires final electronic/ionic convergence for all three calculations.
Energy convention: final free-energy TOTEN (eV), matching the existing collector.
"""
import csv,re,functools,collections,json
from pathlib import Path
import argparse, importlib.util
parser = argparse.ArgumentParser(description="Extract and audit Perlmutter adsorption energies without changing calculation files")
parser.add_argument('--vasp-root', type=Path, required=True)
parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'dft_perlmutter_energy_audit.csv')
args = parser.parse_args()
BASE = args.vasp_root.resolve()
spec = importlib.util.spec_from_file_location('collector', BASE / 'calc_binding_energy.py')
collector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collector)
collector.read_energy_from_outcar = functools.lru_cache(None)(collector.read_energy_from_outcar)
collector.count_element_atoms = functools.lru_cache(None)(collector.count_element_atoms)
@functools.lru_cache(None)
def check(rel):
 p=BASE/rel
 if not p.exists(): return 'missing'
 with p.open('rb') as f:
  head=f.read(100000).decode(errors='replace')
  f.seek(max(0,p.stat().st_size-1000000));tail=f.read().decode(errors='replace')
 nsw=re.search(r'\bNSW\s*=\s*(\d+)',head)
 if 'General timing and accounting' not in tail: return 'unfinished'
 if not nsw: return 'unknown_NSW'
 energies=list(re.finditer(r'free  energy   TOTEN',tail))
 if not energies: return 'no_final_energy'
 start=energies[-2].end() if len(energies)>1 else 0
 if 'aborting loop because EDIFF is reached' not in tail[start:energies[-1].start()]: return 'electronic_unconverged'
 if int(nsw.group(1))>0 and 'reached required accuracy' not in tail: return 'ionic_unconverged'
 return 'converged'
rows = []
for functional in ['PBE', 'PBE_D3', 'r2scan', 'beef_vdw']:
 rows.extend(collector.calc_binding_energies([BASE/'dft_jobs'], BASE/'vasp_slab', BASE/'vasp_mol', functional=functional, molecule_first=True))
 print('Extracted', functional, flush=True)
for r in rows:
 f={'pbe':'PBE','pbe_d3':'PBE_D3','r2scan':'r2scan','beef_vdw':'beef_vdw'}[r['functional']]
 mol={'C2H5OH':'CH3CH2OH'}.get(r['molecule'],r['molecule'])
 paths=[f'dft_jobs/{r["system"]}/{f}/OUTCAR',f'vasp_slab/{r["surface"]}/{f}/OUTCAR',f'vasp_mol/{mol}/{f}/OUTCAR']
 statuses=[check(p) for p in paths]
 r['convergence']='; '.join(f'{k}: {v}' for k,v in zip(['complex','slab','molecule'],statuses))
 r['publishable']=str(r['status']=='ok' and all(s=='converged' for s in statuses)).lower()
 r['source_dir']='dft_jobs'
 r['note']=(r['note'] or '').replace(str(BASE)+'/', '')
 r['complex_outcar'],r['slab_outcar'],r['molecule_outcar']=paths
with args.output.open('w', newline='') as out:
 w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('collector status',collections.Counter(r['status'] for r in rows))
print('publishable by functional',collections.Counter(r['functional'] for r in rows if r['publishable']=='true'))
print('unpublishable numerical results',collections.Counter(r['convergence'] for r in rows if r['status']=='ok' and r['publishable']!='true'))
