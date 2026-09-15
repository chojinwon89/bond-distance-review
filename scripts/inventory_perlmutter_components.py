#!/usr/bin/env python3
"""Inventory every prepared calculation in the three requested Perlmutter trees."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from extract_kestrel_singlepoint import component


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,default=Path(__file__).resolve().parents[1])
    args=p.parse_args();project=args.project_root.resolve();rows=[];queue={}
    try:
        tasks=(project/'perlmutter/spe_runs/ml_poscar_20260910/joblist.txt').read_text().splitlines()
        output=subprocess.check_output(['squeue','-r','-j','58164924','-h','-o','%K|%T'],text=True)
        for line in output.splitlines():
            task,state=line.split('|')
            if task.isdigit() and int(task)<len(tasks):queue[tasks[int(task)]]=(state,'58164924_'+task)
    except (OSError,subprocess.CalledProcessError):pass
    for tree in ['vasp_mol','vasp_slab','dft_jobs']:
        root=project/tree
        paths=subprocess.check_output(['rg','--files','--hidden','--no-ignore',str(root),'-g','INCAR'],text=True).splitlines()
        for path in sorted(paths):
            directory=Path(path).parent
            # Do not follow a calculation link outside the requested data tree.
            if root not in directory.resolve().parents:continue
            r=component(directory);state,job=queue.get(str(directory),('',''))
            mode='single-point' if r['nsw']==0 else 'relaxation' if r['nsw'] and r['nsw']>0 else 'unknown'
            rows.append(dict(tree=tree,system=directory.relative_to(root).parts[0],directory=str(directory),
                functional=r['functional'],mode=mode,NSW=r['nsw'],status=r['status'],
                E_TOTEN_eV=r['energy'],queue_state=state,queue_job=job,
                composition=json.dumps(r['composition'],sort_keys=True),potentials=json.dumps(r['potentials'],sort_keys=True),
                OUTCAR=r['outcar'],outcar_mtime_ns=r.get('outcar_mtime_ns',''),note=r.get('note','')))
        print('Inventoried',tree,len(paths),flush=True)
    with (args.output_root/'dft_perlmutter_components.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(extracted_at=datetime.now(timezone.utc).isoformat(),project_root=str(project),
        policy='All INCAR paths discovered without ignore rules or recursive symlink following. NSW=0 identifies SPE; NSW>0 identifies relaxation. A raw TOTEN is not an adsorption energy. Nonconverged totals are recorded as diagnostics only.',
        calculations=len(rows),by_tree=dict(Counter(r['tree'] for r in rows)),
        by_mode=dict(Counter(r['mode'] for r in rows)),by_status=dict(Counter(r['status'] for r in rows)),
        converged_by_tree_and_mode=dict(Counter(r['tree']+'/'+r['mode'] for r in rows if r['status']=='converged')))
    (args.output_root/'dft_perlmutter_components_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
