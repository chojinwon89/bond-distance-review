"""Audit the pre-policy candidate set; retain mismatched energies as diagnostics."""
import csv,json
from collections import Counter
from pathlib import Path
from component_store import read_store,latest,reference_candidates
from potential_matching import potential_match

ROOT=Path(__file__).resolve().parents[1]

def audit(root=ROOT):
    records=read_store(root/'dft_component_store.jsonl.gz');by={r['snapshot_id']:r for r in records};current=latest(records)
    with (root/'dft_potential_candidates_before.csv').open() as f:before=list(csv.DictReader(f))
    rows=[]
    for entry in before:
        c,s,g=[by[entry[k+'_snapshot_id']] for k in ['complex','slab','molecule']]
        bad=[role for role,r in [('slab',s),('molecule',g)] if not potential_match(c['calculation'],r['calculation'])]
        if not bad:continue
        refs=reference_candidates(c,current)['candidate_reference_ids']
        rows.append(dict(entry,mismatched_roles=','.join(bad),complex_potentials=json.dumps(c['calculation']['potentials'],sort_keys=True),
            slab_potentials=json.dumps(s['calculation']['potentials'],sort_keys=True),molecule_potentials=json.dumps(g['calculation']['potentials'],sort_keys=True),
            compatible_slab_count=len(refs['slab']),compatible_molecule_count=len(refs['molecule']),
            compatible_slab_ids=json.dumps(refs['slab'])))
    with (root/'dft_potential_mismatches.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=dict(policy='Matching PAW TITEL identities (element, variant and dataset date) required for slab and gas references. No scaling or offsets. Existing mismatched values remain only in historical audit downloads.',
        audited_candidates=len(before),mismatched_candidates=len(rows),by_surface=dict(sorted(Counter(r['surface'] for r in rows).items())),
        mismatch_cells=len({(r['surface'],r['molecule'],r['functional'],r['mode']) for r in rows}),
        candidates_with_compatible_slab=sum(int(r['compatible_slab_count'])>0 for r in rows),
        provenance='Live Perlmutter scan plus saved Kestrel exports; fresh Kestrel authentication unavailable.')
    with (root/'dft_perlmutter_singlepoint_audit.csv').open() as f:fresh=list(csv.DictReader(f))
    summary['current_perlmutter_SPE_scan']=dict(candidates=len(fresh),statuses=dict(Counter(r['status'] for r in fresh)),source='dft_perlmutter_singlepoint_audit.csv')
    (root/'dft_potential_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':audit()
