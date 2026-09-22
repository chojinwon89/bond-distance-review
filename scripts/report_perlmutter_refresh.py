"""Summarize a website refresh against a pinned public Git revision."""
import csv
import gzip
import io
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from component_store import latest, read_store
from potential_matching import potential_match

ROOT=Path(__file__).resolve().parents[1]
BEFORE='23048650b883c5118871d1829a317187a65f0989'


def previous(name):
    return subprocess.check_output(['git','show',BEFORE+':'+name],cwd=ROOT)


def csv_rows(path):
    with path.open() as f:return list(csv.DictReader(f))


def number(value):
    try:return float(value)
    except (ValueError,TypeError):return None


def main():
    before=list(csv.DictReader(io.StringIO(previous('dft_completion_coverage.csv').decode())))
    after=csv_rows(ROOT/'dft_completion_coverage.csv')
    key=lambda r:(r['surface'],r['molecule'],r['functional'])
    old={key(r):r for r in before};changes=[];counts=Counter()
    for row in after:
        for mode in ['relaxed','SPE']:
            a=number(old.get(key(row),{}).get('E_ads_'+mode));b=number(row['E_ads_'+mode])
            state='unchanged'
            if a is None and b is not None:state='newly_filled'
            elif a is not None and b is None:state='withheld_unverified_or_incompatible'
            elif a is not None and b is not None and abs(a-b)>0.00051:state='updated'
            counts[mode+'_'+state]+=1
            if state!='unchanged':changes.append(dict(surface=row['surface'],molecule=row['molecule'],functional=row['functional'],mode=mode,previous_eV=a,current_eV=b,change=state,current_status=row[mode+'_status']))
    with (ROOT/'dft_refresh_20260922_changes.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(changes[0]));w.writeheader();w.writerows(changes)
    records=read_store(ROOT/'dft_component_store.jsonl.gz');byid={r['snapshot_id']:r for r in records};current=latest(records)
    earlier=latest([json.loads(line) for line in gzip.decompress(previous('dft_component_store.jsonl.gz')).decode().splitlines()])
    old_calcs={(r['cluster'],r['calculation']['directory']):r for r in earlier}
    transitions=[]
    for r in current:
        if r['cluster']!='perlmutter':continue
        c=r['calculation'];prev=old_calcs.get((r['cluster'],c['directory']))
        if c['status']=='converged' and (prev is None or prev['calculation']['status']!='converged'):
            transitions.append(dict(role=r['role'],directory=c['directory'],functional=c['functional'],previous_status=prev['calculation']['status'] if prev else 'not in previous store',status='converged',snapshot_id=r['snapshot_id']))
    selections=csv_rows(ROOT/'dft_cluster_selection.csv');verified=0
    for r in selections:
        if r.get('application_status')!='selected candidate applied':continue
        c,s,g=[byid[r[role+'_snapshot_id']]['calculation'] for role in ['complex','slab','molecule']]
        assert all(x['status']=='converged' for x in [c,s,g])
        assert potential_match(c,s) and potential_match(c,g)
        assert abs(float(r['E_ads'])-(c['energy']-s['energy']-g['energy']))<1e-7
        verified+=1
    summary=dict(refreshed_at=datetime.now(timezone.utc).isoformat(),previous_public_commit=BEFORE,
        scope='Fresh Perlmutter files; archived Kestrel records retained with provenance. No new DFT jobs submitted.',
        page_changes=dict(counts),numeric_cells_before={m:sum(number(r['E_ads_'+m]) is not None for r in before) for m in ['relaxed','SPE']},
        numeric_cells_after={m:sum(number(r['E_ads_'+m]) is not None for r in after) for m in ['relaxed','SPE']},
        newly_available_converged_components_by_role=dict(Counter(r['role'] for r in transitions)),
        previously_recorded_now_converged_by_role=dict(Counter(r['role'] for r in transitions if r['previous_status']!='not in previous store')),
        newly_available_converged_components=transitions,
        verified_selected_component_sets=verified,all_selected_energies_subtraction_and_PAW_identity_verified=True,
        matched_potential_recovery_status=dict(Counter(r['calculation']['status'] for r in current if '/completion_runs/matched_potential_slabs_20260915/' in r['calculation']['directory'])),
        dme_gas_recovery=json.loads((ROOT/'dft_gas_reference_recovery.json').read_text())['jobs'],
        policy='Perlmutter first; matching potentials, functional, composition and full slab cell; no scaling or offsets. Unusual valid energies retain review labels. Historical unverified/mixed-potential values are not promoted as validated data.')
    (ROOT/'dft_refresh_20260922.json').write_text(json.dumps(summary,indent=2)+'\n')
    text=(f"<b>Perlmutter refresh · 22 September 2026:</b> {counts['relaxed_newly_filled']} relaxed and {counts['SPE_newly_filled']} SPE blanks filled. "
          f"Current validated tables contain {summary['numeric_cells_after']['relaxed']} relaxed and {summary['numeric_cells_after']['SPE']} SPE energies, including review-marked values. "
          '<a href="dft_refresh_20260922_changes.csv">Cell-by-cell changes</a> · <a href="dft_refresh_20260922.json">Convergence and reference audit</a>.')
    import re
    for name in ['index.html','dft_comparison.html']:
        page=(ROOT/name).read_text()
        page=re.sub(r'<!-- refresh-20260922 -->.*?<!-- /refresh-20260922 -->\n?', '',page,flags=re.S)
        note='<!-- refresh-20260922 -->\n<div class="sub" style="margin:12px 0;padding:10px;border:1px solid #60758a;border-radius:8px">'+text+'</div>\n<!-- /refresh-20260922 -->\n'
        marker='<div id="count">' if name=='index.html' else '<h2>Per-system structure'
        assert marker in page
        (ROOT/name).write_text(page.replace(marker,note+marker,1))
    print(json.dumps({k:v for k,v in summary.items() if k not in ['newly_available_converged_components','dme_gas_recovery']},indent=2))


if __name__=='__main__':main()
