#!/usr/bin/env python3
"""Describe each remaining page gap and its current calculation/source dependency."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import re
import subprocess
from bs4 import BeautifulSoup
from extract_perlmutter_singlepoint import canonical
from completion_support import completion_jobs

ROOT=Path(__file__).resolve().parents[1]
FUNCTIONALS={'PBE':'pbe','PBE+D3':'pbe_d3','r²SCAN':'r2scan','BEEF-vdW':'beef_vdw'}


def numeric(text):
    try:return float(text)
    except ValueError:return None


def rows(path):
    if not path.exists():return []
    with path.open() as f:return list(csv.DictReader(f))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--site-root',type=Path,default=ROOT)
    args=p.parse_args();root=args.site_root;project=args.project_root
    page=(root/'dft_comparison.html').read_text();soup=BeautifulSoup(page,'html.parser')
    spe={}
    for r in rows(root/'dft_perlmutter_singlepoint_audit.csv'):
        key=(r['surface'],r['molecule'],r['functional'])
        spe.setdefault(key,[]).append(r)
    relax={}
    for r in rows(root/'dft_perlmutter_energy_audit.csv'):
        relax.setdefault((r['surface'],canonical(r['molecule']),r['functional']),[]).append(r)
    image_audit=json.loads((root/'dft_structure_sources.json').read_text())
    dft_sources={(r['surface'],r['molecule']):r['sources']['CONTCAR']['path']
                 for r in image_audit['candidates'] if r.get('selected')}
    queued={};states=Counter()
    joblist=project/'perlmutter/spe_runs/ml_poscar_20260910/joblist.txt'
    try:
        tasks=joblist.read_text().splitlines()
        output=subprocess.check_output(['squeue','-r','-j','58164924','-h','-o','%K|%T'],text=True)
        for line in output.splitlines():
            index,state=line.split('|')
            if index.isdigit() and int(index)<len(tasks):queued[tasks[int(index)]]=state;states[state]+=1
    except (OSError,subprocess.CalledProcessError):states['unavailable']=1
    completion=completion_jobs(project);recovery={};campaigns=[]
    for manifest in sorted((project/'completion_runs').glob('*/manifest.json')):
        submission=manifest.parent/'submission.json'
        if submission.exists():
            record=json.loads(submission.read_text());job=record['job_id']
            try:state=subprocess.check_output(['squeue','-j',job,'-h','-o','%T'],text=True).strip() or 'finished; see results'
            except (OSError,subprocess.CalledProcessError):state='unknown'
            campaign=dict(job_id=job,state=state,manifest=str(manifest))
            export_file=manifest.parent/'export_submission.json'
            if export_file.exists():campaign['export']=json.loads(export_file.read_text())
            campaigns.append(campaign)
    for j in completion:
        result_file=Path(j['directory'])/'completion_result.json'
        if result_file.exists():
            result=json.loads(result_file.read_text())
            status='recovery completed; awaiting usable adsorption result' if result['completed'] else 'recovery failed convergence; review outputs'
        else:
            campaign=next((c for c in campaigns if Path(c['manifest']).parent in Path(j['directory']).parents),None)
            status='recovery '+campaign['state'].lower() if campaign else 'recovery prepared; not submitted'
        for system in j['target_systems']:
            m=re.fullmatch(r'(.+)_([A-Z][a-z]?\d+)(?:_.+)?',system)
            if m:recovery[(m.group(2),canonical(m.group(1)),j['functional'])] = dict(directory=j['directory'],status=status)
    report=[];geometry=[]
    for card in soup.select('.g'):
        s,m=card['data-surf'],card['data-mol'];key2=(s,m)
        info=card.select_one('.m').get_text(' ',strip=True)
        dft_match=re.search(r'DFT contact:\s*([A-Z][a-z]?-[A-Z][a-z]?)\s*([\d.]+)',info)
        base='published base' in info
        geometry.append(dict(surface=s,molecule=m,mlip_image=bool(card.select_one('.pair > div:first-child img')),
                             mlip_variant='published base GOAD' if base else 'SevenNet-OMNI',dft_image=bool(card.select_one('.pair > div:nth-child(2) img')),
                             dft_contact_available=bool(dft_match),missing_geometry_source=dft_sources.get(key2,'') if not dft_match else ''))
        for tr in card.select('table tr'):
            cells=tr.select('td')
            if not cells:continue
            f=FUNCTIONALS[cells[0].get_text()];key=(s,m,f)
            er,em,es=[numeric(cells[i].get_text()) for i in [1,2,4]]
            candidates=spe.get(key,[])
            if es is not None:status='available'
            elif key in recovery:status=recovery[key]['status']
            elif any(r['complex_directory'] in queued for r in candidates):status='SPE '+next(queued[r['complex_directory']] for r in candidates if r['complex_directory'] in queued).lower()
            elif candidates:status='; '.join(sorted({r['status'] for r in candidates}))
            else:status='no matching Perlmutter calculation; check Kestrel sources'
            rr=relax.get(key,[])
            report.append(dict(surface=s,molecule=m,functional=f,E_ads_ML=em,E_ads_relaxed=er,E_ads_SPE=es,
                relaxed_status='available' if er is not None else '; '.join(sorted({r['status'] for r in rr})) or 'no matching Perlmutter calculation',
                SPE_status=status,completion_job_directory=recovery.get(key,{}).get('directory',''),
                dft_contact_available=bool(dft_match),required_CONTCAR=dft_sources.get(key2,'') if not dft_match else ''))
    with (root/'dft_completion_coverage.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(report[0]));w.writeheader();w.writerows(report)
    with (root/'dft_missing_geometry_sources.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(geometry[0]));w.writeheader();w.writerows(r for r in geometry if not r['dft_contact_available'])
    summary=dict(snapshot=datetime.now(timezone.utc).isoformat(),systems=len(geometry),functional_entries=len(report),
                 MLIP_images=sum(r['mlip_image'] for r in geometry),DFT_images=sum(r['dft_image'] for r in geometry),
                 DFT_contacts=sum(r['dft_contact_available'] for r in geometry),relaxed_energies=sum(r['E_ads_relaxed'] is not None for r in report),
                 SPE_energies=sum(r['E_ads_SPE'] is not None for r in report),SPE_status=dict(Counter(r['SPE_status'] for r in report)),
                 main_SPE_array='58164924',main_SPE_queue=dict(states),recovery_campaigns=campaigns)
    (root/'dft_completion_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    page=re.sub(r'<!-- completion-coverage -->.*?<!-- /completion-coverage -->\n?', '',page,flags=re.S)
    jobs=', '.join(c['job_id']+' ('+c['state']+')' for c in campaigns) or 'none'
    note=f'''<!-- completion-coverage -->
<div class="note info"><b>Completion status ({summary['snapshot'][:10]}):</b>
{summary['relaxed_energies']} / 1660 relaxed-energy entries and {summary['SPE_energies']} / 1660 SPE entries are available.
All 415 systems have images on both sides; 15 DFT contact measurements await their source CONTCARs.
<br><b>Recovery work:</b> two clean-slab references and two electronic SPE retries target five missing entries. Job {html.escape(jobs)}.
The existing SPE array continues separately. Pending results are not displayed as completed data.
<br><a href="dft_completion_coverage.csv" download>Every energy gap and its status</a> &middot;
<a href="dft_missing_geometry_sources.csv" download>15 required Kestrel CONTCAR paths</a> &middot;
<a href="dft_completion_summary.json">Coverage and scheduler snapshot</a>.</div>
<!-- /completion-coverage -->
'''
    page=page.replace('<h2>Per-system structure',note+'<h2>Per-system structure')
    (root/'dft_comparison.html').write_text(page)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
