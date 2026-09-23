"""Publish timestamped missing-component job status; never infer DFT convergence from Slurm."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json,shutil,subprocess
from pathlib import Path
from extract_kestrel_singlepoint import component

ROOT=Path(__file__).resolve().parents[1]

def report(campaign,root=ROOT):
    manifest=json.loads((campaign/'manifest.json').read_text());submission=json.loads((campaign/'submission.json').read_text());job_id=submission['job_id']
    accounting=subprocess.check_output(['sacct','--array','-X','-j',job_id,'-n','-P','-o','JobID,State,ExitCode,Elapsed'],text=True)
    states={}
    for line in accounting.splitlines():
        values=line.strip().split('|')
        if len(values)>=4:states[values[0]]=dict(zip(['job_id','state','exit_code','elapsed'],values))
    queue=subprocess.run(['squeue','--array','-j',job_id,'-h','-o','%i|%T'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    for line in queue.stdout.splitlines():
        jid,state=line.strip().split('|');states.setdefault(jid,dict(job_id=jid))['state']=state
    jobs=[]
    for j in manifest['jobs']:
        jid=job_id+'_'+str(j['array_index']);scheduler=states.get(jid,dict(job_id=jid,state='UNKNOWN'))
        result=component(Path(j['directory']));converged=result['status']=='converged'
        status='converged' if converged else 'running' if scheduler['state']=='RUNNING' else 'queued' if scheduler['state']=='PENDING' else 'awaiting status' if scheduler['state']=='UNKNOWN' else 'needs review'
        jobs.append(dict(job_id=jid,role=j['role'],surface=j['surface'],molecule=j['molecule'],functional=j['functional'],directory=j['directory'],target_cells=j['target_cells'],status=status,calculation_status=result['status'],scheduler=scheduler))
    folder=root/'data/adsorption_results'/campaign.name;folder.mkdir(parents=True,exist_ok=True)
    for name in ['manifest.json','submission.json','gap_plan.json','preflight_validation.json','preflight_components.jsonl.gz','joblist.txt','run_array.slurm']:
        dest=folder/name
        if dest.exists() and dest.read_bytes()!=(campaign/name).read_bytes():raise ValueError('Archived file changed: '+name)
        if not dest.exists():shutil.copyfile(campaign/name,dest)
    status=dict(schema_version=1,checked_at=datetime.now(timezone.utc).isoformat(),array_job_id=job_id,submitted_at=submission['submitted_at'],
                summary=dict(jobs=len(jobs),roles=dict(Counter(j['role'] for j in jobs)),statuses=dict(Counter(j['status'] for j in jobs))),
                archive=str(folder.relative_to(root)),jobs=jobs,
                note='A saved scheduler snapshot. New DFT structures starting from published MLIP coordinates are not assumed identical to unavailable Kestrel inputs. Previously published energies remain until verified results are imported.')
    for p in [root/'missing_component_status.json',folder/'status.json',campaign/'status.json']:p.write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status['summary'],indent=2))
    return status

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);a=p.parse_args();report(a.campaign)
