#!/usr/bin/env python3
"""Archive restart provenance and publish a timestamped Slurm/result snapshot."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from extract_kestrel_singlepoint import component

ROOT = Path(__file__).resolve().parents[1]


def report(campaign, root=ROOT):
    manifest = json.loads((campaign/'manifest.json').read_text())
    submission = json.loads((campaign/'submission.json').read_text())
    audit_bytes = (campaign/'source_audit.json').read_bytes()
    if hashlib.sha256(audit_bytes).hexdigest() != manifest['audit_sha256']:
        raise ValueError('Original audit checksum mismatch')
    audit = json.loads(audit_bytes)
    observed = datetime.now(timezone.utc).isoformat()
    job_id = submission['job_id']
    accounting = subprocess.check_output(['sacct','--array','-X','-j',job_id,'-n','-P',
        '-o','JobID,State,ExitCode,Elapsed,Timelimit,Start,End'], text=True)
    states = {}
    fields = ['job_id','state','exit_code','elapsed','time_limit','start','end']
    for line in accounting.splitlines():
        values = line.strip().split('|')
        if len(values)>=len(fields):states[values[0]]=dict(zip(fields,values))
    # Queue expansion is authoritative for pending/running tasks, including
    # array elements that have not yet acquired individual accounting rows.
    queue = subprocess.check_output(['squeue','--array','-j',job_id,'-h','-o','%i|%T'],text=True)
    for line in queue.splitlines():
        jid,state=line.strip().split('|')
        states.setdefault(jid,dict(job_id=jid))['state']=state
    jobs = {}
    for job in manifest['jobs']:
        jid = job_id+'_'+str(job['array_index'])
        scheduler = states.get(jid,dict(job_id=jid,state='UNKNOWN'))
        calculation = component(Path(job['directory']))
        # Slurm COMPLETED alone never establishes scientific convergence.
        converged = calculation['status']=='converged' and (calculation.get('nsw') or 0)>0
        display = ('converged' if converged else
                   'running' if scheduler['state']=='RUNNING' else
                   'queued' if scheduler['state']=='PENDING' else
                   'awaiting status' if scheduler['state']=='UNKNOWN' else
                   'needs review')
        jobs[job['source_directory']+'/CONTCAR'] = dict(
            source_sha256=job['source_sha256']['CONTCAR'],surface=job['surface'],
            molecule=job['molecule'],functional=job['functional'],
            previous_cause=job['cause'],restart_directory=job['directory'],
            job_id=jid,restart_status=display,calculation_status=calculation['status'],
            scheduler=scheduler,changes=job['changes'])
    folder = root/'data/adsorption_results'/campaign.name
    folder.mkdir(parents=True,exist_ok=True)
    # Keep immutable audit/inputs/submission copies; only the status snapshot moves.
    for name in ['manifest.json','submission.json','run_array.slurm','joblist.txt','cross_cluster_review.json']:
        source = campaign/name;dest = folder/name
        if dest.exists() and dest.read_bytes()!=source.read_bytes():
            raise ValueError('Refusing to replace archived '+name)
        if not dest.exists():shutil.copyfile(source,dest)
    compressed = gzip.compress(audit_bytes,mtime=0)
    archived = folder/'source_audit.json.gz'
    if archived.exists() and archived.read_bytes()!=compressed:
        raise ValueError('Refusing to replace archived audit')
    if not archived.exists():archived.write_bytes(compressed)
    result = dict(schema_version=1,checked_at=observed,campaign=campaign.name,
        array_job_id=job_id,submitted_at=submission['submitted_at'],
        summary=dict(total=len(jobs),causes=dict(Counter(r['cause'] for r in audit)),
                     statuses=dict(Counter(j['restart_status'] for j in jobs.values())),
                     scf_warning_jobs=sum(r['scf_limit_warning'] for r in audit)),
        kestrel=json.loads((campaign/'cross_cluster_review.json').read_text()),
        audit_directory=str(folder.relative_to(root)),jobs=jobs)
    (root/'relaxation_restart_status.json').write_text(json.dumps(result,indent=2)+'\n')
    (folder/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    (campaign/'status.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checked_at=observed,summary=result['summary']),indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',type=Path,required=True)
    args=parser.parse_args()
    report(args.campaign.resolve())
