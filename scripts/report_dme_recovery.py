"""Export the status of the two isolated DME gas-reference retries."""
import json
from datetime import datetime,timezone
from pathlib import Path
import subprocess
from extract_kestrel_singlepoint import component

ROOT=Path(__file__).resolve().parents[1]
CAMPAIGN=Path('/pscratch/sd/j/jcho5/VASP/completion_runs/dme_references_davidson_20260915')


def main():
    manifest=json.loads((CAMPAIGN/'manifest.json').read_text())
    submission=json.loads((CAMPAIGN/'submission.json').read_text());job_id=submission['job_id'];states={}
    result=subprocess.run(['squeue','-r','-j',job_id,'-h','-o','%i|%T'],capture_output=True,text=True)
    for line in result.stdout.splitlines():
        job,state=line.split('|');states[job]=state
    accounting=subprocess.run(['sacct','-j',job_id,'--noheader','--parsable2','--format=JobID,State'],capture_output=True,text=True)
    for line in accounting.stdout.splitlines():
        fields=line.split('|')
        if len(fields)>=2 and '.' not in fields[0]:states.setdefault(fields[0],fields[1].strip())
    rows=[]
    for index,job in enumerate(manifest['jobs']):
        c=component(Path(job['directory']));task=job_id+'_'+str(index)
        rows.append(dict(functional=job['functional'],molecule='DME',job_id=task,
            scheduler_state=states.get(task,'not in queue'),calculation_status='not_started' if c['energy'] is None and states.get(task)=='PENDING' else c['status'],
            energy_TOTEN_eV=c['energy'],directory=job['directory'],input_changes=job['changes'],
            geometry_source=job['geometry_source'],target_systems=job['target_systems']))
    previous=CAMPAIGN.parent/'dme_references_20260915/stopped_for_retry.json'
    report=dict(snapshot=datetime.now(timezone.utc).isoformat(),jobs=rows,
        previous_attempt=json.loads(previous.read_text()) if previous.exists() else None,
        policy='Gas-only PBE/PBE+D3 retries. Original functional and potentials retained. A queued job is not an energy; completed outputs must pass electronic and ionic convergence checks.')
    (ROOT/'dft_gas_reference_recovery.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
