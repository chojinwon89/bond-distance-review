#!/bin/bash
#SBATCH -A m5281
#SBATCH -C cpu
#SBATCH -q shared
#SBATCH -n 32
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH -J dft_page_gapfill
set -euo pipefail
JOBLIST="${1:?Pass campaign joblist}"
JOBDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$JOBLIST")
[[ -n "$JOBDIR" ]] || exit 1
cd "$JOBDIR"
mkdir .completion_running
trap 'rmdir .completion_running' EXIT
[[ ! -e OUTCAR && ! -e vasprun.xml ]] || { echo 'Existing outputs; refusing overwrite'; exit 1; }
python3 - <<'PY'
import hashlib,json
from pathlib import Path
record=json.load(open('completion_inputs.json'))
for name, expected in record['staged_sha256'].items():
 assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==expected, name
PY
module load vasp-tpc/6.4.2-cpu
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
set +e
srun --exclusive --exact -N 1 -n 32 -c 2 --mem=120G --cpu-bind=cores vasp_std
rc=$?
set -e
python3 - "$rc" <<'PY'
import json,re,sys
from pathlib import Path
record=json.load(open('completion_inputs.json'))
p=Path('OUTCAR'); text=''
if p.exists():
 with p.open('rb') as f:
  f.seek(max(0,p.stat().st_size-1000000));text=f.read().decode(errors='replace')
e=list(re.finditer(r'free  energy   TOTEN\s*=\s*([-+\d.Ee]+)',text))
start=e[-2].end() if len(e)>1 else 0
ok=int(sys.argv[1])==0 and bool(e) and 'General timing and accounting' in text and 'aborting loop because EDIFF is reached' in text[start:e[-1].start()]
if record['role'] in ('slab','relaxed'):ok=ok and 'reached required accuracy' in text
result=dict(completed=ok,role=record['role'],vasp_exit_code=int(sys.argv[1]),energy_TOTEN_eV=float(e[-1].group(1)) if e else None)
Path('completion_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
sys.exit(0 if ok else 1)
PY
