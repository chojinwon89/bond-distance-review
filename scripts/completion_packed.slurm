#!/bin/bash
#SBATCH -A m5281
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH -J dft_page_gapfill
set -uo pipefail
CAMPAIGN="${1:?Pass absolute campaign directory}"
WORKER="$CAMPAIGN/completion_array.slurm"
pids=()
# At most three concurrent MPI applications, respecting the Slingshot VNI limit.
for index in 0 1 2; do
    SLURM_ARRAY_TASK_ID="$index" bash "$WORKER" "$CAMPAIGN/joblist.txt" > "$CAMPAIGN/logs/task_${SLURM_JOB_ID}_${index}.out" 2>&1 &
    pids+=("$!")
done
rc=0
wait "${pids[2]}" || rc=1
SLURM_ARRAY_TASK_ID=3 bash "$WORKER" "$CAMPAIGN/joblist.txt" > "$CAMPAIGN/logs/task_${SLURM_JOB_ID}_3.out" 2>&1 &
pids[2]=$!
for pid in "${pids[@]}"; do wait "$pid" || rc=1; done
exit "$rc"
