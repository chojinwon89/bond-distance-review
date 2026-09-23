# Unfinished relaxation audit and restarts — 23 September 2026

All 169 selected Perlmutter geometries labelled unfinished were checked against the complete original OUTCAR files. None contains the normal-completion or ionic-convergence marker. The original output hashes and Slurm log evidence are preserved in `source_audit.json.gz`.

| Cause | Jobs |
|---|---:|
| Slurm time limit | 146 |
| ZBRENT geometry line-search failure | 22 |
| Interrupted without identified scheduler cause | 1 |

Eight jobs also contain electronic self-consistency warnings; these overlap the cause counts. The functional totals are 143 r²SCAN, 22 BEEF-vdW, 3 PBE and 1 PBE+D3.

## Kestrel comparison

No completed same-system relaxation was found in the merged component archive or the separate saved Kestrel relaxation audit. The component archive includes 2,225 Kestrel complex records, mostly single-point calculations; the separate audit contains 120 relaxation candidates, of which 55 were accepted. A completed single-point calculation does not establish geometry convergence.

Canonical molecule names, surface and functional were compared first. No completed same-system candidates remained for coordinate-fingerprint comparison. This is an archived-data check, not proof that Kestrel has no completed copies. Fresh SSH to `jcho5@kestrel.nlr.gov` failed authentication. `cross_cluster_review.json` records the archive hashes and this limitation.

## Restart campaign

Slurm array **58782557**, indices 0–168, maximum four simultaneous jobs; account m5281, CPU shared queue, 32 MPI ranks, two hardware threads per rank, 120 GB memory, 12-hour limit per job. No automatic resubmission.

Campaign directory: `/pscratch/sd/j/jcho5/VASP/completion_runs/unfinished_relaxations_20260923`.

- Each restart copies the last saved CONTCAR byte-for-byte to a new POSCAR. Original calculation directories are unchanged.
- Potentials, k points, atom order, cell, fixed-atom constraints, functional, cutoff, smearing, spin settings and force-convergence threshold are preserved.
- All jobs start electronic self-consistency from atomic charges (`ISTART=0`, `ICHARG=2`), without a reused WAVECAR. Wavefunction/charge-file writing is disabled.
- The 22 ZBRENT failures use `EDIFF=1e-6`, ionic `IBRION=1`, and `POTIM=0.2`.
- The eight jobs with electronic-convergence warnings use `ALGO=Normal` and `NELM=300`, removing the conflicting original `IALGO` entry.
- `manifest.json` records every input hash, source hash and setting change. The runner checks staged hashes, refuses existing outputs, and requires both electronic and ionic convergence before reporting completion.

The optimizer choice follows the [VASP ionic-relaxation documentation](https://vasp.at/wiki/IBRION); the electronic solver is described in the [VASP ALGO documentation](https://vasp.at/wiki/ALGO). These are targeted restart settings, not a guarantee of convergence.

## Status and future imports

`status.json` and the top-level `relaxation_restart_status.json` are timestamped snapshots, not live scheduler feeds. Slurm COMPLETED alone does not mean DFT convergence. Current figures remain the original saved geometry until a new verified geometry is imported.

Refresh the restart status with:

```bash
/pscratch/sd/j/jcho5/envs/bond-review/bin/python scripts/report_relaxation_restarts.py --campaign /pscratch/sd/j/jcho5/VASP/completion_runs/unfinished_relaxations_20260923
```

The component collector now discovers manifest-listed `relaxed` restart jobs. Subsequent exports and merges include their final energies and coordinates while retaining the original Perlmutter and Kestrel snapshots. The ordinary energy and geometry refresh must still run after import; this submission does not fabricate converged values or change the current tables.
