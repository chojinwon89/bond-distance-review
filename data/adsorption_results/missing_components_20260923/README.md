# Missing-component calculations on Perlmutter — 23 September 2026

Slurm array **58795983** contains **354 jobs**: 9 gas-molecule references, 2 clean-slab references, 34 fixed-geometry single-point calculations, and 309 geometry relaxations. It is separate from the existing 169-job restart array 58782557. Each array is limited to four simultaneous tasks; each new task requests 32 MPI ranks, two hardware threads per rank, 120 GB memory and a 12-hour walltime on m5281 / CPU shared. No automatic retries are enabled.

## Why these jobs

The saved page had 398 blank energy cells (262 relaxed, 136 SPE) and 309 absent per-functional geometries. Their union produces 593 calculation targets. Fresh Perlmutter outputs plus archived Kestrel records classified them as:

- 343 targets needing a new complex calculation: 309 relaxations and 34 SPE calculations.
- 70 targets already covered by the earlier restart campaign; no duplicate jobs submitted.
- 180 targets with completed complexes. Of these, 15 already have usable references in the fresh snapshot, 143 lack only a molecular reference, 20 lack only a slab reference, and 2 lack both. These are retained for the next verified energy refresh.

The nine gas calculations are CO2/PBE+D3, NO/r²SCAN, O2/all four functionals, acetic acid/PBE and PBE+D3, and glycerol/BEEF-vdW. The two slab calculations are Cu100/BEEF-vdW and Ir100/r²SCAN. Reference inputs are deduplicated by identity, composition, functional, potential identities, and (for slabs) cell. The gas references serve 145 currently blank cells; slab references serve 22, with two cells overlapping.

## Starting structures and constraints

54 complex jobs reuse native local input sets (25 relaxations and 29 SPE jobs). The other 289 complex jobs (284 relaxations and 5 SPE jobs) start from the published MLIP structures of 71 systems. These are explicitly new Perlmutter validations, not claimed reproductions of unavailable Kestrel input files.

For each new structure, the surface template has the same metal, metal atom count and cell. Its constraint mask is checked to fix only lower metal layers, and that layer convention is transferred to the MLIP structure. No molecular atoms are fixed. A layer boundary must have at least 0.5 Å separation. Potentials follow available same-system archived identities where present; otherwise they come from the chosen local surface template and recorded molecular donors. POTCAR block order is validated against POSCAR species order.

Native SPE jobs retain the original input POSCAR and disable ionic motion (`NSW=0`, `IBRION=-1`). Relaxations retain the reference force threshold (normally 0.05 eV/Å) and use `IBRION=1`, `POTIM=0.2`, `EDIFF<=1e-6`. All new jobs use the electronic Davidson solver (`ALGO=Normal`, removing `IALGO`) with `NELM>=300`. Functional definitions, including dispersion parameters, are inherited from same-functional templates. These choices do not guarantee convergence.

All source files remain unchanged. Preflight checked 512 source-file hashes, every staged input, atom counts, potential ordering, and nonoverlap with the existing restart campaign. The runner checks input hashes, refuses to overwrite outputs, and requires electronic convergence plus ionic convergence for every relaxation, including gas molecules.

## Archived records

- `manifest.json`: per-job source paths, input hashes, potential-donor hashes, constraints, settings changes, and affected table entries.
- `gap_plan.json`: every target, its existing source/status, and whether it needs a new calculation.
- `preflight_components.jsonl.gz`: immutable fresh Perlmutter component observations used for planning. Existing Kestrel observations remain in the repository's component store; live Kestrel authentication was unavailable.
- `preflight_validation.json`: input-integrity and duplicate checks.
- `submission.json`, `joblist.txt`, `run_array.slurm`: exact submission and task mapping.
- `status.json`: timestamped scheduler and VASP status, also exposed as the top-level `missing_component_status.json`.

Runtime directory: `/pscratch/sd/j/jcho5/VASP/completion_runs/missing_components_20260923`.

Refresh status with:

```bash
/pscratch/sd/j/jcho5/envs/bond-review/bin/python scripts/report_missing_components.py --campaign /pscratch/sd/j/jcho5/VASP/completion_runs/missing_components_20260923
```

The existing manifest-aware component collector and coordinate export can ingest these jobs as they finish. A scheduler completion alone is insufficient: electronic convergence and, for relaxations, ionic convergence must be verified. Existing energies, figures, and structural measurements are preserved until verified output import and regeneration.
