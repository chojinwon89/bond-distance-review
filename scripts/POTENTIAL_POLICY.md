> **Publication preservation update (2026-09-22):** Previously published numbers now remain visible with historical audit labels when no validated replacement exists. They stay excluded from validated plots. This supersedes earlier display-withholding wording below; scientific compatibility checks are unchanged. See [local result archive](../data/adsorption_results/README.md).

# Matching-potential reference audit

The latest user instruction requires the same potential for each species in the
complex and its references. This supersedes the earlier instruction permitting
POTCAR variants. A physically plausible energy or a converged calculation does
not validate a mixed-potential subtraction.

The component matcher now requires equal recorded PAW TITEL identities, including
element, variant and dataset date. Missing identities fail closed. This checks
archived dataset identity; it does not claim byte-identical POTCAR files where
the original files are unavailable. Newly staged slabs copy the exact metal
POTCAR block from a local complex and retain SHA256 hashes of their inputs.
Gas reference potentials must also match the complex's adsorbate potentials.

The frozen `dft_potential_candidates_before.csv` retains 1,719 candidate energies
from before this policy. Of these, 1,120 use mismatched metal potentials, spanning
999 distinct relaxed/SPE cells across 14 surfaces of Ag, Cu, Pd, Pt and Rh.
These are candidate counts, not counts of independent calculations. The mismatch
audit records the original energies, component IDs, potential identities and
compatible alternative reference IDs. It includes Perlmutter and archived Kestrel
records. Perlmutter was rescanned; no new authenticated Kestrel read is claimed.

Every numeric energy in the active comparison must now come from a selected,
validated component set. `potential mismatch` labels known mismatched references
without a validated replacement. `refs unverified` labels historical numbers for
which no validated component set is selected; it does not assert a proven
potential mismatch. Before/after CSVs and immutable component snapshots retain
the original numbers. Figures use the resulting validated table and retain the
existing exclusion of magnitude-review values.

## Recovery

51 missing clean-slab calculations were staged and submitted as array 58411454,
indices 0–50, with at most four running concurrently. Each requests one node,
32 MPI ranks, 120 GB and a two-hour wall limit in the shared queue. Inputs live at
`/pscratch/sd/j/jcho5/VASP/completion_runs/matched_potential_slabs_20260915`.
`dft_potential_slab_recovery.json` records every source, input hash, setting change,
target and submission command. POTCAR contents are not published.

Slabs are metal-only subsets of source complex POSCARs, preserving the cell and
constraints. The exact metal POTCAR block and k-point file are copied. Jobs are
deduplicated by surface, atom count, cell, constraints, potential hash, functional
and electronic settings. IALGO=38 uses the Davidson electronic solver; the
original electronic tolerance and functional parameters remain. Existing
BEEF-vdW parameter choices are not silently changed: correcting those requires a
separate consistent set of complex, slab and molecule calculations.

The worker refuses to overwrite outputs and verifies staged hashes. Scheduler
completion alone is insufficient; electronic and ionic convergence are checked.
Collectors discover these jobs through their campaign manifests, and matching
slabs can be reused across adsorbates after validation.

To refresh completed jobs and regenerate:

```bash
python scripts/component_store.py collect --cluster perlmutter --project-root /pscratch/sd/j/jcho5/VASP
python scripts/extract_perlmutter_singlepoint.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/audit_potential_references.py
python scripts/cluster_selection.py
python scripts/update_comparison_energies.py
python scripts/plot_energy_decomposition.py
python scripts/build_completion_report.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/component_store.py report
```

Publishing remains a separate step after reviewing and testing the generated data.
