> **Publication preservation update (2026-09-22):** Previously published numbers now remain visible with historical audit labels when no validated replacement exists. They stay excluded from validated plots. This supersedes earlier display-withholding wording below; scientific compatibility checks are unchanged. See [local result archive](../data/adsorption_results/README.md).

# Perlmutter-first source selection

The user's current matching-potential policy supersedes the initial fill-only and variant-permissive update policies:
prefer Perlmutter; consider Kestrel when a result is missing or unusual, using
the same structure. Historical totals and previous displayed values remain
available in the component store and before/after CSVs.

The implemented selection order is:

1. A completed, referenced Perlmutter candidate inside the existing ±5 eV review
   band takes priority. This is a screening convention, not a guarantee of
   physical accuracy or numerical convergence.
2. For a Perlmutter candidate outside that band, an unflagged Kestrel candidate
   can replace it only when a structure fingerprint verifies equivalence.
   Selecting an energy because it is lower or closer to ML is not allowed.
3. Without a verified better alternative, retain the unusual Perlmutter result
   with its review marker, provided all references use matching potentials.
4. Kestrel can fill a missing result when Perlmutter has no usable candidate.
   If a corresponding Perlmutter structure exists, the geometry match is required.
   Historical numbers without a validated matching-potential component set are withheld; their previous values remain in the audit downloads.

The threshold is explicit and configurable in the existing energy auditors.
Positive or unusual total energies are never discarded based on sign. POTCAR
identities must match for all reference species, including variant and dataset date.
No energy offsets, scaling or atom-count corrections are applied.

## Same-structure evidence

The collector now records raw POSCAR/CONTCAR SHA256 hashes and normalized geometry
fingerprints. A normalized fingerprint contains ordered chemical symbols,
Cartesian positions rounded to 1e-4 Å, cell entries rounded to 1e-6 Å, and
selective-dynamics constraints. Identical raw files or identical fingerprints
provide conservative matching evidence. Atom reorderings or coordinate shifts
can produce false nonmatches; they are held for review rather than assumed equal.
For SPE, input geometry must match. For completed relaxations, final geometry
must match. A shared formula, surface, cell, directory name or similar total energy
does not establish coordinate equivalence.

Older Kestrel exports did not include these fingerprints. Their total energies
remain usable as archived component records under the existing reference rules,
but they cannot certify a cross-cluster complex-geometry replacement. The new
collector adds the required geometry evidence on the next Kestrel export.

## Reproduce

Collect on Kestrel using the user's active root:

```bash
python scripts/component_store.py --store kestrel-components.jsonl.gz collect \
  --cluster kestrel --project-root /scratch/jcho5/goad-global-optimization
```

Transfer that export to the website workspace, then:

```bash
python scripts/component_store.py merge --input /path/to/kestrel-components.jsonl.gz
python scripts/cluster_selection.py
python scripts/update_comparison_energies.py
python scripts/plot_energy_decomposition.py
python scripts/build_completion_report.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/component_store.py report
```

Run `extract_perlmutter_singlepoint.py` and import its source JSON after new local
SPE jobs finish. The collector itself can also scan Perlmutter. Both extractors now
capture structure identities. `update_comparison_energies.py` reapplies the cluster
policy after the legacy and shared-reference updates.

Downloads:

- `dft_cluster_candidates.csv`: every matched component-derived candidate.
- `dft_cluster_duplicates.csv`: Perlmutter/Kestrel energies and geometry evidence.
- `dft_cluster_selection.csv`: candidate choice, reason and actual application
  status for all 3,320 relaxed/SPE table cells. A candidate and an applied value
  are distinct when an old number is withheld pending compatible references.
- `dft_cluster_selection_before.csv`: table values before the first policy pass.
- `dft_cluster_selection_changes.csv`: numeric changes with their reasons.
- `dft_missing_reference_priorities.csv`: missing references grouped by affected
  page cells, allowing one gas or clean-slab calculation to resolve several gaps.

## Ag111–DME findings and recovery

The saved PBE and PBE+D3 adsorbed complexes are converged on both clusters. Their
saved gas references failed with ZBRENT. A size-matched clean slab alone does not
complete those binding energies. The BEEF-vdW complex SPE totals differ by roughly
1e-6 eV across clusters, so switching that complex would leave the unusual
binding energy essentially unchanged. No better, geometry-verified Kestrel
substitute is available in the archived records.

Two isolated Perlmutter gas-reference retries were submitted as array 58376018,
tasks 0 and 1, under
`/pscratch/sd/j/jcho5/VASP/completion_runs/dme_references_davidson_20260915`.
They use the converged BEEF-vdW gas geometry as a starting geometry, then optimize
with the original PBE/PBE+D3 functional, cutoff, potentials and k points. They do
not substitute a BEEF energy for either functional. RMM-DIIS replaces the failed
conjugate-gradient line search, with POTIM=0.2, EDIFF=1e-6 and the original force
threshold; all changes and input hashes are recorded in the campaign manifest.
The initial retry array 58374583 hit NELM=200 on both first ionic steps with
electronic IALGO=48 and was stopped. Its outputs remain in the original recovery
campaign. The new array additionally uses electronic IALGO=38 (blocked Davidson),
the robust default described in the [VASP ALGO documentation](https://vasp.at/wiki/index.php/ALGO).
The choice of ionic RMM-DIIS near a converged geometry follows the
[VASP optimizer documentation](https://vasp.at/wiki/index.php/IBRION).

Each task requests 16 MPI ranks, 60 GB and 30 minutes in the shared queue. Original
outputs are preserved. A scheduler completion alone is insufficient: the output
must pass electronic and ionic convergence checks before reuse.

The Kestrel external SSH endpoint was reached during this investigation, but
authentication was unavailable from the Perlmutter session. No fresh Kestrel
filesystem inspection is claimed by this publication.
