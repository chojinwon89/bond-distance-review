# Total energies, binding-energy gaps, and cluster imports

The component store separates a completed calculation's VASP free energy TOTEN
(eV) from a derived adsorption energy. A missing binding energy never deletes
the available complex total energy. Unfinished totals are retained with their
status for diagnosis; they are not promoted to converged references.

`dft_component_store.jsonl.gz` is the portable source of component observations.
`dft_component_energies.csv` projects the most recent output metadata for each
cluster/path into a spreadsheet. Both relaxed and NSW=0 calculations are included.
The full store retains superseded snapshots. Archived Kestrel records are labeled
as archived, not claimed as newly read files.

Each schema-v1 observation includes:

- Cluster and absolute calculation path, role, surface, original name and
  canonical molecule. Explicit aliases unify DME/CH3OCH3 without confusing DME
  and ethanol or other isomers with the same gross formula.
- Functional, INCAR settings, NSW, final available TOTEN, output convergence,
  composition, potential identities, and POSCAR cell when readable.
- Available output size/mtime and, for new reads, a head/tail SHA256 fingerprint.
  This is explicitly **not** a whole-OUTCAR hash. Existing archived exports keep
  their available metadata without invented hashes or extraction dates.
- A SHA256 snapshot ID over normalized data, with separate import provenance.
  Reimporting an identical observation is idempotent. Different clusters, source
  directories, and changed output snapshots are retained independently.

The CSV projection ranks snapshots by output mtime, then observation time.
Deleted/unmounted outputs do not erase previously captured results. Consult the
provenance and full snapshot history before treating an archived output as current.
No importer changes a published adsorption energy or overwrites the other cluster.

## Collect on Kestrel

Copy/clone this repository to Kestrel, activate Python 3.8+ with ASE, NumPy and
BeautifulSoup4, and ensure `rg` is installed. Use the repository's shared scripts;
an external `mol_canon.py` is no longer required.

```bash
python scripts/component_store.py --store kestrel-components.jsonl.gz collect \
  --cluster kestrel \
  --project-root /kfs3/scratch/jcho5/goad-global-optimization
```

The read-only collector discovers every INCAR beneath existing `vasp_mol`,
`vasp_slab`, `vasp_slab_kestrel`, `dft_jobs`, and `poscar/best` roots. Optional roots
may be absent. It includes nested `fully_relaxed` and `singlepoint` layouts;
NSW/settings, rather than directory spelling, establish the calculation type.
It does not recursively traverse symlink directories or submit calculations.
An empty/unavailable collection fails without replacing the store.

Use `--cluster perlmutter --project-root /pscratch/sd/j/jcho5/VASP` for Perlmutter.
The parser reads INCAR, POSCAR and OUTCAR/OUTCAR.gz. Settings and composition
checks are the same as those used by the existing SPE audit.

## Merge on the website host

Transfer the exported gzip file using your established cluster-transfer method,
then run these commands in the website repository:

```bash
python scripts/component_store.py merge --input /path/to/kestrel-components.jsonl.gz
python scripts/component_store.py report
```

The merge preserves old observations, validates snapshot checksums, and writes
atomically. The report regenerates the component CSV, gap explanation page,
summary JSON and `dft_kestrel_search_requests.json` from the current website
coverage snapshot. It does not promote candidate totals to binding energies.
Existing source-export JSON can also be imported:

```bash
python scripts/component_store.py import-sources --cluster kestrel \
  --source dft_kestrel_singlepoint_sources.json
```

## Search requests and reference matching

Each missing page cell has a request with surface, functional, relaxed/SPE mode,
canonical molecule, all recognized spellings, and available complex snapshot IDs.
For completed complexes, requests also contain required slab and molecule atom
counts, the complex cell, and IDs of converged reference candidates from either
cluster. Slab candidates must match the full cell within 1e-5 Å; molecule
candidates must match chemical identity and composition. Potential variants are
recorded and do not exclude candidates under the project's current policy.

These are **search candidates**, not automatic cross-cluster substitutions.
Check geometry provenance, atom constraints, k points, cutoff, spin, smearing,
potential choices and the intended reference convention before calculating a
new binding energy. The store does not currently capture every one of those
compatibility inputs. No energy scaling, offset or atom-count normalization is
performed. A request list is not a job submission list: search existing completed
outputs before considering new calculations.

The existing `extract_kestrel_singlepoint.py` remains the SPE derivation/export
path for its documented `poscar/best` layouts; the new collector covers a broader
set of nested relaxed and SPE directories for discovery and reuse. Once a fresh
Kestrel source export is available, its audited derived records can feed the
existing page updater. Preserve existing page values when new references fail.

## September 15, 2026 investigation

Ag111–DME PBE complex totals are −138.67977025 eV (relaxed) and
−138.64591398 eV (SPE), both converged. The standard Perlmutter clean slab has
64 Ag atoms, versus 36 in the complex. Both DME/CH3OCH3 gas reference outputs
ended with ZBRENT failure. A saved Kestrel Ag111_n36 PBE reference has the same
cell as the complex and TOTEN −91.56082669 eV; its saved PBE DME gas references
are also unfinished. Therefore a completed gas reference still needs locating.
The public explanation gives the component values and distinguishes all blank
categories. Run `python -m unittest discover -s scripts -p 'test_*.py'` to check
alias handling, merge history, reference candidates and existing energy logic.
