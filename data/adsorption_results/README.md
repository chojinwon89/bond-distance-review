# Local adsorption result archive

The 2026-09-22 recovery separates saved observations from the current validation decision. A failed or incomplete future scan must not delete a published observation or turn a previously numeric cell into a blank.

- `publications.jsonl`: append-only, SHA-256 checked observations of published relaxed and SPE adsorption energies. Includes canonical molecule name, surface, functional, mode, eV value, original display precision, source Git commit, original cell/tooltip, ML display value, and original SPE source metadata where available. Missing historical provenance stays unknown; a source commit is not proof of a cluster or reference match.
- [Component snapshots](../../dft_component_store.jsonl.gz): immutable calculation records from Perlmutter and archived Kestrel, including component energies, available settings, potential identities, paths and structure metadata. This is the existing local component archive; it was not erased by the display refresh.
- `manifest.json`: paths, hashes and recovery counts. Preserve this entire repository when transferring/backing up the archive, including the component store at the repository root.
- [Current selection](../../dft_cluster_selection.csv): replaceable projection, including which values are currently validated or retained as historical. This is not the immutable history.

## Recovery and future updates

The bootstrap imports both public commits `23048650b883c5118871d1829a317187a65f0989` (before the broad withholding) and `9d88d3812666d5d0822f160999b9bdf4f0925a02` (fresh validated results). No raw HPC files were deleted. Old and corrected numbers coexist; corrected validated results win on the page.

`cluster_selection.apply()` now checkpoints numeric observations before selection, retains historical values when no validated replacement exists, and checkpoints validated selections afterward. These steps are part of the existing website refresh pipeline. Historical cells carry visible audit labels and record IDs linking them to this ledger. They are excluded from validated figure statistics. A known potential mismatch is not a physically valid adsorption energy merely because its old number remains visible.

Future Perlmutter and Kestrel collection should merge into the existing `component_store.py` archive, never replace it with a single-cluster export. Run the existing selection/render pipeline against the merged archive. Perlmutter remains preferred; Kestrel fallback still requires compatible references and, where applicable, verified structure equivalence. Failed scans remain observations and cannot remove published history. This recovery reused saved Kestrel records; it does not claim a fresh Kestrel login.

To import another historical publication explicitly:

```sh
python scripts/publication_store.py --commit <full-git-commit>
```

The import is idempotent and refuses a ledger with a checksum error. Back up the repository with Git and a separate copy; a directory on the same HPC alone is not an independent backup. Full OUTCAR/CONTCAR files remain on their source clusters; this archive stores extracted results and available provenance, not every raw VASP file.

After rendering a future refresh, run `python scripts/report_result_preservation.py` to regenerate the manifest and current counts. The manifest records the archive state when it was generated; individual ledger records are always independently checksum-verified by the reader.

## Current permissive display and geometry archive

The user's subsequent instruction enables `display_policy.json`: potential/reference verification flags no longer suppress numeric results. Matching available references remain preferred. Component totals, molecule identity, functional and atom counts are still required to derive a new value. See [functional geometry and acquisition workflow](../../scripts/FUNCTIONAL_GEOMETRY.md).

`geometry_snapshots/` saves immutable, compressed, content-addressed copies of the actual structure coordinates, atom mappings and measured geometry. The current projection is `functional_geometry.json` at the repository root. Future authenticated Kestrel bundles merge both component totals and coordinates through `cluster_result_bundle.py`; archived Kestrel metadata alone cannot supply missing atomic positions.
