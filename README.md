# Bond distance review

Adsorption structures, MLIP predictions, and DFT validation across metal surfaces.

- [Structure gallery](https://chojinwon89.github.io/bond-distance-review/)
- [DFT energies and geometry comparison](https://chojinwon89.github.io/bond-distance-review/dft_comparison.html)
- [Research update — 2 October 2026](https://chojinwon89.github.io/bond-distance-review/research_updates/2026-10-02/): four-functional single-point comparisons, selected independent pymatgen structures, the paired Ag(111)–glycerol initialization experiment, and Perlmutter/Kestrel campaign status.

The dated update includes downloadable numerical data, final structures, protocols,
source provenance, and an artifact inventory with SHA-256 hashes. Calculation status
is a timestamped snapshot, not a live cluster monitor. Existing comparison tables
are retained separately from these new campaign results.

## Interpreting the comparisons

DFT is the reference for assessing prediction accuracy. Keep DFT single-point
energies at MLIP geometries separate from DFT-relaxed energies, and distinguish
SevenNet-Omni results from newer SevenNet-MF experiments. Compare adsorption
energies using matched slab and molecule references, cells, potentials, functional,
and energy conventions. Existing review flags remain visible; an unflagged point
is not a new verification of its calculation.

## Maintaining the research snapshot

The dated snapshot is assembled with an explicit artifact allowlist:

```sh
python scripts/publish_research_snapshot_20261002.py /path/to/charge-guided-adsorption-search /path/to/current_beef_queue_snapshot.json
```

First refresh the source experiment's report and scheduler status. The publisher
copies figures, numerical records, final structures, and metadata; it does not
submit jobs or change the existing energy tables. For subsequent campaigns,
create a new dated snapshot and preserve the earlier records. Licensed POTCAR
files and private transfer archives must not be committed.

GitHub Actions publishes the static site when changes reach `main`.
