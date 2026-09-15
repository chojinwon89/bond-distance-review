# Data and figure completion — 2026-09-15 UTC

The comparison page now has 415 pairs of structure images. For 67 systems without
a SevenNet-OMNI image, the added image and minimum contact distance come from the
exact published base GOAD PNG and its `bond_distances.csv` record. Captions identify
this variant, and `mlip_base_structure_sources.json` records each source and image
hash. These additions do not change any energy value or its geometry assignment.
The base records do not contain coordination sites or original coordinates.

The new relaxed-energy audit accepts 682 results, 17 more than the previous audit.
The page preserves previously published numbers where the current audit cannot
replace them. The energy decomposition figures are rebuilt from the refreshed
table. The SPE audit selects 87 Perlmutter candidates; the two additional candidates
overlap retained published values, so the page still displays 759 SPE values.
`dft_completion_summary.json` and `dft_completion_coverage.csv` give the complete
coverage and queue snapshot. Scheduler states are snapshots, not live indicators.

## Recovery calculations

Campaign: `/pscratch/sd/j/jcho5/VASP/completion_runs/site_gapfill_20260914`.
Four isolated calculations were submitted in regular CPU job **58336200**:

| Calculation | Missing adsorption entries targeted |
| --- | --- |
| Cu100 clean slab, BEEF-vdW, 64 Cu atoms | C2H4_Cu100 and C2H6_Cu100 SPE |
| Ir100 clean slab, r2SCAN, 36 Ir atoms | C2H4_Ir100 SPE |
| C2H6_Ir100, r2SCAN, electronic SPE retry | C2H6_Ir100 SPE |
| C2H6_Rh100, BEEF-vdW, electronic SPE retry | C2H6_Rh100 SPE |

`prepare_completion_jobs.py` stages these cases only. Slabs retain the source
complex cell and metal-atom constraints, with adsorbates removed and the matching
metal POTCAR block. Slabs are relaxed; retries keep the exact original ML POSCAR
and NSW=0. NELM increases to 300 while the original electronic tolerance and
functional settings are retained. Input hashes, settings, target systems, and
submission details are recorded in the campaign manifests. Original calculations
are not modified. VASP pseudopotentials remain on Perlmutter.

The packed job uses one CPU node for at most six hours, with at most three
concurrent 32-rank calculations, following the
[NERSC guidance for multiple parallel jobs](https://docs.nersc.gov/jobs/examples/#multiple-parallel-jobs-while-sharing-nodes).
Each worker checks input hashes, refuses to overwrite existing outputs, and
records electronic convergence and, for slabs, ionic convergence. A queued or
failed calculation cannot supply a new published energy. The original SPE array
**58164924** continues independently.

Export job **58336271**, dependent on the recovery job ending, rebuilds an isolated
site snapshot at `campaign/reports/site` using the persistent Python environment
`/pscratch/sd/j/jcho5/envs/bond-review`. It runs even if a recovery calculation fails,
so failures are recorded. It does **not** publish to GitHub automatically. Inspect
the export log, `completion_result.json` files, energy audit, and rebuilt figures
before importing and publishing the resulting snapshot. The completion extractor
requires the normal convergence, composition, cell, geometry provenance, and
energy-screening checks; existing published SPE values retain priority over retries.

## Remaining source gaps

Fifteen DFT images have no corresponding contact measurement in the local archive.
`dft_missing_geometry_sources.csv` lists the exact Kestrel CONTCAR paths used to
render those images. The Kestrel tree is not mounted here. A matching geometry is
required to measure these contacts; a different relaxation or a distance estimated
from a PNG would not establish the missing value.

Other empty energy cells retain their specific reason in the coverage CSV:
queued SPE calculations, missing references, review-threshold failures, or no
matching Perlmutter calculation. This update does not represent those gaps as
completed results.

## Refresh commands

From this repository, using Python with ASE, NumPy, Matplotlib, and BeautifulSoup:

```bash
python scripts/extract_perlmutter_energies.py --vasp-root /pscratch/sd/j/jcho5/VASP --output dft_perlmutter_energy_audit.csv
python scripts/extract_perlmutter_singlepoint.py --project-root /pscratch/sd/j/jcho5/VASP --output-root .
python scripts/update_comparison_energies.py
python scripts/fill_missing_gallery_images.py
python scripts/plot_energy_decomposition.py
python scripts/build_completion_report.py --project-root /pscratch/sd/j/jcho5/VASP
python -m unittest discover -s scripts -p 'test_*.py'
```

`completion_support.py` discovers only jobs listed in an isolated campaign
manifest. It rejects directories outside that campaign. New successful slab
references are restricted to their declared target systems, functional, cell,
and atom count. SPE retries must retain the original ML POSCAR source and staged
input hashes; originals have priority when both are usable.
