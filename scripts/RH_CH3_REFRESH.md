# Rh–CH3 and full Perlmutter refresh — 2026-09-15 UTC

The Rh–CH3 BEEF-vdW relaxations were already electronically and ionically
converged. Their adsorption values were hidden by the historical ±5 eV review
screen, not by missing outputs or failed relaxation.

| Surface | Complex TOTEN (eV) | Relaxed slab TOTEN (eV) | Relaxed CH3 TOTEN (eV) | Raw adsorption energy (eV) |
| --- | ---: | ---: | ---: | ---: |
| Rh100 | -149.58705116 | -104.08110193 | -17.11758718 | -28.38836205 |
| Rh110 | -97.85787493 | -61.02501679 | -17.11758718 | -19.71527096 |
| Rh111 | -258.65683547 | -192.72273486 | -17.11758718 | -48.81651343 |

Sources are `dft_jobs/CH3_Rh*/beef_vdw/OUTCAR`,
`vasp_slab/Rh*/beef_vdw/OUTCAR`, and `vasp_mol/CH3/beef_vdw/OUTCAR`, relative
to `/pscratch/sd/j/jcho5/VASP`. The complexes use `PAW_PBE Rh 04Feb2005`, while
the slabs use `PAW_PBE Rh_pv 25Jan2005`. The project policy does not reject
potential variants, but that policy does not establish physical reference
equivalence. No energy offset, scaling, or correction is applied.

The page now displays these raw values with amber cells and a dagger (†), meaning
**converged calculation, adsorption energy needs review**. This treatment applies
to other converged, composition-matched Perlmutter results outside the review
range, including SPE results with verified original-POSCAR provenance. It does not
admit incomplete calculations, unmatched compositions, or missing references.
Previously published values retain priority. A later screened SPE result can
replace an explicitly provisional review value.

The three Rh–CH3 BEEF-vdW SPE inputs have no OUTCAR yet and are queued as tasks
`58164924_759`, `58164924_763`, and `58164924_767` at the inventory snapshot.
Their table cells say `queued`; relaxed energies are not substituted for SPE.

## Search and coverage

`inventory_perlmutter_components.py` searches every INCAR under `vasp_mol`,
`vasp_slab`, and `dft_jobs`, without ignore-file exclusions or recursive symlink
following. The snapshot contains 4,554 prepared calculations: 557 molecule,
189 slab, and 3,808 adsorbed-system calculations. NSW classifies 2,650 as
relaxations and 1,904 as SPE. All prepared molecule and slab calculations found
in these trees have NSW>0; no separate reference SPE results were found there.
The existing adsorption-SPE convention therefore continues to subtract converged
relaxed molecule and slab references.

`dft_perlmutter_components.csv` records total energy, convergence status, NSW,
potential identities, composition, output path/time, and queue state. Unconverged
total energies are diagnostic records, not published adsorption energies.
The inventory and adsorption audits have their own snapshot timestamps because
the batch continues to run during extraction.

The four isolated recovery calculations and export job have completed. Their
clean slabs can be reused for another adsorbate on the same surface only after
checking the functional, metal composition/count, full cell, and unchanged staged
input hashes. This removes the earlier restriction to the original target list;
the source calculations and manifests are unchanged. In particular, the completed
Ir100 slab also supplies the reference needed by the C2H6_Ir100 r2SCAN SPE retry.

The refreshed page has 1,311 relaxed entries (186 marked for energy review) and
893 SPE entries (43 marked for energy review). The screened relaxed audit accepts
715 candidates; the SPE audit selects 248 candidates before preservation of
existing published values. The main SPE array continues to run.

## Figures and reproducibility

Newly displayed review values carry `relaxed_review` and `SPE_review` flags in
`dft_comparison_singlepoint.csv`. Paired figures exclude rows with either flag;
the functional-SPE comparison excludes SPE review values. This changes no
historical number and adds no new magnitude filter to unflagged historical data.
The paired analysis contains 673 complete triples; the common four-functional
SPE cohort contains 125 systems. Missing values are never replaced with zero.

```bash
python scripts/inventory_perlmutter_components.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/extract_perlmutter_energies.py --vasp-root /pscratch/sd/j/jcho5/VASP --output dft_perlmutter_energy_audit.csv
python scripts/extract_perlmutter_singlepoint.py --project-root /pscratch/sd/j/jcho5/VASP --output-root .
python scripts/update_comparison_energies.py
python scripts/plot_energy_decomposition.py
python scripts/build_completion_report.py --project-root /pscratch/sd/j/jcho5/VASP
python -m unittest discover -s scripts -p 'test_*.py'
```

The CSV flags and amber styling survive regeneration. Regression checks cover
review display versus screened exports, preservation of published values,
replacement of a provisional SPE by a screened result, queue labels remaining
non-numeric, reference compatibility, and exclusion from figure statistics.
