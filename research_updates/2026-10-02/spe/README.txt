GOAD + SevenNet-Omni versus DFT single-point adsorption energies

Source: https://chojinwon89.github.io/bond-distance-review/dft_comparison.html
Snapshot: 2 October 2026. Download times, SHA256 hashes and HTTP modification
dates are recorded in sources/download_manifest.json.

The plotted x column is E_ads_DFT_SP; y is E_ads_ML_SP from the website's
dft_comparison_singlepoint.csv. These are adsorption energies, not raw VASP
total energies. The x value is NOT the DFT-relaxed adsorption energy or the
website's unusual ML-minus-SPE-minus-relaxed-DFT column. Original source
calculator metadata identifies SevenNet-Omni; this is not SevenNet-MF.

All 1450 finite pairs were checked against the displayed HTML values at
their published precision, and ML minus SPE was checked against delta_SP.
There are no duplicate surface/molecule/functional keys. The 74 missing
MLIP entries are exported separately and are not backfilled from other
datasets or geometries. Website rounding of MLIP values is preserved.

PBE: 370 pairs; PBE+D3: 351; r2SCAN: 363; BEEF-vdW: 366.
The full-range plot includes every pair. The zoom plot covers -5 to +1 eV
on both axes and leaves 1 PBE+D3 and 107 BEEF-vdW points outside the view.
Hollow markers mean SPE_review=true in the source CSV. No observations
were removed for audit flags. analysis_eligible is preserved in the CSV
but not used to filter this SPE comparison: it also depends on relaxed
comparison availability, which is not required for these axes.

The snapshot contains published values with unresolved reference/settings
issues, including BEEF-vdW extremes down to -42.800 eV. These plots reproduce
the website; they do not revalidate reference compatibility, geometry hashes,
VASP convergence, or historical BEEF functional settings. Differences cannot
all be attributed to intrinsic MLIP error. Archived energy/provenance columns
are retained in plotted_pairs.csv. No HPC jobs or website files were changed.

The website describes its SPE adsorption energies as complex single-point
energies on MLIP structures minus relaxed slab and gas-molecule references.
The reference convention therefore matters in interpreting MLIP-DFT differences.

PNG: raster figures; PDF/SVG: vector figures; interactive HTML: self-contained
Plotly display with every molecule/surface identified on hover and a full-range
button. Cu(111)/DME is annotated as a common example, not an exclusive subset.
DME means CH3OCH3. Statistics in validation.json describe all published pairs,
including flagged values, and are not screened model-accuracy estimates.

Reproduce with plot_spe_comparison.py using Python with numpy, matplotlib,
BeautifulSoup4, and plotly. Downloaded inputs are stored in sources/.
