# Functional-specific adsorption geometry and permissive energy display

The user requested that potential mismatches and unverified references no longer hide numerical results. `data/adsorption_results/display_policy.json` activates this policy in the existing cluster selection pipeline. Matching references remain preferred when available. Potential, slab-cell and core-setting differences remain attached to the downloadable source record; they do not suppress the displayed number. Complex and reference energies must exist, their calculations must pass the stored completion check, and functional, molecule identity and atom counts must agree. No slab energy is scaled and no reference total is invented. Earlier published values remain available as historical fallbacks. Flagged values remain excluded from the validated energy plots.

## Structure comparison

`functional_geometry.py` reads actual final CONTCAR coordinates for each of the four functionals, preferring converged Perlmutter structures and the selected energy's complex where possible. An available last geometry from an unfinished run is explicitly labeled. A different geometry source from the energy source is also labeled. It never uses one functional's image as another functional's structure.

MLIP coordinates come from the published SevenNet-OMNI or GOAD five-model CIF chosen by the existing card's model description. Formula/common-name aliases are resolved. Five static panels show the original MLIP image and separate PBE, PBE+D3, r²SCAN and BEEF-vdW images. `render_functional_figures.py` uses the original ASE renderer (rotation -70x,20y,10z; unit-cell outline; white background) on unmodified archived coordinates. Images are packed into one PNG atlas per surface for loading efficiency; `functional_figures.json` maps each panel to its exact image. Bond distances, angles and torsions remain in the tables. The archived JSON retains the cell and Cartesian coordinates for download; it is a geometry record, not a ready-to-run VASP input. Original published DFT images remain accessible in an expandable section.

Distances use ASE minimum-image periodic boundaries. Each molecular atom's nearest surface atom is reported, with both atom labels; nearest overall and nearest heavy-atom contacts are summarized. A geometric contact candidate uses 1.25 times the sum of covalent radii. This does not establish chemical bond order. All distances remain available, including atoms outside that cutoff.

Cross-model comparisons match heavy-atom molecular graphs, not raw atom indices. Graph edges use a 0.35 Å minimum and the same covalent-radius cutoff. Symmetry-equivalent mappings are resolved by centered rigid-fit RMSD over at most 256 mappings. Different molecular connectivity is not forced into a correspondence. Hydrogens with a single covalent attachment are paired within the same mapped heavy atom by minimum bond-direction difference. These symmetry-equivalent choices are geometric conventions, not unique chemical identities. Molecular labels are element plus 1-based adsorbate order; metal labels use the source file's 1-based global atom index, so metal labels across calculations do not assert a shared site identity.

Bond-angle and torsion tables use corresponding molecular bond paths. Heavy atoms are preferred for torsion endpoints; a terminal attached H provides a torsion when no heavy endpoint exists. H-containing angles are included for heteroatoms and carbon centers with fewer than two heavy neighbors. Surface angles are metal–anchor–neighbor, with each structure's own nearest metal. Torsions and their differences wrap to [-180,180) degrees. Undefined or unavailable quantities remain absent. Symmetric equivalent atoms can give different equally defensible paths; the selected mapping is downloadable.

`functional_geometry.json` contains raw and display coordinates, atom maps, distances, angles, source paths and SHA-256 hashes. Each generation saves a content-addressed compressed coordinate snapshot under `data/adsorption_results/geometry_snapshots/`, preserving structures locally even if HPC files later change. `functional_geometry_contacts.csv` provides the full distance table; `functional_geometry_audit.json` describes coverage and selected sources.

## Fresh Kestrel acquisition

This session could read Perlmutter directly. The official external endpoint `kestrel.nlr.gov` was reachable but did not accept noninteractive authentication; internal Kestrel DNS was unavailable. Therefore current Kestrel energy rows come from retained exports, not a fresh remote scan. Unavailable Kestrel coordinate files are explicitly marked; images are not reconstructed from energies.

After authenticated access is available, run the following from a checkout of this repository on Kestrel with ASE/NumPy available:

```sh
python scripts/cluster_result_bundle.py export \
  --cluster kestrel \
  --project-root /scratch/jcho5/goad-global-optimization \
  --output /scratch/jcho5/kestrel-results.json.gz
```

Transfer that bundle to the local checkout and merge it:

```sh
python scripts/cluster_result_bundle.py import --input /path/to/kestrel-results.json.gz
python scripts/cluster_selection.py
python -c 'import sys; sys.path.insert(0,"scripts"); import cluster_selection; cluster_selection.apply()'
python scripts/functional_geometry.py
python scripts/render_functional_figures.py
python scripts/plot_energy_decomposition.py
python scripts/build_completion_report.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/report_result_preservation.py
```

The export reads `vasp_mol`, `vasp_slab`, `vasp_slab_kestrel`, `dft_jobs`, `poscar/best`, and manifest-listed recovery calculations. It does not submit jobs. Import validates component and coordinate hashes, merges immutable observations, and caches coordinates by content hash without overwriting earlier snapshots.
