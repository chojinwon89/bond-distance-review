# Missing Kestrel DFT images

This update adds only the 15 previously missing DFT structure images to
`dft_comparison.html`. Existing PNGs, geometry metrics, energy tables, energy
CSVs, single-point data and energy-filtering code are preserved.

Sources are recursively discovered under
`/kfs3/scratch/jcho5/goad-global-optimization/poscar/best/`, including carbon-group
subdirectories and the `singlepoint/<functional>/fully_relaxed` layout.
Only exact surface/molecule directory names matching a page placeholder qualify.
Formate is CHO2; formic acid is CH2O2. They are never aliases.

## Image selection

- Ag111 formate uses the requested BEEF-vdW example.
- Other targets prefer PBE, then PBE+D3, r2SCAN, BEEF-vdW among verified runs.
- Equal-functional candidates prefer the shallowest source path, then lexical order.
- Selection never uses an adsorption energy or changes its filtering.
- All available candidates, rejection reasons, selected sources, source SHA-256
  hashes and the selected image hashes are recorded in `dft_structure_sources.json`.
- PBE+D3 inputs retain their actual IVDW and damping parameters in the audit.
  IVDW=11 is D3 zero damping; IVDW=12 is D3(BJ), per the
  [VASP documentation](https://vasp.at/wiki/index.php/IVDW).

## Verification and rendering

The updater checks nonempty CONTCAR, OUTCAR, INCAR and POSCAR; an exact INCAR
SYSTEM match; active ionic relaxation; input and executed functional/settings;
normal completion; final electronic convergence; the ionic convergence marker;
matching atom counts and molecular composition; finite positions; agreement
between CONTCAR and final OUTCAR cell/positions (periodic tolerance 0.00002 A);
and final forces after applying selective-dynamics constraints against EDIFFG.
Surface identity comes from the exact system path and SYSTEM label, with the
metal composition checked independently. The calculation's crystallographic
orientation is not independently re-indexed.

Rendering uses ASE's PNG writer with `rotation='-70x,20y,10z'` and
`show_unit_cell=2`, matching the existing `render_dft_structures.py` in the GOAD
repository. Coordinates and constraints come directly from CONTCAR. No AI image
generation or geometry modification is used. The functional appears in the image
caption and alt text; the caption links to the source audit.

This verifies completion at the calculation's configured tolerances, not
convergence with respect to numerical parameters. Images do not imply a new
adsorption-energy result or a refresh of the older contact-distance metrics.

## Reproduce

Dependencies: ASE, NumPy, Matplotlib, BeautifulSoup4.

```bash
python scripts/add_missing_dft_images.py \
  --source-root /kfs3/scratch/jcho5/goad-global-optimization/poscar/best \
  --audit-output /tmp/dft-image-audit.json
python scripts/add_missing_dft_images.py \
  --source-root /kfs3/scratch/jcho5/goad-global-optimization/poscar/best --apply
python -m unittest discover -s scripts -p 'test_*.py'
```

The default run only audits. `--apply` writes website files only, refuses to
overwrite an existing image, and asserts that all energy tables are unchanged
byte for byte. Once all placeholders are filled, another run makes no changes.
To reproduce the original addition, run against the parent website revision with
these scripts present. No VASP files are written and no jobs are submitted.
