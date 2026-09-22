> Current policy (22 September 2026): matching PAW identities are required. The counts and permissive-potential rules below describe an older snapshot and are superseded by [the current refresh](REFRESH_20260922.md) and [potential policy](POTENTIAL_POLICY.md).

# Perlmutter SPE refresh — 2026-09-14

The Perlmutter batch `58164924` evaluates the original ML POSCAR with NSW=0,
not the relaxed CONTCAR. Input provenance is recorded in each calculation's
`spe_inputs.json`. The extraction reads inputs/outputs without changing or
submitting calculations.

Snapshot: 1,904 prepared calculations, 226 OUTCARs present, 221 completed with
final electronic convergence. The audit selects 85 adsorption energies; 58
already have published SPE values, which are preserved, and 27 fill missing
entries. The page now displays 759 SPE values: 418 original published, 314 Kestrel,
and 27 Perlmutter. New values by functional: PBE 8, PBE+D3 9, r2SCAN 8, BEEF-vdW 2.

Audit outcomes: 85 screened, 118 without suitable relaxed references, 18 beyond
the existing |E_ads| <= 5 eV review screen, 5 incomplete/invalid/electronically
unconverged complex outputs, and 1,678 prepared inputs without OUTCARs. These are
snapshot counts; the batch was still running during extraction. Completed total
SPE energies remain in the full audit even where adsorption energies cannot be
published. Slab and molecular reference statuses are in the source JSON.

Energy convention matches the existing page:

```
E_ads(SPE) = final TOTEN(complex, NSW=0)
           - final TOTEN(relaxed clean slab)
           - final TOTEN(relaxed gas molecule)
```

Reuse the Kestrel extractor's completion, final electronic/ionic convergence,
INCAR/OUTCAR settings agreement, functional, POSCAR/OUTCAR composition, and
reference compatibility checks. The complex POSCAR SHA256 must match the recorded
original and staged POSCAR hashes; verify the source/destination, functional and
staged INCAR/KPOINTS hashes. No relaxation result is relabeled as single-point.
Both references must be converged ionic relaxations; missing gas or mismatched
slab references are not substituted or scaled.

The existing project policy permits POTCAR variants; names are retained in the
audit and no energy offset is applied. Positive absolute energies are allowed.
The 5 eV screen is a publication review threshold, not a convergence criterion.
Formula aliases are explicit; formic acid and formate stay distinct. For site
duplicates prefer the system without a site suffix, then lexical path, never
lowest energy. Prefer a matching slab in-plane cell and the original gas-molecule
spelling among compatible references; cell mismatch is not a new exclusion rule.

All published SPE energies and deltas are retained, including the 58 overlapping
Perlmutter candidates. The existing combined-column expression and relaxed cells
remain unchanged. New deltas use the displayed relaxed ML value. The descriptive
energy-decomposition plots and input checksums are regenerated from the updated
CSV; this does not claim a new ML/SPE reference-equivalence audit.

## Reproduce

Use Python with ASE, NumPy, BeautifulSoup4 and Matplotlib:

```bash
python scripts/extract_perlmutter_singlepoint.py --project-root /pscratch/sd/j/jcho5/VASP
python scripts/update_kestrel_singlepoint.py
python scripts/plot_energy_decomposition.py
python -m unittest discover -s scripts -p 'test_*.py'
```

`update_kestrel_singlepoint.py` now also loads the Perlmutter dataset. It fills gaps
only, with already-published data first, then Kestrel candidates, then Perlmutter
candidates. The relaxed-page updater reapplies this same merged SPE view.

Artifacts: `dft_perlmutter_singlepoint.csv` (85 selected results),
`dft_perlmutter_singlepoint_audit.csv` (all 1,904 candidates and component energies),
`dft_perlmutter_singlepoint_sources.json` (settings, composition, potentials,
convergence, reference paths, hashes and extraction timestamp), plus the updated
displayed CSV, HTML, decomposition data and plots.
