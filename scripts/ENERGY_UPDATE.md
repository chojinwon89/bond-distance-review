# Comparison-page energy audit

The 2026-09-10 reference-compatibility correction checks 1,904 Perlmutter jobs.
206 pass the current checks: PBE 69, PBE+D3 64, r2SCAN 25, BEEF-vdW 48.
The page also retains 418 published single-point results across 133 systems;
those single-point source calculations have not been re-audited here.

## Why converged results can be excluded

`dft_jobs/CH3_Ag100/beef_vdw`, its slab, and its molecule all converged.
However, the complex OUTCAR reports `PAW_PBE Ag 02Apr2005`, while
`vasp_slab/Ag100/beef_vdw/OUTCAR` reports `PAW_PBE Ag_pv 09Dec2005`.
These use different valence configurations and cannot be mixed for adsorption
energy subtraction. Their raw difference is -26.05120794 eV. Positive component
total energies alone are not a reason to reject a VASP calculation.

The full audit finds 917 pseudopotential mismatches (some also have atom-count
or convergence issues). In BEEF-vdW, 204 entries have potential mismatches;
188 of those have all three calculations converged. The same mismatch occurs
for Ag/Ag_pv, Cu/Cu_pv, Pd/Pd_pv, Pt/Pt_pv and Rh/Rh_pv.

The previous audit checked convergence and metal counts but missed POTCAR
identity. 459 of its 665 verified entries fail the new compatibility check.
The page now withholds incompatible or unverified relaxed energies and their
ML differences, displaying an explicit status instead of `ref broken` or an
old numerical value. Status tooltips give reference details and convergence.
Historical values remain in Git history and the original published datasets.
No VASP inputs, outputs, reference energies or jobs are modified by this audit.

## Checks and limitations

- Final free-energy TOTEN (eV), not sigma-to-zero energy.
- E_ads = complex - clean slab - gas molecule; no slab scaling or offsets.
- Completion and final electronic convergence for every component, and ionic
  convergence when NSW > 0.
- Matching OUTCAR TITEL potential identities for each element between complex
  and references; matching slab metal counts and full component composition.
  TITEL is a recorded identity check, not a byte-level POTCAR hash comparison.
- Missing metadata or component energies exclude the result.
- Absolute component energy signs are unrestricted. |E_ads| > 5 eV is retained
  as an energy-review flag, not labeled failed convergence or a broken file.
- The check does not establish full convergence with respect to numerical
  settings or enforce identical slab lattice vectors. Geometry images and
  distances are from the earlier extraction.
- Multiple sites use the lowest verified energy. Explicit molecule aliases
  are in the update script; single-point joins retain exact dataset labels.

## Reproduce (read-only calculation access)

```bash
python3 scripts/extract_perlmutter_energies.py --vasp-root /pscratch/sd/j/jcho5/VASP
python3 scripts/update_comparison_energies.py
python3 -m unittest discover -s scripts -p 'test_*.py'
```

The extractor imports only the project's directory-discovery and naming helpers
from `calc_binding_energy.py`; the website's audit owns its validation rules.
It reads OUTCAR header/tail data and caches shared references. OUTCAR.gz is
supported when the uncompressed file is absent.

Outputs:

- `dft_perlmutter_energy_audit.csv`: all jobs, component and raw energies,
  convergence, statuses, OUTCAR paths, potential identities and atom counts.
  `E_ads_raw` is diagnostic only, including incompatible subtractions.
- `dft_comparison_perlmutter.csv`: only energies passing the current checks.
- `dft_comparison_singlepoint.csv`: displayed published single-point values.
- `dft_comparison.html`: comparison tables, diagnostic labels and audit links.
