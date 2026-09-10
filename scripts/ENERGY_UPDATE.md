# Comparison-page energy audit

The current project policy ignores POTCAR identity differences when screening
results, including Ag versus Ag_pv. Actual potential identities remain in the
audit. This is a requested filtering assumption, not a claim of physical
energy-reference equivalence. No energy offset or correction is applied.

The audit covers 1,904 jobs. 665 pass the remaining checks: PBE 240, PBE+D3 227,
r2SCAN 150 and BEEF-vdW 48. This restores 459 entries excluded by the earlier
potential-identity check. The page also retains 418 published single-point
results across 133 systems; those source calculations have not been re-audited.

CH3/Ag100 BEEF-vdW is no longer excluded for Ag versus Ag_pv. Its raw energy
subtraction is -26.05120794 eV, so it still displays "Energy needs review"
because the independent |E_ads| <= 5 eV screen remains enabled. All three
component calculations converged. Raw and component energies are downloadable
in the audit even when excluded from the screened comparison.

No VASP inputs, outputs, reference energies or jobs are modified by this audit.

## Checks and limitations

- Final free-energy TOTEN (eV), not sigma-to-zero energy.
- E_ads = complex - clean slab - gas molecule; no slab scaling or offsets.
- Completion and final electronic convergence for every component, and ionic
  convergence when NSW > 0.
- Matching slab metal counts and full component composition. OUTCAR TITEL
  potential identities are recorded but differences do not exclude results.
- Missing metadata or component energies exclude the result.
- Absolute component energy signs are unrestricted. |E_ads| > 5 eV is retained
  as an energy-review flag, not labeled failed convergence or a broken file.
- The check does not establish full convergence with respect to numerical
  settings or enforce identical slab lattice vectors. Geometry images and
  distances are from the earlier extraction.
- Multiple sites use the lowest screened energy. Explicit molecule aliases
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
  `E_ads_raw` is diagnostic only, including out-of-range subtractions.
- `dft_comparison_perlmutter.csv`: only energies passing the current checks.
- `dft_comparison_singlepoint.csv`: displayed published single-point values.
- `dft_comparison.html`: comparison tables, diagnostic labels and audit links.
