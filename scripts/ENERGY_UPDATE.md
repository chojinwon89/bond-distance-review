# Comparison-page energy refresh

The 2026-09-09 refresh adds the existing published single-point dataset to the
415 structure cards using exact (surface, molecule, functional) keys: 418
functional results across 133 systems. Unmatched single-point rows remain
available in the full source CSV; missing matches display a dash.

The 2026-09-10 Perlmutter extraction covers 1,904 functional jobs. Of these, 665 passed the
project collector's reference checks and additional final electronic/ionic
convergence checks for the complex, slab, and gas molecule. Relative to the
published `dft_vs_mlip_pairs.csv`, 94 are additional functional entries and 26
have changed energies. Previously published entries that cannot currently be
verified are retained with explicit tooltips, rather than replaced by values
from incomplete calculations. Images and geometry measurements are the earlier
published extraction; this update concerns energies.

The 2026-09-10 refresh adds the newly verified N2/Ir111 BEEF-vdW result
(-0.43191224 eV). The previous 664 verified energies are unchanged. Twenty
audit rows changed as calculations progressed; incomplete results remain
excluded from the verified dataset.

Energy convention: final VASP free-energy TOTEN in eV, as in the existing
`calc_binding_energy.py`; not sigma-to-zero energy. Adsorption energy is complex
minus clean slab minus gas molecule. The collector checks matching metal atom
counts, negative component total energies, and |E_ads| <= 5 eV. These are the
existing workflow's screening rules, not universal physical bounds. No slab
energy scaling is performed. Same-functional reference paths and convergence
statuses are included in the audit CSV. Multiple sites, if present, use the
lowest verified adsorption energy. Explicit formula/name aliases are in the
page-update script; single-point joins do not alias distinct dataset labels.

Reproduce on Perlmutter (read-only access to calculations):

```bash
python3 scripts/extract_perlmutter_energies.py --vasp-root /pscratch/sd/j/jcho5/VASP
python3 scripts/update_comparison_energies.py
```

The extractor loads `calc_binding_energy.py` from the supplied VASP project.
The extraction used the local 2026-08-28 version. Outputs:

- `dft_perlmutter_energy_audit.csv`: all extracted jobs and their screening status.
- `dft_comparison_perlmutter.csv`: verified relaxed energies displayed on the page.
- `dft_comparison_singlepoint.csv`: displayed single-point values.
- `dft_comparison.html`: refreshed tables, download links, and provenance note.

The published single-point source CSV and other website pages are unchanged.
