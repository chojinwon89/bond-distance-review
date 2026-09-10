# Energy decomposition

Run `python scripts/plot_energy_decomposition.py` after refreshing comparison
energies. Requires NumPy and Matplotlib. It reads the displayed single-point CSV
and never changes source energies, filtering, images or existing table cells.

For adsorption energies M (ML), S (SPE), R (relaxed DFT):

```
A = M - S
B = S - R
T = M - R = A + B
T - A = B
```

M-S-R is not a corrected error. Even a perfect M=S=R=-1 eV gives +1 eV
under that expression, instead of zero. The legacy column is retained only for
publication continuity and excluded from this analysis.

Use one M consistently, not historical deltas with potentially different ML
references. No energy thresholds are added. Complete finite triples alone enter
the distributions, scatter plots and summary. Cross-functional shifts use only
systems with SPE values for all four functionals, without requiring ML or relaxed
energies for this comparison. These are descriptive
sample statistics, not uncertainty estimates or independent observations.

Interpret A as same-geometry model error only after verifying ML/SPE geometry
and reference conventions. To isolate functional effects, evaluate different
functionals on identical complex, molecule and slab structures using compatible
potentials and converged settings. Different-functional adsorption energies
alone do not isolate exchange-correlation differences if other settings vary.

Interpret B as relaxation lowering only after verifying linked initial/final
structures, identical functional/settings and identical adsorption references.
Relaxing molecule/slab references as well adds their deformation contributions;
it is no longer just complex relaxation. New reference SPEs must be explicitly
joined and audited before replacing reference energies. A negative B is flagged
descriptively at -0.01 eV, never removed or silently corrected.

Current data mix historical and audited results and allow POTCAR variants under
the publication policy. No new geometry/reference equivalence audit is claimed.
Published ML/relaxed rounding limits precision. Missingness and the existing
energy screen affect coverage. Formate and formic acid are never merged.

Generated artifacts include PNG plots, row-level data, statistics, and metadata
with an input SHA256 checksum. The section is idempotently refreshed in
`dft_comparison.html#energy-decomposition`.
