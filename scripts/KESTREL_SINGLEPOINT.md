# Kestrel single-point additions

The comparison page retains all 418 published SPE energies and their published
single-point ML differences, and fills 314 missing functional entries. It now
shows 732 SPE energies. The published source dataset remains unchanged; the
original displayed data are also archived as
`dft_comparison_singlepoint_published.csv`.

The single-point ML column is replaced by the requested literal expression:

```text
E_ads(ML) - E_ads(SPE) - E_ads(DFT, relaxed)
```

This uses the ML and DFT values displayed under relaxed comparison and the
full-precision SPE value, then rounds to three decimals. It is blank if any
input is missing. For example, Ag111 formate PBE gives
`-1.98 - (-3.2039268) - (-3.21) = +4.4339268 eV`, displayed as `+4.434`.
The existing single-point delta remains. Historical deltas retain their original
ML references, which sometimes differ from the relaxed comparison's ML values;
the download preserves those original references. New deltas use the displayed
relaxed-comparison ML value. Existing relaxed cells and all images are preserved.

## Sources and reference convention

Complexes are discovered under `poscar/best/<system>` and
`poscar/best/C<n>/<system>`, including the functional and `singlepoint/functional`
layouts. NSW is verified in both INCAR and OUTCAR; directory names alone never
establish single-point status.

```text
E_ads(SPE) = E(complex, NSW=0) - E(relaxed slab) - E(relaxed molecule)
```

References come from `vasp_mol`, `vasp_slab_kestrel`, and existing size-specific
or matching slabs in `vasp_slab`. No reference energy is scaled by atom count.
All 314 new entries select matching references from `vasp_slab`; the standard
`vasp_slab_kestrel` references were audited as well. Selection requires the same
surface, functional and composition. Molecule aliases use the project's
`mol_canon.py`; formate and formic acid remain distinct.

The extractor audited 393 reference calculations and 2,156 complex candidates.
706 candidates pass the checks, yielding 665 distinct screened SPE entries.
Of these, 314 fill gaps and 351 overlap preserved published values. Remaining
candidates: 1,227 unavailable/unconverged complex results, 100 without suitable
relaxed references, and 123 with adsorption energies outside the existing
5 eV review screen. Counts include duplicate directory layouts.

Requirements for new entries:

- Completed OUTCAR with final electronic convergence and finite final TOTEN.
- Complex INCAR and executed OUTCAR NSW both zero.
- Relaxed references with NSW greater than zero, active ionic optimization,
  electronic convergence, and the ionic convergence marker.
- INCAR/OUTCAR settings agreement, recognized and matching functional,
  and POSCAR/OUTCAR composition agreement.
- Exact complex composition equals slab plus molecule composition.
- Existing `|E_ads| <= 5 eV` publication screen. Larger values stay in the audit
  as review findings, not assertions of failed convergence.
- Potential names remain recorded, but POTCAR variants do not exclude results
  under the existing project policy. No absolute-energy sign check, gas-energy
  fallback, offset or correction is used.

For duplicate complexes, prefer `singlepoint`, then the shallowest path, then
lexical path; do not choose by energy. Slab references prefer matching in-plane
cells, then `vasp_slab_kestrel`, then lexical path. Cell differences are not a
new exclusion rule. Molecule references prefer the original molecule spelling.
The full selection rules, settings and component energies are recorded in the
source JSON. These are completion checks at the configured calculation settings,
not a numerical-parameter convergence study.

## Reproduce

Use Python with ASE, NumPy, Matplotlib and BeautifulSoup4:

```bash
python scripts/extract_kestrel_singlepoint.py \
  --project-root /kfs3/scratch/jcho5/goad-global-optimization
python scripts/update_kestrel_singlepoint.py
python -m unittest discover -s scripts -p 'test_*.py'
```

The calculation files are read only. All generated artifacts stay in the website
repository. No jobs are submitted. The existing relaxed-energy updater reapplies
this SPE layout when the Kestrel dataset is present, without changing its energy
filtering rules.

Artifacts:

- `dft_kestrel_singlepoint.csv`: all 665 selected screened results, component
  energies, NSW values and reference paths; overlapping historical SPE cells
  are not overwritten by these results.
- `dft_kestrel_singlepoint_audit.csv`: every examined complex candidate,
  selection flag and status, including rejected results and raw subtractions.
- `dft_kestrel_singlepoint_sources.json`: component settings, convergence,
  composition, potential identities, energies, source paths and file metadata.
- `dft_comparison_singlepoint.csv`: all displayed SPE energies, original ML
  references/deltas, the new expression and its relaxed-comparison inputs.
