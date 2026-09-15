# Shared molecule notation and energy availability

`molecule_names.py` is the shared alias registry for relaxed energies, SPE
energies, reference lookup, coverage reports, and the comparison-page search.
`molecule_aliases.csv` publishes the exact equivalences. The 85 registry entries
include names as well as aliases; Unicode subscript digits normalize to ASCII.

Examples: DME = CH3OCH3 = dimethyl ether; ethanol = CH3CH2OH = C2H5OH;
methanol = CH3OH; acetaldehyde = CH3CHO; acetic acid = CH3COOH;
formaldehyde = H2CO; formic acid = HCOOH; ethene = ethylene = C2H4.
Molecular formulas alone do not establish identity: C2H6O is not automatically
assigned to ethanol or DME. Formate/formic acid, methoxy/hydroxymethyl, and
propanol/isopropanol remain distinct.

The page shows condensed formulas in headings and the molecule selector. An
exact recognized alias matches that molecule only; other queries retain partial
surface/name search. In particular, `ethanol` does not match `methanol`, and
`CH3OCH3`, `CH₃OCH₃`, `DME`, and `dimethyl ether` return the same DME cards.
`update_molecule_notation.py` preserves all energy tables and structure images.

`molecule_name_audit.csv` maps 613 source directory names to their canonical names
and indicates whether the molecule appears on the current comparison page.
Relaxed reference lookup preserves a usable original reference. When it is
missing or failed, it can use a converged equivalent-name reference with the
same functional and composition. Reference selection never substitutes an
isomer or chooses an alias by its energy. SPE reference lookup uses the same
canonical identities and retains its geometry/provenance checks.

## DME findings at the 2026-09-15 refresh

The relaxed `dft_jobs/CH3OCH3_*/beef_vdw` complexes are converged. Their unusual
binding energies remain visible with amber review markers wherever the slab and
molecule references are converged and composition matched. Magnitude alone does
not suppress a completed binding energy. These marked values remain outside
the figure statistics, and previously published values are preserved.

Some remaining DME gaps have separate causes:

- The DME/CH3OCH3 PBE gas outputs stop with `ZBRENT: fatal error in bracketing`.
  Neither spelling provides a converged alternative; their final total energy
  is available as a diagnostic in the component inventory.
- The PBE+D3 gas outputs are also unfinished.
- Several surfaces have a different number of metal atoms in the complex and
  the available clean slab. A direct subtraction is not a composition-matched
  adsorption energy.
- Some SPE complexes have not converged or have not run yet.

`dme_energy_status.csv` contains both relaxed and SPE DME candidates, raw energy
subtractions where available, convergence/status, and source/reference paths.
An unusual converged value is distinguished from a failed or unmatched reference.

Run `python scripts/update_molecule_notation.py` after the existing energy,
figure, and coverage refresh commands. Tests cover alias equivalence, isomer
separation, reference preference, and idempotent notation updates without energy
changes. The refreshed SPE batch supplies three further table entries, bringing
the displayed total to 896.
