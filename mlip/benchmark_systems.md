# MLIP geometry benchmark -- curated systems

28 systems. Columns: surface, adsorbate, expected site, key metal-adsorbate bond, internal bond, method, source.


| id | surface | ads | ref site | d(M-X) | internal | method | source |
|----|---------|-----|----------|--------|----------|--------|--------|
| CO_Pt111_atop | Pt(111) | CO | atop | Pt-C 1.85 | C-O 1.15 | LEED expt + PBE | Ogletree 1986 SurfSci 173,351; Feibelman 2001 JPCB 105,4018 |
| CO_Pt111_fcc | Pt(111) | CO | fcc | Pt-C 2.05 | C-O 1.16 | PBE | Feibelman 2001 JPCB 105,4018 |
| CO_Cu111_atop | Cu(111) | CO | atop | Cu-C 1.83 | C-O 1.15 | expt/PBE | Hollins 1983; DFT |
| CO_Pd111_fcc | Pd(111) | CO | fcc hollow | Pd-C 2.05 | C-O 1.17 | PBE | Loffreda 1998; Honkala DFT |
| O_Pt111_fcc | Pt(111) | O | fcc hollow | Pt-O 2.0 | - | PBE | Todorova 2006; Hafner DFT |
| H_Pt111_fcc | Pt(111) | H | fcc hollow | Pt-H 1.87 | - | PBE | Olsen 1999; Kresse DFT |
| NO_Pt111_fcc | Pt(111) | NO | fcc hollow | Pt-N 1.95 | N-O 1.18 | PBE | Aizawa 1999; DFT |
| NH3_Pt111_top | Pt(111) | NH3 | atop (N lone pair down) | Pt-N 2.1 | - | PBE | Offermans 2007; DFT |
| H2O_Pt111_top | Pt(111) | H2O | atop, near-flat | Pt-O 2.4 | - | PBE+vdW | Michaelides 2003; Carrasco 2012 |
| N2_Pt111_top | Pt(111) | N2 | atop, end-on / weak | Pt-N 2.1 | N-N 1.1 | indicative | NIST CCCBDB (N-N 1.098); indicative surface geometry |
| C2H4_Pt111_bri | Pt(111) | C2H4 | di-sigma (bridge) | Pt-C 2.1 | C-C 1.49 | PBE | Watson 2001; Cremer 1996 |
| CO2_Cu111_top | Cu(111) | CO2 | atop / physisorbed | Cu-O 2.3 | C-O 1.16 | indicative | NIST CCCBDB (C=O 1.162); indicative surface geometry |
| CH4_Pt111_top | Pt(111) | CH4 | atop / physisorbed | Pt-C 3.1 | C-H 1.09 | indicative | NIST CCCBDB (C-H 1.087); indicative surface geometry |
| CH3_Pt111_top | Pt(111) | CH3 | atop | Pt-C 2.08 | C-H 1.09 | indicative | NIST CCCBDB (C-H 1.079); indicative surface geometry |
| OH_Pt111_top | Pt(111) | OH | atop / bridge | Pt-O 2.0 | O-H 0.97 | indicative | NIST CCCBDB (O-H 0.970); indicative surface geometry |
| CH3OH_Cu111_top | Cu(111) | CH3OH | atop (O lone pair down) | Cu-O 2.1 | C-O 1.43 | indicative | NIST CCCBDB (C-O 1.428); indicative surface geometry |
| H2CO_Pt111_top | Pt(111) | H2CO | atop (eta1-O) / eta2 | Pt-O 2.1 | C-O 1.21 | indicative | NIST CCCBDB (C=O 1.208); indicative surface geometry |
| HCOOH_Cu111_top | Cu(111) | HCOOH | atop (carbonyl O down) | Cu-O 2.2 | - | indicative | NIST CCCBDB geometry; indicative surface geometry |
| CH3O_Cu111_fcc | Cu(111) | CH3O | fcc hollow (O down) | Cu-O 1.95 | C-O 1.42 | indicative | NIST CCCBDB (C-O ~1.42); indicative surface geometry |
| C2H6_Pt111_top | Pt(111) | C2H6 | atop / physisorbed | Pt-C 3.2 | C-C 1.53 | indicative | NIST CCCBDB (C-C 1.535); indicative surface geometry |
| C2H2_Pd111_bri | Pd(111) | C2H2 | di-sigma / bridge | Pd-C 2.05 | C-C 1.32 | indicative | NIST CCCBDB (C#C 1.203); indicative surface geometry |
| CH3CHO_Pt111_top | Pt(111) | CH3CHO | atop (eta1-O) | Pt-O 2.1 | C-O 1.22 | indicative | NIST CCCBDB (C=O 1.216); indicative surface geometry |
| CH3COOH_Pt111_top | Pt(111) | CH3COOH | atop (carbonyl O down) | Pt-O 2.1 | C-C 1.5 | indicative | NIST CCCBDB (C-C 1.50); indicative surface geometry |
| C2H5OH_Pt111_top | Pt(111) | CH3CH2OH | atop (O lone pair down) | Pt-O 2.2 | C-O 1.43 | indicative | NIST CCCBDB (C-O 1.431); indicative surface geometry |
| DME_Pt111_top | Pt(111) | CH3OCH3 | atop (ether O down) | Pt-O 2.3 | C-O 1.41 | indicative | NIST CCCBDB (C-O 1.410); indicative surface geometry |
| H2S_Pt111_top | Pt(111) | H2S | atop (S down) | Pt-S 2.3 | S-H 1.34 | indicative | NIST CCCBDB (S-H 1.336); indicative surface geometry |
| SO2_Pt111_top | Pt(111) | SO2 | atop / eta2 (S down) | Pt-S 2.2 | S-O 1.44 | indicative | NIST CCCBDB (S-O 1.431); indicative surface geometry |
| HCN_Pt111_top | Pt(111) | HCN | atop, end-on (N down) | Pt-N 2.0 | C-N 1.16 | indicative | NIST CCCBDB (C#N 1.156); indicative surface geometry |

## Notes per system

- **CO_Pt111_atop** (Pt(111), CO): CO/Pt(111) SITE PUZZLE: expt=atop, plain PBE wrongly prefers fcc hollow -- a classic MLIP stress test
- **CO_Pt111_fcc** (Pt(111), CO): same system started in the fcc hollow; compare final energy vs the atop start to read out the MLIP's site preference
- **CO_Cu111_atop** (Cu(111), CO): CO binds atop Cu(111), weak-moderate
- **CO_Pd111_fcc** (Pd(111), CO): Pd(111) favours the fcc hollow at low coverage
- **O_Pt111_fcc** (Pt(111), O): atomic O in the fcc hollow, O-Pt ~2.0
- **H_Pt111_fcc** (Pt(111), H): light atom, very flat PES; fcc~atop near-degenerate. Pt-H(hollow)~1.87, Pt-H(atop)~1.55
- **NO_Pt111_fcc** (Pt(111), NO): N-down NO in fcc hollow at low coverage
- **NH3_Pt111_top** (Pt(111), NH3): molecular NH3 datives through N atop
- **H2O_Pt111_top** (Pt(111), H2O): weak; O nearly atop, molecular plane ~parallel
- **N2_Pt111_top** (Pt(111), N2): N2 is a WEAK binder; on non-magnetic Pt it adsorbs end-on only weakly (or desorbs) -- a long distance can be CORRECT. Internal N-N (1.098) should be preserved
- **C2H4_Pt111_bri** (Pt(111), C2H4): rehybridisation test: C=C (1.33) stretches toward 1.49 as it binds di-sigma
- **CO2_Cu111_top** (Cu(111), CO2): CO2 physisorbs weakly on Cu(111) and stays near-linear; it only bends to CO2^d- when activated -- a large height here can be CORRECT
- **CH4_Pt111_top** (Pt(111), CH4): methane physisorbs; expect ~gas-phase geometry at a large height (a big d(M-C) is correct, not a failure)
- **CH3_Pt111_top** (Pt(111), CH3): methyl radical chemisorbs C-down atop Pt(111), Pt-C ~2.1
- **OH_Pt111_top** (Pt(111), OH): hydroxyl binds O-down, Pt-O ~2.0; O-H preserved near 0.97
- **CH3OH_Cu111_top** (Cu(111), CH3OH): methanol datives through O; weak, molecular, C-O ~1.43 preserved
- **H2CO_Pt111_top** (Pt(111), H2CO): formaldehyde; carbonyl C=O (1.21) elongates if it binds eta2-CO
- **HCOOH_Cu111_top** (Cu(111), HCOOH): molecular formic acid; may deprotonate to formate (classic on Cu)
- **CH3O_Cu111_fcc** (Cu(111), CH3O): methoxy binds O-down in a hollow on Cu(111); a textbook oxygenate fragment
- **C2H6_Pt111_top** (Pt(111), C2H6): ethane physisorbs; large height + gas-phase C-C (1.535) is the correct answer
- **C2H2_Pd111_bri** (Pd(111), C2H2): acetylene binds di-sigma; the C#C triple bond (1.203) stretches toward ~1.32
- **CH3CHO_Pt111_top** (Pt(111), CH3CHO): acetaldehyde binds through the carbonyl O; C=O ~1.22
- **CH3COOH_Pt111_top** (Pt(111), CH3COOH): molecular acetic acid; may deprotonate to acetate on more reactive metals
- **C2H5OH_Pt111_top** (Pt(111), CH3CH2OH): ethanol datives through O; weak, molecular
- **DME_Pt111_top** (Pt(111), CH3OCH3): dimethyl ether; ether O to the surface, weakly bound
- **H2S_Pt111_top** (Pt(111), H2S): H2S binds S-down atop Pt; can dissociate to SH+H on reactive metals
- **SO2_Pt111_top** (Pt(111), SO2): SO2 binds S-down or eta2-S,O; S-O ~1.43
- **HCN_Pt111_top** (Pt(111), HCN): hydrogen cyanide adsorbs end-on through N; C#N ~1.16

EMT-runnable subset (no torch needed): CO_Pt111_atop, CO_Pt111_fcc, CO_Cu111_atop, CO_Pd111_fcc, O_Pt111_fcc, H_Pt111_fcc, NO_Pt111_fcc, NH3_Pt111_top, H2O_Pt111_top, N2_Pt111_top, C2H4_Pt111_bri, CO2_Cu111_top, CH4_Pt111_top, CH3_Pt111_top, OH_Pt111_top, CH3OH_Cu111_top, H2CO_Pt111_top, HCOOH_Cu111_top, CH3O_Cu111_fcc, C2H6_Pt111_top, C2H2_Pd111_bri, CH3CHO_Pt111_top, CH3COOH_Pt111_top, C2H5OH_Pt111_top, DME_Pt111_top, HCN_Pt111_top