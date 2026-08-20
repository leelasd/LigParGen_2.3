# BOSS vs. OpenMM/GROMACS single-point energy validation, and a documented residual

The reusable methodology, scripts, Dockerfiles, and full gotcha list for reproducing or extending this validation live in [`tools/energy_validation/`](../../tools/energy_validation/README.md) -- read that first if you're running this again rather than just reading about what it already found.

Replicated the original webserver's own validation practice -- generate parameters for real molecules, then check the output actually reproduces BOSS's own energy -- against 10 real reference Zmatrices from this machine's BOSS install (`molecules/small/`), fed through LigParGen via `-z` and compared to a from-scratch OpenMM single-point evaluation of the generated XML+PDB. Both sides are genuine single-point evaluations (BOSS's own `xSPM`, zero accepted/rejected moves; OpenMM's `Context.getState(getEnergy=True)` with no minimization) -- not an optimize-then-compare mismatch, which is its own well-known trap (see the cross-code discussion at [ParmEd#907](https://github.com/ParmEd/ParmEd/issues/907) and [openmm#1463](https://github.com/openmm/openmm/issues/1463), which also established that a total-energy match can hide compensating per-term errors, so every comparison here is decomposed by force group -- bond/angle/torsion/nonbonded -- not just totaled).

**Result: 8 of 10 match closely** (benzene, phenol, furan, aniline, acetamide, ammonia within BOSS's own ~2-decimal print-precision floor on its nonbonded term; methane and ethane's BOSS energy is captured correctly but can't be cross-checked against OpenMM at all -- BOSS's own `xSPM` writes an empty `plt.pdb` for these two files' older Zmatrix header style, a BOSS-side issue, not a LigParGen one).

**2 of 10 (toluene, anisole) show a real, unexplained nonbonded-term residual** (+0.13 and +0.27 kcal/mol respectively) that exhaustive investigation could not attribute to any LigParGen transcription bug: per-atom charges/LJ parameters match BOSS's own table exactly; the bond graph (and therefore OpenMM's auto-generated 1-4 exceptions) is correct; no missing torsions (verified both by hand-checking the actual dihedral geometry and by confirming OpenMM's class-based `<Proper>` matching already covers every physical instance of a rotatable-bond pattern, not just the one explicitly listed in the XML); no BOSS-side synonym/fallback parameter substitution in play; coordinate-precision truncation (BOSS's 5 internal decimals vs. the PDB writer's 3) is a real effect but two orders of magnitude too small (0.0011 kcal/mol) to be the cause. Both affected molecules share a structural motif the 8 matching ones don't: a ring carbon bonded to two symmetric ring neighbors plus one extra substituent. Best remaining guess -- unverified, and unverifiable from here -- is some difference in BOSS's own internal Fortran nonbonded summation specific to that motif; BOSS is a closed, proprietary binary, so this can't be chased past its own printed output.

Documented as a known, small (~2-4% relative) residual for substituted-ring molecules rather than treated as fixed. If it needs to move forward, it needs BOSS's own source or its maintainers, not anything on this side.

## Follow-up: ruling out "bad reference Zmatrix" as the cause

Sharpened one open question directly rather than leaving it as a guess: is the residual a property of this specific *molecular motif* (ring carbon with two symmetric neighbors plus a substituent), or an artifact of `toluen.z`/`anisol.z` specifically being old, possibly poorly-constructed reference files? There was good reason to suspect the latter -- inspecting `toluen.z`'s own BOSS-declared Fourier Coefficients table directly showed it has **zero ring-torsion parameters at all** (only 8 declared rows, all methyl-rotation-related, all zero-valued), unlike benzene's same table (24 rows, all real nonzero values) or a freshly-generated toluene's (36 rows, real `k2=15.167` ring torsions and `k2=10.46` planarity-enforcing impropers throughout). That's a genuine, confirmed data-quality gap in that one legacy file -- not a guess.

Tested it directly: submitted toluene by SMILES to the live production Space (`lsdodda/ligpargen`, the actual deployed app, not the local CLI), which builds its own fresh Zmatrix from scratch rather than reusing `toluen.z` at all. Pulled the resulting `.z` file back out of the Space's own output, re-ran it through BOSS's `xSPM` locally, and compared against a from-scratch OpenMM evaluation of the same run's XML+PDB, exactly as above:

| | BOSS | OpenMM | diff |
|---|---|---|---|
| bond | 0.2372 | 0.2339 | -0.003 |
| angle | 0.0677 | 0.0654 | -0.002 |
| torsion | 0.0011 | 0.0013 | +0.0002 |
| nonbonded | 3.60 | 3.7432 | **+0.143** |
| total | 3.9107 | 4.0438 | +0.133 |

Torsion now matches closely (both sides genuinely small and real, unlike the `toluen.z` case where it was trivially zero on both sides because nothing was declared) -- confirming the missing-parameter issue really was specific to that one legacy file, not a pipeline-wide problem. But the nonbonded gap is essentially unchanged (+0.143 vs. +0.129 with the old reference file) on a completely different geometry from a completely different code path. That rules out "bad reference Zmatrix" as the explanation for the nonbonded residual specifically, even though it was a real, separate, worth-fixing-if-anyone-revisits-that-file problem in its own right. The residual is reproducible across geometry sources, which points back toward the molecular motif itself (or BOSS's own nonbonded summation for it) rather than any one input file's quality.

## Extending to GROMACS: a real bug found and fixed

Extended the same methodology to `BOSS2GMX.py`'s `.itp`/`.gro` output (GROMACS is freely available via `apt-get install gromacs`, no license needed) to cross-check that the OpenMM-focused work above hadn't missed something writer-specific. It had -- a serious one, and not something introduced by this session's other changes.

**Bug found**: `boss2gmx()` gated writing the *entire* `[ dihedrals ]` section (both proper and improper torsions) on `len(tor_df.index) != len(full_tor.index)` -- i.e. only wrote anything at all if BOSS2GMX's own type-based torsion deduplication happened to remove at least one row. The deduplicated set (`tor_df`) was never actually used for writing either block; that length comparison only ever fed this gate. Confirmed directly for a freshly-generated toluene (36 real torsions, `k2=15.167` ring terms and `k2=10.46` planarity-enforcing impropers, all with distinct type-class names): `tor_df` and `full_tor` came out the same length, the condition was false, and the `.itp` had **no `[ dihedrals ]` section at all** -- not even an empty one. Not a rare case: phenol, benzene, and acetamide all hit the identical omission before the fix. Fixed by changing the gate to `len(full_tor.index) > 0` (write whenever there's anything to write, matching the "skip only when genuinely torsion-free" pattern already used elsewhere in this codebase). No other writer has the identical broken condition (checked directly via `grep`), so this was specific to `BOSS2GMX.py`.

**Post-fix energy comparison, toluene** (kcal/mol; same Space-generated molecule as the table above):

| | BOSS | OpenMM | GROMACS |
|---|---|---|---|
| bond | 0.2372 | 0.2339 | **0.3623** |
| angle | 0.0677 | 0.0654 | 0.0855 |
| torsion | 0.0011 | 0.0013 | 0.0085 |
| nonbonded | 3.60 | 3.7432 | 3.5111 |
| **total** | **3.9107** | 4.0438 | **3.9673** |

GROMACS's total actually lands closer to BOSS than OpenMM's does. But its bond term specifically shows a real ~0.125 kcal/mol gap, much larger than OpenMM's near-exact bond match. Per-bond force constants and equilibrium lengths were verified byte-identical between the `.xml` and the (fixed) `.itp` -- this isn't a data transcription bug. The most plausible explanation, not fully proven: the `.gro` coordinate format only stores 3 decimal places in *nanometers* (0.01 Å resolution), a full order of magnitude coarser than the `.pdb` format's 3 decimals in *Ångströms* (0.001 Å) that the OpenMM path reads -- and bond terms, having by far the largest force constants of any term, are the most sensitive to exactly this kind of positional truncation. Documented as a known, GROMACS-specific, format-precision-driven residual rather than a code bug -- see `tools/energy_validation/README.md`'s gotcha list for the full reasoning.
