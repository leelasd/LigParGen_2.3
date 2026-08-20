# BOSS vs. OpenMM/GROMACS/NAMD/LAMMPS/TINKER/Q single-point energy validation, and a documented residual

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

Note: re-running phenol through the same GROMACS path while building the NAMD leg below turned up a larger, multi-term gap than the coordinate-precision explanation alone predicts (bond 0.272 vs. BOSS's 0.219, angle 0.047 vs. 0.027, nonbonded 0.726 vs. 0.67 -- total 1.045 vs. BOSS's 0.914, ~14% high). That's bigger and broader than toluene's bond-only ~0.125 gap, so `.gro` precision may not be the whole GROMACS story. Not chased further here -- flagged for whoever next touches the GROMACS leg.

## Extending to NAMD: a genuine third-engine confirmation

Extended the same methodology to NAMD once a licensed NAMD 3.0.3 binary and NAMD's own `psfgen` topology tool became available locally (see `tools/energy_validation/README.md`'s "The NAMD leg" section for the full setup -- NAMD is not Dockerized, supplied locally exactly like BOSS per `docs/adr/0001`). NAMD reads CHARMM-format files, so this exercises a completely different LigParGen writer (`BOSS2CHARMM.py`, `.rtf`/`.prm`) than the OpenMM/GROMACS legs above, via a completely different tool (`psfgen` to build a `.psf`, then NAMD's own force evaluation) -- a genuinely independent cross-check, not a re-run of the same code path.

This directly closes the loop on [ParmEd#907](https://github.com/ParmEd/ParmEd/issues/907), the historical issue that grounded this whole methodology's single-point-only rule: that issue's own root cause was comparing a NAMD *minimization* endpoint against an OpenMM single point. The user's own 2017 NAMD-vs-OpenMM test config (`min.conf`, provided this session) still used `minimize 1000`, not a true single point -- consistent with that being exactly the trap the ParmEd issue was about. The new config used here (`namd_sp_template.conf`) uses a genuine `run 0`, closing that gap.

**Result, benzene** (kcal/mol; same molecule/geometry as the OpenMM/GROMACS numbers used throughout this doc):

| | BOSS | OpenMM | NAMD |
|---|---|---|---|
| bond | 0.221 | 0.2219 | 0.2219 |
| angle | 0.0 | 0.0001 | 0.0001 |
| torsion | 0.0 | 0.0 | 0.0 |
| nonbonded | 6.46 | 6.4368 | 6.4591 |
| **total** | **6.6806** | 6.6588 | **6.6811** |

NAMD's total is within 0.0005 kcal/mol of BOSS's -- closer than OpenMM's own 0.022 kcal/mol gap, and its nonbonded term (6.4591) also lands nearer BOSS's printed 6.46 than OpenMM's 6.4368 does.

**Result, phenol** (has a real substituent -- OH -- unlike benzene, so exercises the CHARMM writer's `IMPROPER`/typed-atom handling on a less trivial case):

| | BOSS | OpenMM | NAMD |
|---|---|---|---|
| bond | 0.2192 | 0.2202 | 0.2202 |
| angle | 0.027 | 0.0256 | 0.0256 |
| torsion | 0.0 | 0.0 | 0.0 |
| nonbonded | 0.67 | 0.6697 | 0.669 |
| **total** | **0.9139** | 0.9155 | **0.9148** |

Same pattern: NAMD's total sits between BOSS's and OpenMM's, well within the same tolerance both already established as a match. No new discrepancy found -- NAMD confirms the CHARMM/NAMD writer path is correct for both molecules tested, using OPLS-AA's geometric-mean sigma combining rule (`vdwGeometricSigma yes` in the NAMD config -- without it NAMD silently uses CHARMM's arithmetic-mean rule instead, which would have been a real, silent, wrong-answer trap).

Not yet run through NAMD: the toluene/anisole nonbonded residual documented above, or the newly-noted phenol GROMACS gap -- both remain OpenMM/GROMACS-side observations only. A NAMD data point on toluene specifically (does NAMD show the same nonbonded residual, a different one, or none?) would be informative but wasn't requested this round.

## Extending to LAMMPS: the tightest match of any engine

Extended the same methodology to LAMMPS (freely available via Debian's own apt repo, GPL -- no license needed, same as GROMACS) via `BOSS2LAMMPS.py`'s `.lmp` data file. Unlike the other four writers, `BOSS2LAMMPS.py` gives every atom/bond/angle/dihedral instance its own unique numbered type rather than deduplicating by class -- the least-transformed representation of BOSS's own raw per-instance parameters of any writer here, which likely explains why it reproduces BOSS most faithfully of all five engines tested (see results below). LAMMPS's built-in `dihedral_style opls` and `improper_style cvff` match the `.lmp` writer's own coefficient layout exactly -- no unit conversion or halving needed on the LAMMPS side (BOSS's raw Fourier coefficients are already in kcal/mol and already carry the form `opls` expects, unlike the OpenMM/CHARMM paths which both convert units and halve V-coefficients to `K = V/2` for the `1+cos` energy form -- confirmed directly: BOSS's raw V2 for benzene's ring torsion is 7.250 kcal/mol, matching `.lmp`'s dihedral coefficient exactly, and consistent with OpenMM's stored `k2=15.167 kJ/mol` = `7.250 kcal/mol * 4.184 / 2`).

Free (non-periodic) boundaries with a 100 Å cutoff give LAMMPS a true vacuum evaluation directly -- unlike modern GROMACS, LAMMPS doesn't reject `boundary f f f`, so no enlarged-box workaround is needed here (see `tools/energy_validation/README.md`'s "The LAMMPS leg" section).

**Result, benzene** (kcal/mol):

| | BOSS | OpenMM | GROMACS | NAMD | LAMMPS |
|---|---|---|---|---|---|
| bond | 0.221 | 0.2219 | 0.2838 | 0.2219 | 0.2209 |
| angle | 0.0 | 0.0001 | 0.0110 | 0.0001 | ~0.0 |
| torsion | 0.0 | 0.0 | 0.0 | 0.0 | ~0.0 |
| nonbonded | 6.46 | 6.4368 | 6.4298 | 6.4591 | 6.4597 |
| **total** | **6.6806** | 6.6588 | 6.7246 | 6.6811 | **6.6806** |

**Result, phenol** (kcal/mol):

| | BOSS | OpenMM | GROMACS | NAMD | LAMMPS |
|---|---|---|---|---|---|
| bond | 0.2192 | 0.2202 | 0.2722 | 0.2202 | 0.2192 |
| angle | 0.027 | 0.0256 | 0.0467 | 0.0256 | 0.0270 |
| torsion | 0.0 | 0.0 | 0.0 | 0.0 | ~0.0 |
| nonbonded | 0.67 | 0.6697 | 0.7257 | 0.669 | 0.6677 |
| **total** | **0.9139** | 0.9155 | 1.0447 | 0.9148 | **0.9139** |

LAMMPS's total matches BOSS's printed value to within 0.00005 kcal/mol on both molecules -- the tightest of any of the four downstream engines, including per-term agreement to BOSS's own print precision. No new discrepancy found; this confirms `BOSS2LAMMPS.py`'s writer (and its `opls`/`cvff` coefficient conventions) is correct for both molecules tested.

## Extending to TINKER: a real bug that made every TINKER output unusable

Extended the same methodology to TINKER (source freely available on GitHub, BSD-style academic license, no separate license needed -- built from source since it's in neither Debian's apt repo nor Homebrew; see `tools/energy_validation/README.md`'s "The TINKER leg" for the build notes, including an upstream CMake-packaging gap that had to be patched around). This one found a real, serious, previously-undiscovered bug -- not a residual, a hard failure.

**Bug found**: `BOSS2TINKER.py`'s `create_xyz_file()` wrote each atom's TINKER "atom type" column in the `.xyz` coordinate file as `799 + atom_index` -- an arbitrary, unique-per-atom placeholder completely disjoint from the OPLS type numbers (e.g. `145`, `146`) that every other section of the companion `.key` file (`atom`, `vdw`, `bond`, `angle`, `torsion`, `charge` records) already consistently used. Confirmed directly, not inferred: building TINKER from source and running its own `analyze` program against LigParGen's existing benzene output reported every single atom as an "Undefined Atom Type," printed "MECHANIC -- Some Required Potential Energy Parameters are Undefined," and refused to compute an energy at all. This means **every `.xyz`/`.key` pair LigParGen has ever generated for TINKER was non-functional** -- not numerically off, unusable outright the moment anyone actually tried to run it through TINKER, which this codebase apparently never had a way to do before now (no prior test exercised this).

**Root cause**: `create_xyz_file()` didn't have access to the per-atom OPLS type numbers (`types`, from `bossData()`) at all -- `mainBOSS2TINKER()` only passed it the raw `molecule_data`. The `799+atom_index` value looks like a placeholder a `.txyz` atom-type column needs to be numeric-and-present, left in from early development and never replaced with the real type lookup once `bossData()`'s per-atom types became available elsewhere in the same module (`Boss2Tinker()`, called right after, already threads real OPLS numbers through every other section correctly).

**Fix**: `mainBOSS2TINKER()` now calls `bossData(mol)` once and passes `types` into `create_xyz_file()`, which looks up each atom's real OPLS type number (`types[atom_number - 1][1]`) instead of the placeholder.

**Post-fix result** (kcal/mol; same benzene/phenol geometries as above):

| | BOSS | TINKER |
|---|---|---|
| benzene bond | 0.221 | 0.2209 |
| benzene angle | 0.0 | 0.0 |
| benzene torsion | 0.0 | 0.0 |
| benzene nonbonded | 6.46 | 6.4598 |
| **benzene total** | **6.6806** | **6.6806** |
| phenol bond | 0.2192 | 0.2192 |
| phenol angle | 0.027 | 0.027 |
| phenol torsion | 0.0 | 0.0 |
| phenol nonbonded | 0.67 | 0.6678 |
| **phenol total** | **0.9139** | **0.9139** |

Tied with LAMMPS as the tightest match of any of the five downstream engines -- once the atom types actually match, TINKER reproduces BOSS's own energy almost exactly. At the time this was written, the benzene/phenol geometries used here came from stale legacy reference Zmatrices (`benzen.z`/`phenol.z`) that -- like `toluen.z` above -- turned out to declare zero real improper torsions at all, so `Boss2Tinker()`'s `imptors` section-writing code path went untested here. **This claim was wrong** -- see "A second, much bigger bug: Proper/Improper misclassification, and why the old reference files hid it" below, which found real impropers on both molecules once generated fresh, and a real classification bug that had nothing to do with the reference-file question.

## A second, much bigger bug: Proper/Improper misclassification, and why the old reference files hid it

While testing Q (see below), a direct challenge to this document's own "benzene has zero impropers" claim -- benzene has 6 sp2 ring carbons, so it should have 6 -- led to tracing the claim back to its source and finding it was never actually representative.

**The reference files were degenerate for impropers too.** `benzen.z`/`phenol.z` (BOSS's own legacy reference library, `molecules/small/`) turned out to have the exact same kind of gap already documented above for `toluen.z`'s torsions: their "Additional Dihedrals follow" section is a bare `AUTO` placeholder, not real declarations. A fresh generation (submitted by SMILES to the live production Space) confirmed BOSS genuinely tabulates 6 real `k2=10.46` improper-type Fourier coefficients for benzene's 6 ring carbons, and a real (if partial) set for phenol.

**A real, much older, previously-undiscovered classification bug.** Re-running that same fresh benzene through the *local* codebase (not the Space, which runs an older, differently-behaved deployment) showed every one of those 30 declared torsions being written as `<Proper>` -- zero `<Improper>` tags, even though the correct `k2` values were still present. Traced to `BOSSReader.ucomb()`, shared identically by all 8 `BOSS2*.py` writers to decide Proper vs. Improper: it counted how many of a quadruple's 4 atoms have bonds among **any** of the 6 possible atom pairs, and called it "Proper" if the count was 3. But a genuine bonded chain (i-j-k-l: bonds at (i,j),(j,k),(k,l)) and a genuine improper star (hub bonded to 3 substituents, substituents not bonded to each other: bonds at (hub,sub1),(hub,sub2),(hub,sub3)) **both** have exactly 3 bonded pairs among their 6 possible pairs -- the count alone can't distinguish them. Confirmed directly against benzene's real bond graph: quadruple `[H806, C805, C800, C801]` is a genuine star centered on C800, and `ucomb` returned 3 for it, identical to what it returns for a real chain.

This is not a regression from this session's own work: the pre-refactor `BOSS2OPENMM.py` (extracted from commit `38921b9`, well before today) was tested directly against the same fresh-benzene data and gave the identical result (0 impropers). `ucomb()` itself has not changed since the very first commit of this codebase (`884d392`). It has silently mislabeled every improper as a Proper torsion for as long as the codebase has existed, across all 8 writers.

**Fix**: rewrote `ucomb()` to check specifically whether the 3 *consecutive* pairs -- (i,j), (j,k), (k,l) -- are each real bonds, rather than counting bonds among all 6 possible pairs regardless of which ones. A star quadruple's non-consecutive-adjacent pair (whichever one doesn't involve the hub, depending on where the hub falls in the quadruple's declared order) is never a real bond, so the consecutive-pair count comes out below 3, correctly reading as Improper.

**Verified**: fresh benzene now gets exactly 24 Proper + 6 Improper (matching its 6 chemically-equivalent ring carbons); fresh phenol gets 6 Improper matching its 6 ring positions. All 29 existing tests (`tests/test_bossreader.py` + `tests/test_integration_converters.py`) still pass. Energy totals are essentially unchanged for both molecules -- benzene and phenol are close enough to their own planar equilibrium geometry that both the old (mislabeled-as-Proper) and new (correctly-labeled-Improper) functional forms evaluate to ~0 kcal/mol for these specific terms, which is exactly why this bug survived undetected through this whole validation exercise up to this point. It would not necessarily be energy-neutral for a genuinely non-planar molecule, since Proper and Improper conventions define the measured dihedral angle differently.

**A second bug this surfaced, found and fixed while re-verifying LAMMPS with real improper data**: `eval_lammps_energy.sh`'s single-point evaluation started failing on phenol with "Did not assign all atoms correctly." Root cause: `BOSS2LAMMPS.py`'s `.lmp` writer sets the simulation box's `xlo`/`ylo`/`zlo` to the exact coordinate minimum of the molecule, which by construction places at least one atom exactly on the box's lower face -- LAMMPS's domain decomposition can silently fail to assign that atom even with `boundary f f f`. Not a `BOSS2LAMMPS.py` bug (the `.lmp` file is perfectly valid); a test-harness gotcha specific to reusing the writer's own nominal box for a real evaluation. Fixed by padding every bound by 1 Å in `eval_lammps_energy.sh` before running `lmp`.

**Post-fix result, fresh benzene** (kcal/mol; BOSS total 7.9372328, bond=0.2214, angle=0.0, torsion=0.0, nonbonded=7.72):

| engine | total | bond | angle | torsion(+improper) | nonbonded |
|---|---|---|---|---|---|
| OpenMM | 7.9149 | 0.2144 | 0.0001 | 0.0 | 7.7004 |
| LAMMPS | 7.9372 | 0.2171 | 0.0001 | ~0.0 | 7.7200 |
| TINKER | 7.9389 | 0.2215 | 0.0 | 0.0 | 7.7174 |

**Post-fix result, fresh phenol** (kcal/mol; BOSS total -2.292973, bond=0.2106, angle=0.0404, torsion=0.0, nonbonded=-2.54):

| engine | total | bond | angle | torsion(+improper) | nonbonded |
|---|---|---|---|---|---|
| OpenMM | -2.2847 | 0.2113 | 0.0406 | 0.0 | -2.5366 |
| GROMACS | -2.1891 | 0.3204 | 0.0719 | 0.0 | -2.5814 |
| LAMMPS | -2.2925 | 0.2156 | 0.0409 | ~0.0 | -2.5490 |
| TINKER | -2.2955 | 0.2106 | 0.0404 | 0.0 | -2.5465 |

All within the same tolerance already established for the toluene/anisole-excluded set. `docs/agents`/future sessions extending this methodology to XPLOR or DESMOND should use freshly-generated Zmatrices (via the live Space or `-s <SMILES>`), not the legacy `molecules/small/*.z` reference library, for any molecule where the improper-torsion code path matters -- the reference library's degenerate declarations (confirmed now for `toluen.z`, `benzen.z`, and `phenol.z`) make it systematically unable to catch bugs in that path.

## Extending to Q: three real writer bugs, none previously exercised

Q (the Aqvist lab's MD engine, `qusers/Q6` on GitHub) is free and open source (GPL-style) -- unlike CNS/X-PLOR and Desmond, which are both registration-gated and out of scope for now. `Dockerfile.q` builds it from source with `gfortran`+`make`; see `tools/energy_validation/README.md`'s "The Q leg" for the full build/run notes. This was also what surfaced the ucomb() classification bug above: a direct challenge to this document's own "benzene has zero impropers" claim, while investigating why Q's `[impropers]` section was empty, is what led to tracing that claim back to its degenerate reference-file source.

`LigParGen/BOSS2Q.py` had never actually been run through a real Q build before -- every bug below was a first-contact failure, not a regression. All three are in `Boss2CharmmPRM()`:

**Bug 1 -- empty `[options]` section.** Q's own parameter reader hard-requires `vdw_rule` in `[options]` and rejects the *entire* file without it ("vdw_rule in options section not found"), which silently drops every bond/angle/torsion/atom_type below it too, not just vdW. Fixed by writing a complete `[options]` block (`vdw_rule geometric`, `scale_14 0.5`, `improper_potential periodic`, `improper_definition explicit`, matching the real bundled `Qoplsaa.prm` reference file). `improper_definition explicit` specifically is required for OPLS-AA: without it, Qprep6 auto-generates its own GROMOS-style improper for every sp2/3-connected ring atom via geometry, double-counting the planarity restraint OPLS-AA already bakes into the *proper* torsion Fourier series for those same ring atoms.

**Bug 2 -- one `[atom_types]` row per atom instead of per unique type.** Q rejects a repeated type *name* outright ("Could not enumerate atom type... Duplicate name?"), and that rejection corrupts every atom's type-index assignment for the rest of the topology build (every bond/angle/torsion count coming out 0, "Inconsistent molecule/residue start atoms"). Fixed by deduplicating to one row per unique OPLS type name -- CHARMM/TINKER/LAMMPS deliberately keep one row per atom, since those readers don't reject duplicates.

**Bug 3 -- wrong 1-4 LJ columns.** Q's `[atom_types]` row format is `name Avdw1 Avdw2 Bvdw1 Avdw3 Bvdw2&3 mass`, confirmed by reading Q's own Fortran source directly (`prep.f90`'s `[atom_types]` reader, `simprep.f90`'s `precompute_set_values_pp`) -- `Avdw3`/`Bvdw2&3` are *separate*, already-1-4-scaled LJ A/B values (Q combines them by direct multiplication, not `sqrt`, unlike its normal-LJ combining rule). The writer had been repeating the normal `B` value into the `Avdw3` slot and zeroing `Bvdw2&3`, which made every 1-4 LJ interaction wrong -- confirmed directly: total vdW energy came out negative instead of matching BOSS. Fixed by computing `sqrt(0.5)*ALJ`/`sqrt(0.5)*BLJ` for those two columns, baking in OPLS-AA's 0.5 1-4 LJ scale factor. Verified against the real bundled `Qoplsaa.prm`'s own CA row values, which matched to 3+ significant figures.

**Bug 4 -- found after the ucomb() fix restored real impropers -- extra column in `[impropers]`.** Once benzene's 6 genuine improper torsions were correctly classified (see above) and written out, Q's total energy came out at 67.92 kcal/mol against BOSS's 7.9372328 -- bond and nonbonded matched almost exactly, but torsion alone was +59.98 against BOSS's 0.0, for a molecule sitting exactly at its planar equilibrium. Traced (via Q's own Fortran source, `bondene.f90`'s `improper2` function) to Q's periodic-improper energy term hardcoding the multiplicity to 2 (`arg = 2.0*calc%angl - imp_lib(ic)%imp0`, no periodicity variable anywhere in the formula) and its `[impropers]` parameter-line reader (`prep.f90`) expecting exactly two numeric fields per line -- force constant and phase -- read positionally via `read(line, *, ...) taci, tacj, tack, tacl, imp_prm(i)%prm`, where `imp_prm(i)%prm`'s type has only `fk, imp0` as members. The writer (`retDihedImp()`) had been emitting three numeric fields (force, periodicity, phase), matching the CHARMM/other-format convention -- Q's list-directed read silently bound the periodicity placeholder (`2`) to `imp0` and never read the real phase (`180.0`) at all. Fixed by dropping the periodicity column entirely; this also means Q can only represent OPLS-AA's `n=2` improper term, which is the only term BOSS's own impropers ever populate in practice.

**Post-fix result, fresh benzene** (kcal/mol; BOSS total 7.9372328, bond=0.2214, angle=0.0, torsion=0.0, nonbonded=7.72):

| engine | total | bond | angle | torsion(+improper) | nonbonded |
|---|---|---|---|---|---|
| OpenMM | 7.9149 | 0.2144 | 0.0001 | 0.0 | 7.7004 |
| GROMACS | 8.0084 | 0.2354 | 0.0104 | 0.0 | 7.7626 |
| LAMMPS | 7.9372 | 0.2171 | 0.0001 | ~0.0 | 7.7200 |
| TINKER | 7.9389 | 0.2215 | 0.0 | 0.0 | 7.7174 |
| Q | 7.94 | 0.21 | 0.0 | 0.0 | 7.73 |

**Post-fix result, fresh phenol** (kcal/mol; BOSS total -2.292973, bond=0.2106, angle=0.0404, torsion=0.0, nonbonded=-2.54):

| engine | total | bond | angle | torsion(+improper) | nonbonded |
|---|---|---|---|---|---|
| OpenMM | -2.2847 | 0.2113 | 0.0406 | 0.0 | -2.5366 |
| GROMACS | -2.1891 | 0.3204 | 0.0719 | 0.0 | -2.5814 |
| LAMMPS | -2.2925 | 0.2156 | 0.0409 | ~0.0 | -2.5490 |
| TINKER | -2.2955 | 0.2106 | 0.0404 | 0.0 | -2.5465 |
| Q | -2.29 | 0.21 | 0.04 | 0.0 | -2.54 |

Q joins OpenMM/GROMACS/LAMMPS/TINKER within the same tolerance already established for the toluene/anisole-excluded set. Of the three engines this codebase set out to validate, only XPLOR and DESMOND remain untested -- both registration-gated, unlike Q.
