# Energy validation: BOSS vs. OpenMM vs. GROMACS

Checks that LigParGen's generated parameter files actually reproduce
BOSS's own energy for a given molecule and geometry, across output
formats -- not just that the CLI runs to completion without an exception.
Two real, previously-undiscovered bugs (a missing-atom-mapping bug for
real BOSS reference Zmatrices, and a silently-omitted `[ dihedrals ]`
section in the GROMACS writer) were found this way, not by reading code.
See `docs/adr/0006-boss-vs-openmm-energy-validation.md` for the full
findings; this directory is the reusable *methodology* -- the scripts,
setup, and gotchas -- for running it again.

## Why single-point, and why decomposed by term

Two established lessons this methodology follows directly, both from the
person who wrote this tool's own past cross-code validation work:

- **Compare true single-point energies on both sides, not an optimized
  trajectory's endpoint against the other code's single point.** Comparing
  "NAMD after minimization" against "OpenMM single-point on the same
  starting file" produced a spurious ~0.2% electrostatic-constant-shaped
  discrepancy that took a multi-day GitHub thread to run down before
  realizing NAMD's minimizer was quietly doing something different from a
  genuine `run 0` single point
  ([ParmEd#907](https://github.com/ParmEd/ParmEd/issues/907)). Every
  comparison here is BOSS's own `xSPM` (zero accepted/rejected moves) on
  one side and a from-scratch `Context.getState(getEnergy=True)` /
  `gmx mdrun -nsteps 0` with no minimization on the other -- never an
  optimize-then-compare.
- **Decompose by force-group term (bond/angle/torsion/nonbonded), don't
  just compare totals.** A total-energy match can hide compensating
  errors in different terms
  ([openmm#1463](https://github.com/openmm/openmm/issues/1463)) -- e.g.
  ParmEd's own `energy_decomposition_system` exists specifically because a
  single aggregate number isn't enough to trust a conversion. Every script
  here reports `bond=`, `angle=`, `torsion=`, `nonbonded=` (GROMACS also
  splits nonbonded further into `LJ-14`/`Coulomb-14`/`LJ (SR)`/
  `Coulomb (SR)` if you need finer resolution -- see
  `eval_gromacs_energy.sh`).

## Setup

1. Build the main repo's own CLI image first (from the repo root, needs
   your own licensed BOSS install -- see `docs/adr/0001`):
   ```bash
   ./build.sh /path/to/your/boss/install
   ```
2. Build the two derived validation images (adds OpenMM via pip, and
   GROMACS via `apt-get` -- both freely available, no license needed
   beyond what `ligpargen:dev` already required):
   ```bash
   cd tools/energy_validation
   ./build.sh
   ```
   This produces `ligpargen-openmm:dev` and `ligpargen-gmx:dev`.

## Running a comparison

```bash
./compare.sh <3-letter-resname> <path-to-zmatrix> [charge]

# e.g.
./compare.sh PHN ~/Codes/WLJ/boss/molecules/small/phenol.z 0
```

Prints BOSS's, OpenMM's, and GROMACS's total and per-term energies
(kcal/mol) for that Zmatrix's geometry. The resname **must be exactly 3
characters** -- LigParGen's PDB writer uses a fixed 3-column residue-name
field (`%3s`); anything longer silently overflows into the coordinate
columns and either corrupts them or makes a strict reader like OpenMM's
`PDBFile` raise `ValueError: could not convert string to float`. Hit this
directly while building this tool (used `PHNT`, 4 characters, by mistake)
-- it's an easy trap, not a rare one.

### What it needs as input: a real Zmatrix (`.z` file)

Two sources work:

1. **BOSS's own reference library**, `$BOSSdir/molecules/{small,drugs,peptide}/*.z`
   -- real, pre-parameterized OPLS-AA Zmatrices BOSS ships with. Good for
   quick, varied test cases, but genuinely uneven in quality (see
   "Known BOSS reference-library gotchas" below) -- don't treat a mismatch
   against one of these as automatically meaning something is broken in
   LigParGen.

2. **The current production pipeline's own output**, pulled back out of a
   real submission -- this is how you tell whether a discrepancy is a
   property of one specific (possibly old/degraded) input file, or a
   property of the current code generally:
   ```bash
   # Submit by SMILES to the live Space, or run the local CLI's -s path --
   # either way, download the resulting output zip and pull the .z file
   # out of it. Example against the live Space:
   curl -s -X POST "https://lsdodda-ligpargen.hf.space/gradio_api/call/run_ligpargen" \
     -H "Content-Type: application/json" \
     -d '{"data": ["Cc1ccccc1", null, "0", "1.14*CM1A (neutral or charged)", "0"]}'
   # -> {"event_id": "..."}
   curl -s -N "https://lsdodda-ligpargen.hf.space/gradio_api/call/run_ligpargen/<event_id>"
   # -> event: complete, data: [{"url": ".../<RESID>.zip", ...}, ...]
   curl -s -o out.zip "<that url>"
   unzip out.zip <RESID>.z
   ./compare.sh <RESID or any 3-letter name> <RESID>.z 0
   ```
   This is exactly how the missing-ring-torsion issue in `toluen.z` (see
   below) was confirmed as a property of that one legacy file and not the
   current pipeline: same molecule, submitted fresh, came back fully and
   correctly parameterized.

### Testing more than one molecule

`compare.sh` runs exactly one. Loop it for a batch:
```bash
for f in benzen phenol toluen furan anilin anisol; do
    ./compare.sh "${f:0:3}" ~/Codes/WLJ/boss/molecules/small/"$f".z 0
done
```
(Pick distinct 3-letter names if any of the first-3-characters collide --
`toluen`→`tol` and `toluic`→`tol` would, for instance.)

## What each script does

- **`gen_and_boss_energy.py`** (runs in `ligpargen-openmm:dev`): the
  hardest-won piece. Calls `LigParGen.BOSSReader.BOSSReader(...)` directly
  instead of going through `Converter.convert()`, because `convert()`
  calls `BOSSReader.cleanup()` at the very end, which deletes `/tmp/out`
  -- the only place BOSS's own single-point energy (`NEW E`) and its
  per-term breakdown (`EBNDNE`/`EANGNE`/`EDIHNE`/`ENBNE`) are ever
  printed. Nothing else in the package prints or persists them. Also
  writes the OpenMM XML/PDB and GROMACS itp/gro for that same geometry
  (via `mainBOSS2OPM`/`mainBOSS2GMX`) so the other two scripts have
  something to evaluate.
- **`eval_openmm_energy.py`** (runs in `ligpargen-openmm:dev`): loads the
  XML+PDB into a real OpenMM `System`/`Context`, assigns each `Force` its
  own force group, evaluates `nsteps=0` (no minimization), reports total
  and per-group energy.
- **`eval_gromacs_energy.sh`** (runs in `ligpargen-gmx:dev`): wraps the
  bare `.itp` LigParGen produces into a minimal system `.top` (see
  "GROMACS gotchas" below), runs `gmx grompp` + `gmx mdrun -nsteps 0` +
  `gmx energy`, converts kJ/mol → kcal/mol, reports total and per-term
  energy.
- **`compare.sh`**: orchestrates all three for one molecule in one
  command, in a throwaway temp directory.
- **`vacuum_sp.mdp`**: the GROMACS run parameters for a genuine
  single-point vacuum-equivalent evaluation.

## Known gotchas (each one real, each one hit while building this)

**GROMACS, general:**
- The pre-2020 `cutoff-scheme = group` is removed entirely in modern
  GROMACS (2020+) -- must use `Verlet`.
- The Verlet scheme flatly does not support `pbc = no`
  ("With Verlet lists only full pbc or pbc=xy with walls is supported")
  -- there's no way to ask GROMACS for a truly non-periodic evaluation the
  way OpenMM's `NoCutoff` gives you. `eval_gromacs_energy.sh` works around
  this by enlarging the molecule's box to 5nm and using a 2nm cutoff, so
  periodic images never interact -- the closest available approximation,
  not a mathematically identical one.
- `[ defaults ]` (nbfunc/comb-rule/gen-pairs/fudgeLJ/fudgeQQ) must appear
  before `[ atomtypes ]` anywhere in the assembled topology, or `grompp`
  fails with `Invalid order for directive atomtypes`. LigParGen's own
  `.itp` is a self-contained per-molecule fragment that never includes
  this -- it's conventionally the *wrapper* topology's job, which is why
  `eval_gromacs_energy.sh` writes one.
- `comb-rule = 3` (geometric mean for both sigma and epsilon) and
  `fudgeLJ = fudgeQQ = 0.5` match OPLS-AA's convention and the
  `coulomb14scale`/`lj14scale` values LigParGen's own OpenMM XML writer
  already uses -- keep these two representations consistent if you change
  one.
- `gmx energy`'s term names differ from OpenMM's: GROMACS reports
  `Ryckaert-Bell.` (proper torsions, RB functional form) and
  `Per. Imp. Dih.` (impropers, periodic functional form) as two separate
  categories where OpenMM/BOSS report one combined `torsion` term -- sum
  them for an apples-to-apples comparison (`eval_gromacs_energy.sh`
  already does this).

**`.gro` coordinate precision**: the format stores exactly 3 decimal
places in *nanometers* (0.01 Å resolution) -- a full order of magnitude
coarser than the `.pdb` format's 3 decimals in *Ångströms* (0.001 Å) that
the OpenMM path reads. Bond terms have by far the largest force constants
of any term, so they're the most sensitive to this; a real, small,
GROMACS-specific bond-energy gap traced to this is documented in
`docs/adr/0006`. Not a `BOSS2GMX.py` bug -- a property of the `.gro`
format itself.

**Known BOSS reference-library gotchas** (`$BOSSdir/molecules/*/*.z`):
this is a corpus accumulated over BOSS's multi-decade development history,
not a uniformly-curated, single-pass-validated set -- expect real
unevenness, confirmed directly (not assumed) more than once this session:
- Dummy-atom placement isn't consistent: some files use 2 leading dummies
  (LigParGen's own auto-generated convention), others place a real atom
  first and interleave dummies after it (e.g. `acetam.z`, `ammonia.z`,
  `methan.z`, `ethane.z`). LigParGen now handles both (fixed in
  `boss_common.py`'s `bossData()`), but don't assume every reference file
  follows the convention your own test molecules happen to use.
- Some files (e.g. `methan.z`, `ethane.z`) make BOSS's own `xSPM` write an
  empty `plt.pdb` -- their Zmatrix header uses an older format BOSS's
  script apparently doesn't handle for a bare single-point run. BOSS's own
  energy (`NEW E`) is still printed correctly in `/tmp/out` even when this
  happens; only the coordinate output is affected, which blocks the
  OpenMM/GROMACS side of the comparison specifically, not the BOSS side.
- `toluen.z` has **zero declared ring-torsion parameters at all** in
  BOSS's own Fourier Coefficients table (confirmed directly, not
  inferred) -- an apparently under-parameterized/legacy file, not a
  general property of toluene. A freshly-generated toluene through the
  current pipeline gets the full, correct set (real `k2=15.167` ring
  torsions, `k2=10.46` planarity-enforcing impropers). See
  `docs/adr/0006`'s follow-up section for the full before/after
  comparison that established this.
- The "Additional Dihedrals" section banner text in BOSS's `/tmp/sum`
  output varies -- `"Additional Dihedrals follow (6I4)"` (from a full
  `xZCM1A` optimization run, LigParGen's own auto-generated Zmats always
  produce this) vs. `"Additional Dihedrals (6I4) - Zero types not shown"`
  (from a bare `xSPM` single-point run, which is what feeding a `-z`
  reference file through goes through). `BOSSReader.py`'s
  `find_boss_sections()` now recognizes both.

## What's not covered yet

**NAMD**: not run. Like BOSS itself, NAMD requires a separately licensed
binary from UIUC (not available via any package manager, unlike GROMACS)
-- someone would need to supply one, the same way BOSS is supplied
locally per `docs/adr/0001`. If you have a licensed NAMD install, the
methodology here should extend directly: NAMD reads CHARMM-format
files (`LigParGen`'s own `.rtf`/`.prm` writer, `BOSS2CHARMM.py`), and its
own `run 0` (not a minimization step -- see the ParmEd#907 discussion
above for why that distinction matters) gives a genuine single-point
energy comparable to the others.

**Other output formats** (`.key`/TINKER, `.Q.prm`/Q, `.top`+`.param`/
XPLOR, `.cms`/DESMOND, `.lmp`/LAMMPS): not energy-validated at all yet.
The `[ dihedrals ]`-omission bug this methodology found was specific to
`BOSS2GMX.py`'s own code (confirmed via `grep` that no other writer has
the identical broken condition) -- but that doesn't mean the others are
correct, only that they weren't specifically checked. Each would need its
own `eval_<format>_energy.*` script analogous to the two here, using
whatever engine reads that format natively.
