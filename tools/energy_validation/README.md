# Energy validation: BOSS vs. OpenMM vs. GROMACS vs. NAMD vs. LAMMPS vs. TINKER

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
2. Build the five derived validation images (adds OpenMM via pip, GROMACS
   via `apt-get`, LAMMPS via `apt-get`, TINKER built from source via
   `git clone` + `cmake`, and Q built from source via `git clone` +
   `gfortran`/`make` -- all freely available, no license needed beyond
   what `ligpargen:dev` already required). On a non-amd64 host (e.g.
   Apple Silicon), `build.sh` passes `--platform linux/amd64` for you --
   `ligpargen:dev` itself is amd64-only (BOSS is a 32-bit x86 binary).
   The TINKER image takes a few minutes (compiling ~300 Fortran files
   from source -- see "The TINKER leg" below):
   ```bash
   cd tools/energy_validation
   ./build.sh
   ```
   This produces `ligpargen-openmm:dev`, `ligpargen-gmx:dev`,
   `ligpargen-lammps:dev`, `ligpargen-tinker:dev`, and `ligpargen-q:dev`.
3. (Optional) For the NAMD leg, get your own licensed NAMD install (not
   Dockerized -- see "The NAMD leg" below) and point `NAMD_DIR` at it:
   ```bash
   export NAMD_DIR=/path/to/your/namd/install   # containing namd3, psfgen
   ```
   Leave `NAMD_DIR` unset to skip NAMD and just get the
   BOSS/OpenMM/GROMACS/LAMMPS/TINKER/Q six-way comparison.

## Running a comparison

```bash
./compare.sh <3-letter-resname> <path-to-zmatrix> [charge]

# e.g.
./compare.sh PHN ~/Codes/WLJ/boss/molecules/small/phenol.z 0
# with NAMD_DIR set, this also prints a NAMD_ENERGY/NAMD_TERMS line
```

Prints BOSS's, OpenMM's, GROMACS's, LAMMPS's, TINKER's, and (if `NAMD_DIR`
is set) NAMD's total and per-term energies (kcal/mol) for that Zmatrix's
geometry. The resname **must be exactly 3
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
  writes the OpenMM XML/PDB, GROMACS itp/gro, CHARMM rtf/prm, LAMMPS lmp,
  and TINKER new.xyz/key for that same geometry (via `mainBOSS2OPM`/
  `mainBOSS2GMX`/`mainBOSS2CHARMM`/`mainBOSS2LAMMPS`/`mainBOSS2TINKER`) so
  the other scripts have something to evaluate.
- **`eval_openmm_energy.py`** (runs in `ligpargen-openmm:dev`): loads the
  XML+PDB into a real OpenMM `System`/`Context`, assigns each `Force` its
  own force group, evaluates `nsteps=0` (no minimization), reports total
  and per-group energy.
- **`eval_gromacs_energy.sh`** (runs in `ligpargen-gmx:dev`): wraps the
  bare `.itp` LigParGen produces into a minimal system `.top` (see
  "GROMACS gotchas" below), runs `gmx grompp` + `gmx mdrun -nsteps 0` +
  `gmx energy`, converts kJ/mol → kcal/mol, reports total and per-term
  energy.
- **`eval_lammps_energy.sh`** (runs in `ligpargen-lammps:dev`): runs `lmp`
  directly against the `.lmp` data file with free (non-periodic)
  boundaries (templated from `lammps_sp_template.in`), a genuine `run 0`,
  parses the `thermo_style custom` output line, reports total and
  per-term energy.
- **`eval_tinker_energy.sh`** (runs in `ligpargen-tinker:dev`): copies the
  `.key` file under the `.new.xyz`'s own basename (TINKER auto-associates
  a keyfile by matching filename), runs TINKER's own `analyze` program
  (a pure single-point evaluator -- no simulation, no minimization step
  even exists to accidentally use), parses its `Energy Component
  Breakdown` block, reports total and per-term energy.
- **`eval_namd_energy.sh`** (runs natively, NOT in Docker -- see "The NAMD
  leg" below): builds a `.psf` from the `.rtf`/`.prm`/`.pdb` via `psfgen`
  (templated from `psfgen_template.pgn`), then runs NAMD itself with a
  genuine `run 0` config (templated from `namd_sp_template.conf`), parses
  the `ENERGY:` line, reports total and per-term energy. Only runs when
  `NAMD_DIR` is set.
- **`compare.sh`**: orchestrates all of the above for one molecule in one
  command, in a throwaway temp directory.
- **`vacuum_sp.mdp`**: the GROMACS run parameters for a genuine
  single-point vacuum-equivalent evaluation.
- **`lammps_sp_template.in`**: the LAMMPS input-script template
  `eval_lammps_energy.sh` fills in with the resid.
- **`psfgen_template.pgn`**, **`namd_sp_template.conf`**: templates
  `eval_namd_energy.sh` fills in with the resid (`__RESID__` placeholder)
  for the psfgen and NAMD steps respectively.

## The LAMMPS leg

LAMMPS is Dockerized the same way GROMACS is (`Dockerfile.lammps`, adds
Debian's `apt-get install lammps` on top of `ligpargen-openmm:dev` --
freely available, GPL, no license needed). The installed binary is named
`lmp` (not `lmp_serial`/`lmp_mpi`, if you're used to older LAMMPS
packaging).

`BOSS2LAMMPS.py`'s `.lmp` writer is structurally different from the other
four writers: instead of deduplicating by parameter class (the way the
OpenMM/GROMACS/CHARMM writers all do, giving repeated `<Proper class1=...>`-
style entries a single shared type), every individual atom/bond/angle/
dihedral instance gets its own unique numbered type, with BOSS's raw
per-instance parameters written directly. This makes it the least-
transformed of any writer here, and it shows: LAMMPS's single-point energy
matched BOSS's own printed total to within 0.00005 kcal/mol on both
molecules tested -- see `docs/adr/0006`.

LAMMPS's built-in `dihedral_style opls` and `improper_style cvff` match
the `.lmp` writer's coefficient layout exactly, with **no unit conversion
or halving needed** -- unlike the OpenMM/CHARMM writers, which both
convert kcal→kJ (OpenMM only) and halve BOSS's raw Fourier coefficients
(`K = V_opls / 2`) to match the `1 + cos(...)` energy form those tools
expect. LAMMPS's `opls` dihedral style already bakes the equivalent 0.5
factor into its own formula, so it wants BOSS's raw, unhalved,
kcal/mol-unit V-coefficients directly -- confirmed directly: BOSS's own
raw `V2` for benzene's ring torsion is `7.250` kcal/mol, which is exactly
what `.lmp` writes, and is consistent with OpenMM's stored
`k2 = 15.167 kJ/mol = 7.250 * 4.184 / 2`.

`lammps_sp_template.in` uses free (non-periodic) boundaries
(`boundary f f f`) with a 100 Å cutoff for a true vacuum evaluation --
unlike modern GROMACS (see the GROMACS gotchas above), LAMMPS actually
supports genuinely non-periodic boundaries directly, so no enlarged-box
workaround is needed here. `special_bonds lj/coul 0.0 0.0 0.5` matches
OPLS-AA's 1-4 scaling convention, and LAMMPS's default pair-mixing rule
for `lj/cut` is already geometric-mean, matching OPLS-AA's own combining
rule -- so no `pair_modify mix` override is needed either, since every
atom gets its own unique type (per the writer's own per-instance
convention above) and cross-type LJ parameters are always obtained
through mixing, never given explicitly.

## The TINKER leg

TINKER is freely available (BSD-style academic license) but ships in
neither Debian's apt repo nor Homebrew, so `Dockerfile.tinker` builds it
from source: clones `TinkerTools/tinker` from GitHub and compiles just the
`analyze` target (TINKER's own pure single-point energy-evaluation
program -- no simulation, no minimization even exists in it to
accidentally use, so there's no "compare a trajectory endpoint" trap here
the way there was for NAMD). Two real build gotchas, both patched around
in the Dockerfile rather than waited on:

- TinkerTools' own `cmake/CMakeLists.txt` (an alternative to their
  officially-supported classic Makefile build) has an incomplete source
  file list -- `uatom.f` defines a Fortran module several other files
  `use`, but isn't in the curated `_FILES` list, so the build fails with
  `Cannot open module file uatom.mod` on the first file that needs it.
  Patched with a `sed` that adds it back in before configuring.
- The Fortran module (`.mod`) dependency graph isn't fully expressed as
  CMake target dependencies, so a parallel build (`-jN`, N>1) can try to
  compile a file before the module it needs exists yet, and fails the
  same way. Built with `-j1` (serial) instead -- slower (a few minutes)
  but doesn't hit the race.

**A real bug found here, not just a residual**: `BOSS2TINKER.py`'s
`.xyz`-file atom-type numbering didn't match its own `.key` file's
atom-type declarations at all (a placeholder `799 + atom_index` value
where it needed the real OPLS type number) -- this made every TINKER
output LigParGen has ever generated completely non-functional, confirmed
directly by running the freshly-built `analyze` against the (then-broken)
output and getting "Undefined Atom Type" for every atom. Fixed; see
`docs/adr/0006` for the full story and the post-fix numbers, which now
match BOSS as tightly as LAMMPS does.

`eval_tinker_energy.sh` copies the `.key` file under the `.new.xyz`'s own
basename before running `analyze` -- TINKER auto-associates a keyfile
with a coordinate file purely by matching filename (`<base>.key` next to
`<base>.xyz`), it isn't a `-k`-style command-line flag the way some other
TINKER-family tools use. `analyze`'s interactive "Enter Parameter File
Name" prompt is answered with a blank line (`echo ''`) -- correct, not a
workaround, since the `.key` file already carries every parameter
LigParGen generated with no external `oplsaa.prm` reference to point at.

## The NAMD leg

Unlike OpenMM/GROMACS, NAMD is **not Dockerized** -- like BOSS itself
(`docs/adr/0001`), it's a licensed, proprietary binary (from UIUC, not
available via any package manager) that must be supplied locally at
runtime, never committed to this repo or baked into any image. Point
`NAMD_DIR` at your own install directory (containing `namd3` and
`psfgen`) to include this leg; every script here treats its absence as
"skip NAMD," not an error.

NAMD reads CHARMM-format topology/parameters -- `LigParGen`'s own
`BOSS2CHARMM.py` writer produces the `.rtf`/`.prm` pair, and `psfgen`
(NAMD's own topology tool, ships alongside `namd3`) turns the `.rtf` plus
LigParGen's `.pdb` into a `.psf`. The NAMD config
(`namd_sp_template.conf`) sets `vdwGeometricSigma yes` -- OPLS-AA uses a
geometric-mean combining rule for sigma (as well as epsilon, which NAMD
always combines geometrically), and without this flag NAMD silently falls
back to CHARMM's arithmetic-mean sigma rule, giving the wrong LJ energy
for an OPLS-AA force field. This matches the setting in a real, historical
2017 NAMD-vs-OpenMM comparison config from this codebase's own author.
`run 0` (not `minimize`) is what makes this a genuine single point --
see "Why single-point" above for why that distinction is load-bearing, not
cosmetic (a `minimize 1000` NAMD config was in fact what that 2017
comparison used, and had to be replaced with `run 0` to get a valid
comparison here).

**macOS Gatekeeper**: an unsigned, downloaded NAMD binary (`namd3`,
`psfgen`) will be quarantined the first time you run each one --
"Apple could not verify... is free of malware." Run the binary once (it
will hang or get killed), then go to System Settings > Privacy & Security
and click "Allow Anyway" next to the message about that binary. This is a
security-relevant setting change, so do it yourself rather than having an
agent run `xattr -d com.apple.quarantine` for you -- needs doing once per
binary, per machine.

**Validated**: benzene and phenol (both already in the "8 of 10 match
closely" set from `docs/adr/0006`) were run through this NAMD leg and
matched BOSS/OpenMM/GROMACS to within the same tolerance as the other two
engines -- see `docs/adr/0006`'s NAMD section for the numbers.

## The Q leg

Q (the Aqvist lab's MD engine, `qusers/Q6` on GitHub) is free and open
source (GPL-style), unlike CNS/X-PLOR and Desmond -- `Dockerfile.q` builds
it from source with `gfortran`+`make`. The topology/energy pipeline is
two Q programs, not one: `Qprep6` (`q_prep_template.inp`) reads
LigParGen's `.lib`+`.Q.prm` and a PDB, builds a topology (`mt`), and
writes it out (`wt`); `Qdyn6` (`q_sp_template.inp`) reads that topology
and actually evaluates the energy. `eval_q_energy.sh` drives both and
parses `Qdyn6`'s "Energy summary at step 0" block, printed **before** its
one integration step runs -- a genuine single-point evaluation of the
input geometry.

Several real gotchas, all either worked around in the eval script or
fixed in `LigParGen/BOSS2Q.py` itself (not pre-existing writer bugs --
`BOSS2Q.py` had simply never been exercised against a real Q run before):

- **Qprep6's PDB reader** doesn't understand `TER`/`CONECT`/`END`/`REMARK`
  lines -- their presence makes it miscount "0 molecules" instead of 1,
  corrupting topology assembly. `eval_q_energy.sh` strips the PDB to
  `^ATOM` lines only before handing it to `Qprep6`.
- **Qdyn6 has no true `steps=0` single-point mode** ("Need at least one
  step of dynamics") and refuses `temperature=0` ("No dynamics at zero
  temperature!"). Worked around with `steps=1`, `temperature=0.001` --
  the pre-integration "Energy summary at step 0" block is unaffected by
  the one tiny step that follows it.
- **`BOSS2Q.py`'s `[options]` section was empty.** Q's own parameter
  reader hard-requires `vdw_rule` in `[options]` and rejects the *whole*
  file without it, silently dropping every bond/angle/torsion/atom_type
  below it too. Fixed by writing a complete `[options]` block, including
  `improper_definition explicit` (below).
- **`BOSS2Q.py`'s `[atom_types]` section wrote one row per atom.** Q
  rejects a repeated type *name* outright ("Could not enumerate atom
  type... Duplicate name?"), and that rejection corrupts every atom's
  type-index assignment for the rest of the topology build. Fixed by
  deduplicating to one row per unique OPLS type name (CHARMM/TINKER/
  LAMMPS give every atom its own row deliberately -- fine there, since
  those readers don't reject duplicates).
- **`BOSS2Q.py`'s 1-4 LJ columns were wrong.** Q's `[atom_types]` row
  format is `name Avdw1 Avdw2 Bvdw1 Avdw3 Bvdw2&3 mass`, where
  `Avdw3`/`Bvdw2&3` are *separate*, already-1-4-scaled LJ A/B values (Q's
  own `precompute_set_values_pp` combines them by direct multiplication,
  not `sqrt`). The writer had been repeating the normal `B` value into
  the `Avdw3` slot and zeroing `Bvdw2&3`, which made every 1-4 LJ
  interaction wrong (total vdW energy came out negative instead of
  matching BOSS). Fixed by computing `sqrt(0.5)*ALJ`/`sqrt(0.5)*BLJ` for
  those two columns, baking in OPLS-AA's 0.5 1-4 LJ scale factor -- the
  same value the real bundled `Qoplsaa.prm` reference file uses.
- **`improper_definition explicit` is required for OPLS-AA.** Without
  it, Qprep6 auto-generates its own GROMOS-style improper for every
  sp2/3-connected ring atom via geometry, double-counting the planarity
  restraint OPLS-AA already bakes into the *proper* torsion Fourier
  series for those same ring atoms.
- **`BOSS2Q.py`'s `[impropers]` lines had an extra column.** Q's own
  periodic-improper energy term hardcodes the multiplicity to 2 (its
  source computes `arg = 2*phi - imp0`, with no periodicity term at
  all), and its parameter-file reader expects exactly two numeric fields
  per line -- force constant and phase, read positionally. The writer
  was emitting three (force, periodicity, phase), so Q's read silently
  bound the periodicity placeholder to the phase and never read the real
  phase value -- confirmed directly: ~60 kcal/mol of spurious improper
  energy for benzene's perfectly planar ring, where BOSS's own energy is
  exactly 0. Fixed by dropping the periodicity column. This also means Q
  can only represent OPLS-AA's `n=2` improper term, which is the only
  term BOSS's own impropers ever populate in practice.

**Validated**: benzene and phenol (fresh, non-degenerate Zmatrices --
see `docs/adr/0006`'s note on the legacy-reference-file gap) matched
BOSS/OpenMM/GROMACS/LAMMPS/TINKER to the same tolerance as the other
engines once all of the above were fixed. See `docs/adr/0006`'s Q
section for the numbers.

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

**Other output formats** (`.top`+`.param`/XPLOR, `.cms`/DESMOND): not
energy-validated at all yet, and unlike Q, both are registration-gated
rather than freely downloadable, so validating them needs that access
sorted out first. Bugs this methodology has found so far (GROMACS's
silently-omitted `[ dihedrals ]`, TINKER's mismatched atom-type
numbering, and Q's `[options]`/`[atom_types]`/1-4-LJ/`[impropers]`
issues) were each specific to their own writer's code (confirmed via
`grep`/direct inspection that no other writer has the identical broken
pattern) -- but that doesn't mean XPLOR/DESMOND are correct, only that
they weren't specifically checked yet. Each would need its own
`eval_<format>_energy.*` script analogous to the ones here, using
whatever engine reads that format natively -- and, per the pattern
established by NAMD/TINKER/LAMMPS/GROMACS/Q above, would need that
engine's own licensing/availability checked before assuming it's as easy
to add as the free-and-open ones were.
