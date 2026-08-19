# Phenol regression fixture (SMILES -> BOSS -> OpenMM / GROMACS)

Resolves issue #9. Parent map: issue #4.

This is a **baseline captured from the current, unmodified `LigParGen/*.py`
source**, run against a real BOSS binary. It exists so a future Python-3
port can diff its own output against real, known-good output rather than
against a "we hope this is what it should look like" guess. No source under
`LigParGen/` was edited to produce this fixture.

## Molecule

Phenol, `c1ccc(cc1)O`, resname `PHN`, neutral (charge 0), CM1A-LBCC charges.
This is the README's own worked example. It was used as-is — small, a
single ring, one exocyclic OH (a rotatable bond), well-characterized
chemistry, and it already exercises the LBCC charge-correction path
(`-l` flag). No formal charge is exercised (see "Known limitations" below).

## Exact command

```
LigParGen -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0 -l
```

Run from `/out` (an empty working directory) inside the container described
below. Exit code 0. Full stdout/stderr are captured in `run_stdout.txt` /
`run_stderr.txt` in this directory (a handful of pandas `FutureWarning`s
about `DataFrame.drop`'s positional `axis` argument are expected and
harmless with pandas 1.5.3 — see "Environment" below).

## Environment used to generate this fixture

**Important: these are deliberately OLDER pins than the ones recommended
for the eventual Python-3 port (issue #6's findings, `numpy==2.4.6` /
`networkx==3.5` / `pandas==2.3.3` / `rdkit==2025.9.6`). This baseline must
run the CURRENT source unmodified, and the current source has three known
incompatibilities with those newer libraries (issue #7's `np.int` /
`networkx .node` findings, plus one more found here — see below). Do not
copy these pins into the ported code; do not compare this fixture's exact
warning noise against a run on the newer pins.**

- Base image: `python:3.8-slim-bookworm`, built for `linux/amd64` (per
  issue #5's Docker research), with `dpkg --add-architecture i386` +
  `libc6:i386` / `libstdc++6:i386` / `libgcc-s1:i386` + `csh`/`tcsh` so the
  real BOSS binary and its csh driver scripts run.
- `numpy==1.23.5`, `networkx==2.3` — the pins issue #7 identified as what
  the current code was actually written against (`np.int` alias still
  present; networkx's `.node[...]` graph accessor not yet removed).
- **Python 3.8, not 3.10/3.11 as originally suggested.** `networkx==2.3`
  does `from fractions import gcd` at import time (`networkx/algorithms/dag.py`),
  and `fractions.gcd` was removed in Python 3.9. This import fails
  immediately on Python 3.9/3.10/3.11 regardless of anything LigParGen
  does. Python 3.8 is the newest interpreter `networkx<2.4` actually
  imports on. This was verified empirically (see below) before proceeding,
  as the ticket asked.
- `pandas==1.5.3` — pinned **below 2.0**, not left unpinned. pandas 2.0
  made `DataFrame.drop()`'s `axis` argument keyword-only; the current code
  calls it positionally in several converters (e.g.
  `LigParGen/BOSS2OPENMM.py:83`: `final_df.drop([...], 1)`), which raises
  `TypeError: drop() takes from 1 to 2 positional arguments but 3 were
  given` on pandas >=2.0. This is a **third** real Python-3/library-version
  incompatibility beyond the two issue #7 found (that research scoped only
  Python-2-vs-3 syntax and the two specific numpy/networkx APIs, not the
  full pandas API surface) — noted here for whoever picks up the port.
  pandas 1.5.3 still accepts the positional call (with a `FutureWarning`,
  visible in `run_stderr.txt`).
- `rdkit==2024.3.5` (resolver-picked; no pin needed — just "compatible with
  Python 3.8 and numpy 1.23", per the ticket).
- `openbabel==3.1.1.1` Python bindings, pip-installed from source (needs
  `swig` + the apt `libopenbabel-dev` headers at
  `/usr/include/openbabel3`, and the actual header path had to be passed
  explicitly — the package's own build-time guess was wrong on this base
  image). This satisfies `LigParGen/mol_boss.py`'s `import openbabel`,
  which is a hard import-time dependency of `BOSSReader.py` even on the
  SMILES path where the bound functions themselves aren't called.
- apt `openbabel` (3.1.1, providing the `obabel` CLI) for the actual
  file-conversion shell-outs.
- **`babel` CLI compatibility shim.** Debian's `openbabel` package on
  bookworm ships only `obabel`, not the classic `babel` binary that
  `LigParGen/CreatZmat.py` shells out to (`os.system('babel -i... -o... ')`).
  Beyond the name, the CLI syntax changed: OpenBabel 2.x's `babel` accepted
  `-o<fmt> <outfile>` (format and output filename as one flag plus a bare
  positional argument); OpenBabel 3.x's `obabel` requires an explicit
  `-O <outfile>` — without it, `-o<fmt> <outfile>` is parsed as "format,
  write to stdout" with `<outfile>` treated as an *extra input file*,
  which doesn't exist yet and errors. A small shell wrapper installed as
  `babel` on `PATH` translates old-style invocations to `obabel`'s syntax.
  This is an environment compatibility shim, not a change to any
  `LigParGen/*.py` file — issue #8 (not applied here) is the separate,
  not-yet-done ticket to migrate off the `babel` CLI entirely in favor of
  the openbabel Python API.
- BOSS itself: the real, licensed binary, bind/copied in from
  `/Users/leelasdodda/Codes/WLJ/boss/` at build time (never committed to
  git, never pushed anywhere — see `docs/adr/0001`). `BOSSdir` set to
  point at it. `MCPROdir` deliberately left **unset**, per `docs/adr/0002`
  — though this only matters for PDB input's conformer-clustering
  fallback, not the SMILES path used here.
- The Dockerfile used to build this environment was a throwaway, local-only
  file (not committed — issue #10 is the separate ticket for the real
  production Dockerfile). It was deleted after use, per the ticket. The
  image and container were both removed (`docker rmi`) after the fixture
  was captured; nothing was pushed to any registry.

**Verification before proceeding** (per the ticket): before running the
actual conversion, `np.int`, `networkx`'s `G.node[...]` accessor, and a
full `import LigParGen` / `from LigParGen.Converter import convert` were
confirmed to work cleanly under this exact pin set inside the container.

## What each file is

### `boss/` — raw BOSS output

- `PHN.z` — the final BOSS Zmatrix (LigParGen's own internal
  representation; also embeds the OPLS-AA/CM1A-LBCC partial charges once
  BOSS has computed them). This is the artifact `LigParGen/BOSSReader.py`
  reads back in to drive every downstream converter.
- `sum`, `log`, `olog`, `out`, `plt.pdb`, `optzmat` — BOSS's own raw
  intermediate output from the final single-point calculation
  (`$BOSSdir/scripts/xSPM`), captured just before
  `BOSSReader.cleanup()` deletes them from `/tmp` (LigParGen always
  deletes these at the end of a run; they were captured via a debug-only
  `/bin/rm` shim in the throwaway container, purely for this fixture —
  not something the real Dockerfile needs). `out` contains BOSS's
  Z-matrix-with-Cartesian-coordinates dump and its own banner
  (`BOSS 4.9 Linux Version June 2015`); `sum` contains the final
  "Non-Bonded Parameters for QM (AM1 CM1Ax1.14) Atoms" table — the raw
  per-atom CM1A-LBCC partial charges and OPLS-AA Lennard-Jones parameters
  LigParGen parses to build every converter's output.

### `openmm/` — OpenMM converter output (`LigParGen/BOSS2OPENMM.py`)

Confirmed by reading `boss2opm()`, not guessed from filenames:

- `PHN.xml` — the OpenMM `ForceField` XML (atom types, residue template,
  harmonic bond/angle force, periodic torsion force, nonbonded
  charge/sigma/epsilon per type). This is *the* OpenMM force-field file.
- `PHN.pdb` — companion PDB (`pdb_prep()`) giving OpenMM a `Topology` /
  starting coordinates to pair with the XML.
- `PHN.pqr` — a PDB-with-charges-and-radii "PQR" file (`pqr_prep()`), also
  written by `boss2opm()`'s pipeline; not itself an OpenMM input format,
  kept alongside for completeness.

### `gromacs/` — GROMACS converter output (`LigParGen/BOSS2GMX.py`)

Confirmed by reading `mainBOSS2GMX()`/`boss2gmx()`:

- `PHN.gro` — GROMACS coordinate file (title line, atom count, one
  fixed-width line per atom, box vectors).
- `PHN.itp` — GROMACS include-topology file: `[ atomtypes ]`,
  `[ moleculetype ]`, `[ atoms ]`, `[ bonds ]`, `[ angles ]`,
  `[ dihedrals ]` (both proper and improper), `[ pairs ]`.

### Top level

- `run_stdout.txt` / `run_stderr.txt` — full captured output of the
  `LigParGen` invocation above (includes the pandas `FutureWarning`s noted
  under "Environment"). Named `.txt` rather than `.log` so they aren't
  swept up by this repo's `.gitignore` (`*.log`).

## Known limitations / things worth flagging

- **BOSS's own Cartesian output (`boss/plt.pdb`, and downstream `PHN.gro`)
  is planar** — every atom sits at the same Y coordinate. This is real,
  unmodified BOSS behavior (its internal-coordinate-to-Cartesian builder
  for a genuinely planar molecule like phenol; not a bug introduced by
  this fixture-generation environment) — phenol itself is planar, so this
  is plausible, but it does mean this fixture doesn't exercise any
  out-of-plane / improper-dihedral geometry. Worth keeping in mind if a
  future ticket wants a fixture that stresses 3D geometry handling more.
- This fixture only captures the OpenMM and GROMACS converter outputs
  (this ticket's scope). The same run also produces CHARMM/NAMD, X-PLOR,
  Q, LAMMPS, DESMOND, and TINKER output (`LigParGen/Converter.py`'s
  `convert()` always runs every converter) — none of that was kept here,
  since it's out of scope for issue #9. Regenerate with the command above
  if a future ticket needs those too.
- Phenol has no formal charge and no stereocenter, so this fixture doesn't
  exercise the anion/cation charge path (`-c` other than 0) or the
  MCPRO-absent PDB-input fallback (ADR-0002) — SMILES input never touches
  that code path regardless of `MCPROdir`. A future ticket may want a
  second fixture molecule to cover those.

## Regenerating this fixture

1. Build a throwaway Docker image per the "Environment" section above
   (base `python:3.8-slim-bookworm`, `--platform linux/amd64`, the i386 +
   csh/tcsh + openbabel + swig/libopenbabel-dev system packages, the
   `babel` CLI shim, `pip install numpy==1.23.5 networkx==2.3
   pandas==1.5.3 rdkit`, the source-built `openbabel==3.1.1.1` Python
   bindings, `pip install -e .` for this repo, `BOSSdir` pointed at a
   locally-supplied copy of BOSS, `MCPROdir` left unset).
2. Run `LigParGen -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0 -l` in an empty working
   directory inside the container.
3. Copy `PHN.z`, `PHN.xml`, `PHN.pdb`, `PHN.pqr`, `PHN.gro`, `PHN.itp` out
   of the run directory, and (if wanted) `sum`/`log`/`olog`/`out`/
   `plt.pdb`/`optzmat` out of `/tmp` *before* the run's own cleanup step
   deletes them.
4. Discard the image/container. Never push an image containing BOSS
   anywhere (per `docs/adr/0001`).
