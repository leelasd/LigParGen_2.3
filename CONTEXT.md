# LigParGen

Python tooling that drives BOSS to derive OPLS-AA force-field parameters for a molecule, then converts BOSS's output into the input files a chosen MD/QM engine expects.

## Language

**BOSS**:
W.L. Jorgensen's Monte Carlo / quantum-mechanics simulation program (Biochemical and Organic Simulation System). It computes OPLS-AA parameters and partial charges for a molecule. It is a proprietary, closed-source Linux binary that LigParGen shells out to — not something this project builds or modifies.
_Avoid_: the engine, the simulator, the backend

**BOSSdir**:
The environment variable pointing at a local BOSS installation. LigParGen's shell-outs assert it is set and use it to locate both the BOSS binary and the BOSS driver scripts.

**BOSS driver script**:
One of the csh scripts under `$BOSSdir/scripts` (e.g. `xZCM1A`, `xOPT`, `xSPM`) that LigParGen invokes to run a specific BOSS calculation — charge computation, geometry optimization, or single-point energy. These ship with BOSS and are not edited by this project.
_Avoid_: wrapper script, BOSS command

**Zmatrix**:
BOSS's internal representation of a molecule — atom connectivity plus internal coordinates — read from and written to `.z` files. LigParGen builds the input Zmatrix, hands it to BOSS, and reads the result back out.
_Avoid_: z-file, geometry file

**OPLS-AA**:
The all-atom force field (bonded and non-bonded parameters, indexed by atom type) that BOSS parameterizes a molecule against. LigParGen's output converters translate OPLS-AA parameters into each target engine's own format.
_Avoid_: "the force field" alone — OPLS-UA (united-atom) also exists in BOSS's scope even though this project only targets OPLS-AA

**Partial-charge scheme**:
The method BOSS uses to assign partial atomic charges — CM1A, CM1A-LBCC (a bond-charge-corrected variant of CM1A), or CM5. Each scheme is a distinct, non-interchangeable calculation; LigParGen must know which one was requested to interpret BOSS's output correctly.
_Avoid_: charge model, charges (too generic — always name the specific scheme)

**Resname**:
The short identifier (e.g. `PHN`) that ties together a molecule's Zmatrix, its BOSS working files, and its converted output files. Supplied by the caller via the `-r` flag; analogous to a PDB residue name but scoped to a single LigParGen run, not a multi-residue structure.

**Conformer clustering**:
Grouping a set of candidate 3D geometries for a molecule and picking a representative one. On the PDB input path, LigParGen already falls back to a BOSS-only single-geometry substitute when MCPRO's `clu` utility isn't available — see ADR-0002. No RDKit (or other) replacement is being built for this; that was an earlier plan superseded once the existing fallback was found to already work.

**MCPRO**:
BOSS's sibling program for free-energy perturbation / Monte Carlo work. Being fully dropped as a runtime dependency of LigParGen (its only current use is the `clu` conformer-clustering step). Note: the "MCPRO & BOSS Zmatrix" output format LigParGen can produce is just the shared `.z` file format — it does not invoke the MCPRO binary, and is unaffected by dropping MCPRO.
_Avoid_: assuming any mention of "MCPRO" in the README/output-format list implies the binary is still needed — check whether it means the file format or the program.

### Hugging Face Space (web front end)

**Space**:
A Hugging Face Space — the hosted app unit (git repo, build, and running app together). Distinct from `space/`, the subdirectory in this repo holding the Space's source.
_Avoid_: "the app" alone — ambiguous with the LigParGen CLI/package itself

**Space visibility**:
One of three Hugging Face-defined access levels for a Space: **Public** (source, running app, and built image all fully open), **Protected** (source and image private to owner/collaborators, but the running app is still publicly reachable), **Private** (source, running app, and built image all restricted to owner/collaborators — 404s for anyone else, not listed in search). This project's Space starts at Private.

**BOSS asset store**:
The private Hugging Face Dataset repository holding a copy of the licensed BOSS install, fetched into the Space's container at startup using a Secret token. Never committed to this git repo, never baked into any Docker image layer — the hosted-deployment counterpart to ADR-0001's local-build-time rule.
_Avoid_: "the BOSS repo" — ambiguous with this GitHub repo

**Secret** / **Variable** (Space configuration):
Hugging Face Space configuration values. A **Secret** is a private credential — write-only once set, not copied to duplicated Spaces — used here to authenticate the BOSS asset store fetch. A **Variable** is a public, visible config value. Anything credential-shaped must be a Secret, never a Variable.

**Draw-a-molecule input**:
The alternate input path where a user sketches a 2D structure in an embedded editor instead of typing a SMILES string or uploading a file. Backed by the Ketcher editor (the real Yale site uses JSME instead — a deliberate deviation, not an oversight, since Ketcher already has a proven Gradio integration pattern to build on).

**3D structure preview**:
A rendered 3D view of the BOSS-optimized output geometry, shown before the result files are downloaded. Not present on the original Yale site — a genuine enhancement for this Space, backed by the `gradio_molecule3d` component.
