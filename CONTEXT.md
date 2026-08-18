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
Grouping a set of candidate 3D geometries for a molecule and picking a representative one, done today via MCPRO's `clu` utility during PDB-input processing. Being replaced with an RDKit-based equivalent — MCPRO itself is being dropped as a dependency.

**MCPRO**:
BOSS's sibling program for free-energy perturbation / Monte Carlo work. Being fully dropped as a runtime dependency of LigParGen (its only current use is the `clu` conformer-clustering step). Note: the "MCPRO & BOSS Zmatrix" output format LigParGen can produce is just the shared `.z` file format — it does not invoke the MCPRO binary, and is unaffected by dropping MCPRO.
_Avoid_: assuming any mention of "MCPRO" in the README/output-format list implies the binary is still needed — check whether it means the file format or the program.
