# LigParGen

OPLS-AA/CM1A(-LBCC) force-field parameter generator for organic ligands, from the Jorgensen group at Yale. LigParGen drives [BOSS](https://zarbi.chem.yale.edu/software.html) (Biochemical and Organic Simulation System) to derive OPLS-AA atom types and 1.14\*CM1A(-LBCC) partial charges for a molecule, then writes out ready-to-use parameter/topology files for a wide range of simulation packages.

This is the **3.0** release: a Python 3.11/3.12 port with a Docker-only distribution path. See [CHANGELOG.md](CHANGELOG.md) for the full version history, and the [original webserver](https://zarbi.chem.yale.edu/ligpargen/) this project also powers.

Supported output formats:

- OpenMM (`.xml` + `.pdb`)
- CHARMM/NAMD (`.prm` & `.rtf`)
- GROMACS (`.itp` & `.gro`)
- CNS/X-PLOR (`.param` & `.top`)
- Q (`.Q.prm` & `.Q.lib`)
- DESMOND (`.cms`)
- LAMMPS (`.lmp` & data files)
- TINKER (`.xyz` & `.key`)
- PDB2PQR (`.pqr`)
- BOSS Z-matrix (`.z`)

## Requirements

- Python 3.11 or 3.12
- Docker, plus your own licensed copy of [BOSS](https://zarbi.chem.yale.edu/software.html) — BOSS is proprietary and is never committed to this repo or baked into any published image (see [docs/adr/0001](docs/adr/0001-boss-binary-supplied-locally-never-published.md)). You'll need a working BOSS install on your own machine to build or run either of the two ways to use LigParGen below.

## Usage

There are two ways to run LigParGen; both need BOSS supplied locally at runtime, never committed to git or baked into a published image.

### 1. Command-line, via Docker

Build the image against your own BOSS install with `build.sh`, which stages a trimmed copy of BOSS into the (gitignored) build context and builds the Dockerfile:

```bash
./build.sh /path/to/your/boss/install
```

Then run it, mounting your working directory:

```bash
docker run --rm -v $(pwd):/work -w /work ligpargen:dev -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0 -l
```

Input options (see `LigParGen/Converter.py`'s `--help` for the full flag list):

| Flag | Input |
|---|---|
| `-s`/`--smiles` | SMILES string |
| `-p`/`--pdb` | PDB file (must include hydrogens) |
| `-m`/`--mol` | MDL MOL file (must include hydrogens) |
| `-z`/`--zmat` | BOSS Z-matrix — skips structure generation, uses the Z-matrix as-is |
| `-r`/`--resname` | 3-letter residue name |
| `-c`/`--charge` | net charge (`0`, `±1`, `±2`) |
| `-o`/`--opt` | optimization level (`0`–`3`) |
| `-l`/`--lbcc` | use 1.14\*CM1A-LBCC charges instead of 1.14\*CM1A (neutral molecules only) |

`-s` and `-p` can be combined: `-p file.pdb -s 'SMILES'` uses the SMILES as a trusted template to fix the PDB's bond orders and fill in any missing hydrogens via RDKit's `AssignBondOrdersFromTemplate` (PDB files carry no bond-order information and are often missing hydrogens) — see `LigParGen/mol_boss.py`'s `convert_pdb2mol_with_smiles`.

### 2. Web UI, via a Hugging Face Space

Live at **https://huggingface.co/spaces/lsdodda/ligpargen**. The [`space/`](space/) subdirectory holds a Gradio app mirroring the core of the [original LigParGen webserver](https://zarbi.chem.yale.edu/ligpargen/): SMILES input (typed or drawn with a Ketcher structure editor), PDB/MOL upload, optimization/charge-model options, results with all output formats, plus a 3D preview of the optimized geometry the original site doesn't have. Same BOSS-supply constraint as the CLI, adapted for Hugging Face's infrastructure: BOSS is fetched into the running container at startup from a private HF Dataset asset store via an HF Secret, never committed to git or baked into the Space's image (see [docs/adr/0003](docs/adr/0003-boss-hosted-via-runtime-fetch-not-build-time-copy.md)). The Space is public and includes a job timeout, a basic per-IP rate limit, and anonymous-only usage metrics (see [docs/adr/0004](docs/adr/0004-public-launch-reliability-and-anonymous-usage-metrics.md)). It also runs as an MCP server, so agents can call it as a tool instead of a human using the form. See [`space/README.md`](space/README.md) for details.

## Repo layout

```
LigParGen/          the Python package (CLI entry point: LigParGen.Converter:main)
space/               Hugging Face Space (Gradio web app)
tools/energy_validation/  BOSS vs. OpenMM/GROMACS/NAMD single-point energy validation (see its own README)
docs/adr/            architecture decision records
docs/agents/         issue-tracker/domain-doc config for AI coding agents
tests/fixtures/      regression fixtures (real BOSS-generated reference outputs)
CHANGELOG.md         version history
Dockerfile, build.sh CLI Docker image
```

## Citing LigParGen

If you use LigParGen, please cite:

1. Dodda, L. S.; Cabeza de Vaca, I.; Tirado-Rives, J.; Jorgensen, W. L. **LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands.** *Nucleic Acids Res.* **2017**, *45* (W1), W331-W336. [doi:10.1093/nar/gkx312](https://doi.org/10.1093/nar/gkx312)
2. Dodda, L. S.; Vilseck, J. Z.; Tirado-Rives, J.; Jorgensen, W. L. **1.14\*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations.** *J. Phys. Chem. B* **2017**, *121* (15), 3864-3870. [doi:10.1021/acs.jpcb.7b00272](https://doi.org/10.1021/acs.jpcb.7b00272)
3. Udier-Blagovic, M.; Morales De Tirado, P.; Pearlman, S. A.; Jorgensen, W. L. **Accuracy of free energies of hydration using CM1 and CM3 atomic charges.** *J. Comput. Chem.* **2004**, *25*, 1322-1332. [doi:10.1002/jcc.20059](https://doi.org/10.1002/jcc.20059)
4. Jorgensen, W. L.; Tirado-Rives, J. **Potential energy functions for atomic-level simulations of water and organic and biomolecular systems.** *Proc. Natl. Acad. Sci. USA* **2005**, *102*, 6665-6670. [doi:10.1073/pnas.0408037102](https://doi.org/10.1073/pnas.0408037102)

## Authors

* [Leela S. Dodda](https://github.com/leelasd) — `leela.dodda@yale.edu`
* Israel Cabeza de Vaca — `israel.cabezadevaca@yale.edu`
* Ayan Bhattacharjee — `abhattacharjee.me@gmail.com`
* [Matt Robinson](https://github.com/mc-robinson)
