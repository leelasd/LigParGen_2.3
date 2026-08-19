# Changelog

## Unreleased

- **PDB input can now be paired with a SMILES to fix bond orders.** PDB files
  carry no bond-order information and are often missing hydrogens; supplying
  `-s`/`--smiles` alongside `-p`/`--pdb` uses RDKit's
  `AssignBondOrdersFromTemplate` to recover the correct connectivity and bond
  orders from the SMILES (kept as the 3D coordinates from the PDB), then adds
  any missing hydrogens with RDKit-estimated positions. Mismatches between the
  SMILES and the PDB's heavy-atom connectivity raise a clear error instead of
  silently producing a wrong structure. See `LigParGen/mol_boss.py`'s
  `convert_pdb2mol_with_smiles`. The Hugging Face Space's PDB upload path now
  also accepts an optional accompanying SMILES for the same purpose.

## 3.0

Breaking release: drops Python 2 support, drops the MCPRO dependency, and
introduces a Docker-only distribution path for running LigParGen against BOSS.

- **Requires Python 3.11 or 3.12.** Fixed the handful of constructs that
  broke under modern interpreters/libraries: a removed `numpy` scalar alias
  (`np.int`), a removed `networkx` graph accessor (`Graph.node[...]`), and
  `pandas`'s `DataFrame.drop()` `axis` argument becoming keyword-only.
- **MCPRO dropped entirely.** Its only use (conformer clustering via
  `$MCPROdir/miscexec/clu` on the PDB input path) already had a working
  BOSS-only fallback in the existing code, so no replacement was needed —
  `$MCPROdir` is simply never set. See `docs/adr/0002`.
- **No more shelling out to the `babel` CLI.** `CreatZmat.py` and
  `BOSS2TINKER.py` now call the OpenBabel Python API directly
  (`OBConversion`/`OBOp`/`pybel`), matching the pattern already used in
  `mol_boss.py`. The Docker image needs no system Open Babel package.
- **New: a working Dockerfile.** Bundles a locally-supplied BOSS install
  (never committed to this repo or published in any image — BOSS is
  proprietary; see `docs/adr/0001`) with pinned dependencies
  (`numpy==2.4.6`, `pandas==2.3.3`, `rdkit==2025.9.6`, `networkx==3.5`,
  `openbabel==3.2.1`). See `build.sh` and the Dockerfile itself.
- **Out of scope for this release**: the ORCA/CM5 input path (`-q` flag),
  republishing to PyPI, and native (non-Docker) installation.

See issue [#4](https://github.com/leelasd/LigParGen_2.3/issues/4) for the
full planning history behind this release.
