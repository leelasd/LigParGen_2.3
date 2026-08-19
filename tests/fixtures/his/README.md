# Histidine (Ac-His-NHMe) BOSSReader parsing fixture

Part of #23 (parent map: #19). Captured `out`/`sum` text used by
`tests/test_bossreader.py` to exercise `BOSSReader`'s parsing methods
directly, without Docker/BOSS/a license -- the checked-in text below is all
these tests need.

Unlike `tests/fixtures/phenol/` (LigParGen's own auto-generated Zmat: 2
leading dummy atoms, "C00"/"H0A"-style generated atom names, no ring), this
fixture is one of BOSS's own hand-crafted reference Zmats
(`peptide/his.z` from a real BOSS install, never committed here -- see
"Provenance" below) and exercises everything phenol's fixture doesn't:

- **3 leading dummy atoms** (`DU1`/`DU2`/`DU3`), not LigParGen's usual 2 --
  this is exactly the case that broke `BOSS2OPENMM.py`'s old hardcoded
  `st_no = 3` (see `LigParGen/boss_common.py`'s `bossData()` docstring).
- **Non-"C00"-style atom names** -- real amino-acid names in the `out` file
  (`CB`, `CG`, `ND1`, `CE1`, `NE2`, `CD2`, ...), not LigParGen's own
  generated `C00`/`H0A` convention.
- **A ring (imidazole) with completed `ADD_DIHED`/Additional Bonds/Angles**
  -- see "Ring completion" below for why and how this was necessary.

## Molecule

`Ac-His-NHMe`, BOSS's own OPLS-AA/M reference Zmat for a capped histidine
residue (neutral, `peptide/his.z` in a licensed local BOSS install at
`/Users/leelasdodda/Codes/WLJ/boss/molecules/peptide/his.z` -- never
committed to this repo; see `docs/adr/0001`). It already carries explicit
OPLS atom types (columns 2/3 of each atom line, e.g. `505`, `507`, `511`)
rather than LigParGen's usual "auto-typed from a fresh SMILES" C00-style
Zmat, so it was run with `charge=0, lbcc=False` (no CM1A charge-fitting
pass) -- `optim=0` still runs BOSS's single-point `xSPM` pass
unconditionally (see `BOSSReader.Get_OPT`), which is what produces the
`out`/`sum` text captured here.

## Ring completion

Root cause (confirmed directly on this molecule): a hand-crafted Zmat's
ring bonds/angles are correct in each atom's own NA/NB/NC (bond/angle/
dihedral-partner) tree columns for whatever the tree's own spanning-tree
walk covers, but a ring's closing bond -- imidazole's `CG-CD2` here -- is
never implied by that tree alone, and *nothing else the tree implies about
the ring gets tabulated into BOSS's own Bond Stretching / Angle Bending /
Fourier Coefficients parameter tables either*, until the ring's bonds and
angles are explicitly re-declared in the Zmat's "Additional Bonds follow" /
"Additional Bond Angles follow" sections. This isn't just the one closing
bond: BOSS did not tabulate *any* of the 5-membered ring's own internal
bonds/angles (`CG-ND1`, `ND1-CE1`, `CE1-NE2`, `NE2-CD2`, `CD2-CG`, plus
every angle spanning them and their substituents) from the raw
`peptide/his.z` file as-is, even though the tree's own columns already
encode most of that connectivity.

This was diagnosed empirically, the same way the bug-fixing session earlier
in this repo's history diagnosed it (see `LigParGen/boss_common.py`'s
`bossData()`/`pair_declared_torsions()` docstrings for the general shape of
this class of bug):

1. Ran `BOSSReader(zmat='his.z', optim=0, charge=0, lbcc=False)` once,
   unmodified, from `/tmp` (matching `Converter.py`'s own
   `os.chdir('/tmp/')` convention -- `BOSSReader.Get_OPT` hardcodes
   `/tmp/sum`/`/tmp/out` as both BOSS's own cwd-relative output location
   *and* where `get_ImpDat()` reads them back from). Result:
   `MolData['BONDS']` had only 20 rows, `MolData['ANGLES']` only 36.
2. Built ground-truth connectivity from BOSS's own reported coordinates
   (`MolData['XYZ']`) via RDKit's `rdDetermineBonds.DetermineBonds`
   (geometry-based bond perception, charge=0) -- the same technique
   `LigParGen/Orca2CM5charges.py`'s `xyz_prep()` already uses elsewhere in
   this codebase.
3. Diffed that complete bond/angle set against what BOSS's first pass
   *actually* tabulated (not against the Zmat's own tree redundancy --
   several already-tree-covered ring bonds/angles were still missing from
   BOSS's tables). This found 9 missing bonds and 13 missing angles, all
   involving the imidazole ring atoms (raw indices 14-22:
   `CB CG ND1 CE1 NE2 HE CD2 HNE HD`).
4. Injected those 9 bonds / 13 angles into the unmodified Zmat's
   "Additional Bonds follow" / "Additional Bond Angles follow" sections
   (see `tests/fixtures/his.z` vs `tests/fixtures/his_ring_complete.z` --
   the same pair issue #24's integration test drives BOSS with; a plain
   diff shows exactly this injection and nothing else). Dihedrals were
   deliberately
   left alone for BOSS's own second `xSPM` pass to re-derive from the now-
   complete bond/angle set, rather than injected directly.
5. Ran `BOSSReader(zmat='his.z', optim=0, charge=0, lbcc=False)` a *second*
   time, from a clean `/tmp`, on the ring-completed Zmat. Result:
   `MolData['BONDS']` grew to 29 rows, `MolData['ANGLES']` to 49, and
   `MolData['TORSIONS_BY_DECL_IDX']` now has 59 entries including real ring
   torsions (e.g. declared-index 41: atom 15, type 390/390, coefficients
   `-1.282  1.645  -0.017  0.000`). This second pass's `/tmp/out` and
   `/tmp/sum` -- captured before anything deletes them -- are what's
   checked in here.

The capture script implementing steps 1-5 (`expand_zmat.py`'s
`missing_bonds_angles()`/`inject_additional_sections()`, originally written
during this session's `BOSS2OPENMM.py` bug-fixing pass and reused verbatim
here) is not itself part of this fixture -- it's a one-off generation tool,
kept in scratch, not checked into `LigParGen/`.

## Environment used to generate this fixture

Same `ligpargen:dev` Docker image already built and available locally for
this session (per `Dockerfile` in the repo root) -- **not** rebuilt or
modified for this capture:

```
docker run --rm --platform linux/amd64 --entrypoint python3 \
  -v <scratch capture dir>:/work \
  -v <worktree>/LigParGen:/app/LigParGen:ro \
  ligpargen:dev /work/capture_his.py
```

`/app/LigParGen` was bind-mounted read-only from *this worktree's* current
source (post all of this session's `BOSSReader.py`/`boss_common.py`
hardening/fixes), overriding the image's own baked-in copy, so the captured
text reflects a real run against the current, already-fixed parsing code --
consistent with how `tests/fixtures/phenol/README.md` documents its own
capture. `BOSSdir=/opt/boss` (the image's baked-in BOSS install, per
`docs/adr/0001` -- never committed anywhere).

`rdkit` was imported before anything importing `openbabel`/
`LigParGen.BOSSReader`, per this session's established segfault-avoidance
rule.

## What each file is

The Zmat inputs themselves (`peptide/his.z`, unmodified, and the same file
with the 9 injected "Additional Bonds" / 13 injected "Additional Bond
Angles" lines described above) live at the top level of `tests/fixtures/`
as `his.z` / `his_ring_complete.z` -- not duplicated under this directory --
since issue #24's integration test (`tests/test_integration_converters.py`)
feeds BOSS those same two files directly, and `his_ring_complete.z` is the
actual `BOSSReader` input that produced the `out`/`sum` capture below (a
plain diff between the two `.z` files is exactly the injection and nothing
else).

- `out` -- BOSS's raw Z-matrix-with-Cartesian-coordinates dump plus its own
  banner, captured from `/tmp/out` right after the second
  `BOSSReader(...)` pass completed (same file `BOSSReader.get_ImpDat()`
  reads via `Refine_file('/tmp/out')`, and the same thing
  `tests/fixtures/phenol/boss/out` captures).
- `sum` -- BOSS's Zmat-with-declared-Additional-sections dump plus the
  final "Non-Bonded Parameters for QM" charge table, captured from
  `/tmp/sum` the same way (`Refine_file('/tmp/sum')`).
- `plt.pdb` -- BOSS's own Cartesian-coordinate PDB dump from the same pass,
  kept alongside for completeness (not used by any test in this fixture
  set today).

## Known limitations

- Only the *second* pass's `out`/`sum` are captured -- the first,
  ring-incomplete pass's raw text was captured to a scratch-only
  `pass1_captured/` during generation but not checked in, since nothing in
  this ticket's test scope needs it.
- `charge=0, lbcc=False` -- this fixture doesn't exercise the LBCC
  charge-correction path (phenol's fixture already covers that) or a
  nonzero net charge.
- Net charge parsing (`get_charge`) on this fixture: `Reference Solute` is
  `0.00000`, same as phenol's fixture -- a future ticket may want a fixture
  with a nonzero net charge if `get_charge()`'s non-zero-value handling
  ever needs its own regression coverage.

## Regenerating this fixture

1. Have a working `ligpargen:dev` image built per this repo's own
   `Dockerfile` (needs a locally-supplied, licensed BOSS install -- see
   `docs/adr/0001`; never commit or push the built image).
2. Reimplement (or recover from scratch, if still present)
   `expand_zmat.py`'s `missing_bonds_angles()` / `inject_additional_sections()`
   -- see "Ring completion" above for what they need to do.
3. Write a driver script that: imports `rdkit` first; runs
   `BOSSReader('his.z', optim=0, charge=0, lbcc=False)` from `/tmp` on an
   unmodified copy of `peptide/his.z`; builds ground-truth bonds from
   `MolData['XYZ']` via `rdDetermineBonds.DetermineBonds`; diffs that
   against `MolData['BONDS']`/`MolData['ANGLES']`; injects the missing
   entries into the Zmat's "Additional Bonds follow" / "Additional Bond
   Angles follow" sections; runs `BOSSReader` a second time, from a clean
   `/tmp`, on the patched Zmat; copies `/tmp/out`, `/tmp/sum`, `/tmp/plt.pdb`
   out before anything else touches `/tmp`.
4. Run that script inside the container with `/app/LigParGen` bind-mounted
   read-only from the worktree whose source you want reflected in the
   capture, e.g.:
   ```
   docker run --rm --platform linux/amd64 --entrypoint python3 \
     -v <scratch dir with his.z + the driver script>:/work \
     -v <worktree>/LigParGen:/app/LigParGen:ro \
     ligpargen:dev /work/<driver script>.py
   ```
