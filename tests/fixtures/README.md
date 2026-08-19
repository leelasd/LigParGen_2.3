# Corpus integration-test fixtures (issue #24)

Real BOSS reference Zmats used by `tests/test_integration_converters.py`,
copied from BOSS's own `molecules/peptide/` corpus
(`/Users/leelasdodda/Codes/WLJ/boss/molecules/peptide/` on the machine these
were captured on). Per `docs/adr/0001-boss-binary-supplied-locally-never-published.md`,
committing these Zmats as test data is fine -- only the BOSS binary itself is
restricted.

## Files

- **`his.z`** -- unmodified copy of `peptide/his.z` (`Ac-His-NHMe`). A real
  BOSS reference Zmat for the imidazole-ring histidine side chain, with 3
  leading dummy atoms (`DU1`/`DU2`/`DU3`, raw indices 1-3, first real atom
  at raw index 4) -- the `st_no` bug case `LigParGen/boss_common.py`'s
  `bossData()` now reads from the first real atom's own raw index rather
  than assuming LigParGen's own always-2-dummy convention.

  Its `Additional Bonds follow` section is empty and `Additional Bond
  Angles follow` only partially covers the ring: none of the imidazole ring
  atoms' (`CG`/`ND1`/`CE1`/`NE2`/`CD2`) own parent bonds are declared
  `Variable Bonds`, so BOSS's own `out`/`sum` never tabulates OPLS
  parameters for them even though each ring atom's own NA/NB tree columns
  already encode the connectivity -- confirmed directly: a first
  `BOSSReader(optim=0)` pass on this file tabulates only 20 bonds / 36
  angles.

- **`his_ring_complete.z`** -- the same molecule with the missing ring
  bonds/angles injected into the `Additional Bonds follow` / `Additional
  Bond Angles follow` sections (`diff his.z his_ring_complete.z` shows only
  those two sections changed -- 9 added bond lines, 13 added angle lines).
  A second `BOSSReader(optim=0)` pass on this file tabulates 29 bonds / 49
  angles -- the ring's own bonded terms now included, e.g. real,
  chemically-sensible `CG-ND1` bond/angle parameters where the raw `his.z`
  had none.

  How it was built (RDKit-geometry-based, run once against a real BOSS
  binary, not by hand):
  1. Run `BOSSReader('his.z', optim=0, charge=0, lbcc=False)` once to get
     real Cartesian geometry (`MolData['XYZ']`) and BOSS's own tabulated
     `MolData['BONDS']`/`MolData['ANGLES']`.
  2. Build an RDKit molecule from that geometry (element symbols from
     `LigParGen.boss_common.bossData()`'s periodic-table lookup, positions
     from `MolData['XYZ']`) and call
     `rdkit.Chem.rdDetermineBonds.DetermineConnectivity()` -- pure
     distance-based bond perception, no bond-order/charge guessing needed.
  3. Diff RDKit's complete bond/angle graph against what
     `MolData['BONDS']`/`MolData['ANGLES']` already had (offset by the raw
     Zmat's own first-real-atom index, i.e. `st_no`) to get the missing set.
  4. Inject the missing bonds/angles into the *original* (pre-run) Zmat
     text's `Additional Bonds follow`/`Additional Bond Angles follow`
     sections, in the exact fixed-width format `LigParGen/CreatZmat.py`'s
     `print_ZMAT()` writes (`'%4d%4d\n'` / `'%4d%4d%4d\n'`).
  5. Re-run `BOSSReader` on the result to confirm BOSS now tabulates a
     strictly larger BONDS/ANGLES set (20->29 / 36->49) -- i.e. the fixture
     actually exercises the bonded-term-completion path, not just a
     syntactically-plausible-looking edit.

  `tests/test_integration_converters.py::test_ring_completion_adds_tabulated_bonded_terms`
  re-asserts step 5 (bond/angle counts strictly increase) every run, so this
  stays true even if a future BOSS/RDKit version changes the exact numbers.

- **`ala.z`** -- unmodified copy of `peptide/ala.z` (`Ac-Ala-NHMe`). Plain
  baseline peptide fragment, no ring, no missing bonded terms (its own
  `Additional Bonds`/`Additional Bond Angles` sections already cover
  everything BOSS needs). Also happens to use the same 3-leading-dummy
  convention as `his.z` -- verified directly rather than assumed (BOSS's
  `peptide/` corpus Zmats are consistently 3-dummy; LigParGen's own
  auto-generated Zmats are the ones that are consistently 2-dummy).

## Regenerating `his_ring_complete.z`

Run inside a container with a real `$BOSSdir` and this repo's current
`LigParGen` installed (e.g. `ligpargen:dev`, built via `./build.sh`):

```python
from LigParGen.BOSSReader import BOSSReader
from LigParGen.boss_common import bossData
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
import networkx as nx

original_text = open("his.z").read()
mol = BOSSReader("his.z", 0, 0, False)
types, Qs, num2opls, st_no, num2typ2symb, num2pqrtype = bossData(mol)

xyz = mol.MolData["XYZ"]
rw = Chem.RWMol()
conf = Chem.Conformer(len(types))
for i in range(len(types)):
    rw.AddAtom(Chem.Atom(num2typ2symb[i][3]))
    conf.SetAtomPosition(i, (float(xyz["X"][i]), float(xyz["Y"][i]), float(xyz["Z"][i])))
rw.AddConformer(conf)
rdDetermineBonds.DetermineConnectivity(rw, useHueckel=False)
rdkit_bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rw.GetBonds()]

# diff against mol.MolData['BONDS']/['ANGLES'] (offset by st_no), inject the
# missing entries into original_text's "Additional Bonds follow"/"Additional
# Bond Angles follow" sections -- see git history for this issue's exact
# diff/inject helper.
```
