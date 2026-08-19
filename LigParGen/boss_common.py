"""
Shared, bug-fixed helpers used by every BOSS2*.py converter.

This module exists to de-duplicate two pieces of logic that used to be
copy-pasted, with subtly diverging bugs, into each of the 8 BOSS2*.py
converters (BOSS2OPENMM, BOSS2GMX, BOSS2CHARMM, BOSS2LAMMPS, BOSS2TINKER,
BOSS2XPLOR, BOSS2Q, BOSS2DESMOND):

1. bossData(molecule_data) -- builds the per-atom type/element/mass/charge
   table (num2typ2symb / num2opls / num2pqrtype) and the st_no Zmat-index
   offset, shared verbatim by every converter.
2. pair_declared_torsions(molecule_data, ats) -- pairs a converter's own
   list of declared dihedral quadruples to BOSS's Fourier-coefficient rows
   by BOSS's own declaration-order index, instead of a fragile positional
   zip/assert that breaks the moment BOSS skips a row.

See BOSS2OPENMM.py's git history for the original bug reports/fixes this
was extracted from.
"""

from collections import OrderedDict
from LigParGen.BOSSReader import bossPdbAtom2Element, bossElement2Mass
from rdkit import Chem

_PERIODIC_TABLE = Chem.GetPeriodicTable()


def bossData(molecule_data):
    ats_file = molecule_data.MolData['ATOMS']
    # Elements come from MolData['XYZ']'s own atomic-number column, not
    # bossPdbAtom2Element's name-based guess (strip trailing char, drop
    # digits): that heuristic assumes LigParGen's own "C00"/"H0A"-style
    # generated atom names and mistypes real atom names this package
    # doesn't generate itself, e.g. an amino-acid-style "CH1" comes out as
    # element "CH" and crashes bossElement2Mass. The atomic number BOSS
    # itself reports is unambiguous regardless of naming convention.
    xyz = molecule_data.MolData['XYZ']
    types = []
    for i in enumerate(ats_file):
        types.append([i[1].split()[1], 'opls_' + i[1].split()[2]])
    # st_no offsets a bonded-pair's raw Zmat atom index down to a 0-based
    # index into types/Qs/num2opls (which start at the first REAL atom).
    # LigParGen's own auto-generated Zmats always place exactly 2 leading
    # dummy atoms, so the first real atom's raw index is always 3 -- but
    # BOSS's own reference Zmat library is inconsistent about this (some
    # files use 2 leading dummies, some 3, some place the dummies after
    # the first real atom instead of before), so hardcoding 3 silently
    # mis-indexes bonded pairs for anything that isn't the 2-leading-dummy
    # case, and can even index out of range (reproduced on his.z).
    # Reading it from the first real atom's own raw index is correct for
    # both conventions.
    st_no = int(ats_file[0].split()[0])
    Qs = molecule_data.MolData['Q_LJ']
    assert len(Qs) == len(types), 'Please check the at_info and Q_LJ_dat files'
    assert len(xyz) == len(types), 'Please check the at_info and XYZ data'
    num2typ2symb = {i: types[i] for i in range(len(Qs))}
    for i in range(len(Qs)):
        elem = _PERIODIC_TABLE.GetElementSymbol(int(xyz['at_num'][i]))
        num2typ2symb[i].append(elem + num2typ2symb[i][1][-3:])
        num2typ2symb[i].append(elem)
        num2typ2symb[i].append(bossElement2Mass(elem))
        num2typ2symb[i].append(Qs[i][0])
    num2opls = {}
    for i in num2typ2symb.keys():
        num2opls[i] = num2typ2symb[i][2]
    num2pqrtype = OrderedDict(num2typ2symb)
    for i in range(len(Qs)):
        num2pqrtype[i].append(Qs[i][1])
        num2pqrtype[i].append(Qs[i][2])
    return (types, Qs, num2opls, st_no, num2typ2symb, num2pqrtype)


def pair_declared_torsions(molecule_data, ats):
    """Pair a converter's own list of declared dihedral quadruples (`ats`,
    Variable Dihedrals followed by Additional Dihedrals, in declaration
    order) to BOSS's Fourier-coefficient rows.

    BOSS's own "Angle" column in the Fourier Coefficients table is a
    1-based index into that same declared-dihedral sequence -- see
    BOSSReader.get_ImpDat's TORinit/TORfinal slice. That table is NOT
    guaranteed one row per declared quadruple: BOSS can skip a quadruple
    it can't match to a known torsion type (confirmed: this happens for
    some of BOSS's own re-derived Additional Dihedrals once a molecule's
    ring bonds/angles are completed). Pairing by this declared-order index
    -- MolData['TORSIONS_BY_DECL_IDX'] -- instead of assuming row N always
    means declared quadruple N is what keeps this correct when rows are
    skipped, rather than either mispairing silently or asserting a
    (previously exact-count-required) match.

    Returns (paired_ats, paired_dhd): paired_ats is the subset of `ats`
    BOSS actually tabulated coefficients for (in the same order as
    paired_dhd); paired_dhd is each one's raw [V1, V2, V3, V4] Fourier
    coefficients as floats (still in kcal/mol, unconverted).
    """
    tors_by_idx = molecule_data.MolData['TORSIONS_BY_DECL_IDX']
    paired_ats, paired_dhd = [], []
    for decl_idx, quad in enumerate(ats, start=1):
        coeffs = tors_by_idx.get(decl_idx)
        if coeffs is None:
            continue
        paired_ats.append(quad)
        paired_dhd.append([float(v) for v in coeffs])
    return paired_ats, paired_dhd
