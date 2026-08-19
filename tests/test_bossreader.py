"""BOSSReader-level unit tests (issue #23, parent map #19).

Feeds captured BOSS out/sum text (tests/fixtures/phenol, tests/fixtures/his
-- see each fixture's own README.md for exactly how it was captured)
directly into BOSSReader's individual parsing methods and asserts exact
values pulled from real BOSS output. No Docker/BOSS/license is needed to
run this file: everything these methods touch is plain text checked into
the repo, and none of get_bonds/get_angs/get_tors/get_tors_by_decl_idx/
get_atinfo/get_XYZ/get_pairs/get_charge ever run BOSS or read `self` (see
tests/conftest.py's BossSections helper for how instances are built).

phenol covers the baseline case: 2 leading dummy atoms, LigParGen's own
"C00"/"H0A"-style generated atom names, no ring-completion needed.

his covers the three non-baseline cases that actually broke
BOSS2OPENMM.py earlier in this repo's history (see
LigParGen/boss_common.py's docstrings): 3 leading dummy atoms, real
amino-acid atom names, and an imidazole ring whose bonds/angles/torsions
only get tabulated once its ring-closing internal coordinates are
re-declared as "Additional" entries.
"""
import pytest


# ---------------------------------------------------------------------
# get_atinfo
# ---------------------------------------------------------------------

def test_get_atinfo_phenol(phenol):
    ats = phenol.reader.get_atinfo(phenol.slice('ATMinit', 'ATMfinal'))
    assert len(ats) == 13
    assert ats[0].split()[:2] == ['3', 'C00']
    assert ats[-1].split()[:2] == ['15', 'H0C']


def test_get_atinfo_his(his):
    ats = his.reader.get_atinfo(his.slice('ATMinit', 'ATMfinal'))
    assert len(ats) == 29
    # First real atom is raw index 4 -- 3 leading dummies (DU1/DU2/DU3),
    # not LigParGen's usual 2.
    assert ats[0].split()[:2] == ['4', 'H1']
    assert ats[-1].split()[:2] == ['32', 'HT3']
    # Real amino-acid atom names, not LigParGen's "C00"-style generated
    # names.
    names = [line.split()[1] for line in ats]
    for expected_name in ('CB', 'CG', 'ND1', 'CE1', 'NE2', 'CD2'):
        assert expected_name in names


# ---------------------------------------------------------------------
# get_bonds
# ---------------------------------------------------------------------

def test_get_bonds_phenol(phenol):
    bonds = phenol.reader.get_bonds(phenol.slice('BNDinit', 'BNDfinal'))
    assert len(bonds['cl1']) == 13
    assert bonds['cl1'] == [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 8]
    assert bonds['cl2'] == [3, 4, 5, 6, 3, 6, 3, 4, 5, 7, 8, 9, 7]
    assert bonds['RIJ'][:3] == [1.4, 1.4, 1.4]
    assert bonds['KIJ'][:3] == [469.0, 469.0, 469.0]
    assert bonds['TIJ'][:3] == ['CA-CA', 'CA-CA', 'CA-CA']
    # The ring-closing bond (8-7) is present as its own row, with its own
    # (zero) delta -- it's the last Bond Stretching Parameters row BOSS
    # printed, not folded into an earlier duplicate.
    assert (bonds['cl1'][-1], bonds['cl2'][-1]) == (8, 7)


def test_get_bonds_his(his):
    bonds = his.reader.get_bonds(his.slice('BNDinit', 'BNDfinal'))
    assert len(bonds['cl1']) == 29
    assert bonds['cl1'][:3] == [5, 6, 7]
    assert bonds['cl2'][:3] == [4, 5, 5]
    assert bonds['RIJ'][:3] == [1.09, 1.522, 1.09]
    assert bonds['KIJ'][:3] == [340.0, 317.0, 340.0]
    assert bonds['TIJ'][:3] == ['CT-HC', 'C -CT', 'HC-CT']
    # A single-character trailing atom type ("H") shifts the fixed-width
    # TIJ slice by one column -- get_bonds() still keeps the full type
    # string (leading space and all), not a truncated one.
    assert bonds['TIJ'][-1] == ' NA-H'
    # The 9 ring bonds injected as "Additional Bonds" (see
    # tests/fixtures/his/README.md) are all present and tabulated with
    # real OPLS parameters, e.g. CG-ND1 (15-16) and CB-CG (14-15).
    ring_pairs = list(zip(bonds['cl1'], bonds['cl2']))
    assert (15, 16) in ring_pairs  # CG-ND1
    assert (14, 15) in ring_pairs  # CB-CG
    cg_nd1 = ring_pairs.index((15, 16))
    assert bonds['RIJ'][cg_nd1] == 1.394
    assert bonds['KIJ'][cg_nd1] == 410.0
    assert bonds['TIJ'][cg_nd1] == 'CV-NB'


# ---------------------------------------------------------------------
# get_angs
# ---------------------------------------------------------------------

def test_get_angs_phenol(phenol):
    angs = phenol.reader.get_angs(phenol.slice('ANGinit', 'ANGfinal'))
    assert len(angs['cl1']) == 19
    assert angs['cl1'][:3] == [3, 4, 5]
    assert angs['cl2'][:3] == [4, 5, 6]
    assert angs['cl3'][:3] == [5, 6, 7]
    assert angs['R'][:3] == [120.0, 120.0, 120.0]
    assert angs['K'][:3] == [63.0, 63.0, 63.0]


def test_get_angs_his(his):
    angs = his.reader.get_angs(his.slice('ANGinit', 'ANGfinal'))
    assert len(angs['cl1']) == 49
    assert angs['cl1'][:3] == [4, 6, 6]
    assert angs['cl2'][:3] == [5, 5, 5]
    assert angs['cl3'][:3] == [6, 7, 8]
    assert angs['R'][:3] == [109.5, 109.5, 109.5]
    assert angs['K'][:3] == [35.0, 35.0, 35.0]
    # The 13 ring angles injected as "Additional Bond Angles" are among
    # the tabulated rows.
    assert angs['cl1'][-3:] == [30, 30, 31]
    assert angs['cl2'][-3:] == [28, 28, 28]
    assert angs['cl3'][-3:] == [31, 32, 32]
    ring_angles = set(zip(angs['cl1'], angs['cl2'], angs['cl3']))
    assert (15, 16, 17) in ring_angles  # CG-ND1-CE1


# ---------------------------------------------------------------------
# get_tors / get_tors_by_decl_idx
# ---------------------------------------------------------------------

def test_get_tors_phenol(phenol):
    tors = phenol.reader.get_tors(phenol.slice('TORinit', 'TORfinal'))
    assert len(tors) == 32
    assert tors[0] == ['0.000', '7.250', '0.000', '0.000']
    assert tors[-1] == ['0.000', '5.000', '0.000', '0.000']


def test_get_tors_by_decl_idx_phenol(phenol):
    tors = phenol.reader.get_tors_by_decl_idx(
        phenol.slice('TORinit', 'TORfinal'))
    assert len(tors) == 32
    assert sorted(tors.keys()) == list(range(1, 33))
    assert tors[1] == ['0.000', '7.250', '0.000', '0.000']
    assert tors[10] == ['0.000', '2.060', '0.000', '0.000']
    assert tors[32] == ['0.000', '5.000', '0.000', '0.000']


def test_get_tors_by_decl_idx_his(his):
    tors = his.reader.get_tors_by_decl_idx(his.slice('TORinit', 'TORfinal'))
    assert len(tors) == 59
    assert sorted(tors.keys()) == list(range(1, 60))
    # Real, nonzero Fourier coefficients for backbone and ring torsions
    # (declared-dihedral index -> [V1, V2, V3, V4], pulled directly from
    # the captured "Dihedral ... Fourier Coefficients" table).
    assert tors[9] == ['-0.542', '0.435', '0.000', '0.000']
    assert tors[10] == ['-0.560', '-0.740', '0.349', '0.000']
    # Declared dihedral 41 is a ring torsion (its "Atom" column is 15,
    # CG) -- only tabulated with real coefficients because the ring's
    # bonds/angles were completed before this second BOSSReader pass; see
    # tests/fixtures/his/README.md.
    assert tors[41] == ['-1.282', '1.645', '-0.017', '0.000']
    # A genuinely all-zero row (an as-declared dihedral BOSS didn't find
    # nonzero Fourier terms for) is still a real row, not silently
    # dropped.
    assert tors[11] == ['0.000', '0.000', '0.000', '0.000']


# ---------------------------------------------------------------------
# get_XYZ
# ---------------------------------------------------------------------

def test_get_XYZ_phenol(phenol):
    xyz = phenol.reader.get_XYZ(phenol.slice('XYZinit', 'XYZfinal'))
    assert len(xyz) == 13
    assert list(xyz.columns) == ['at_num', 'X', 'Y', 'Z', 'at_symb']
    row0 = xyz.iloc[0]
    assert row0['at_num'] == 6
    assert row0['at_symb'] == 'C00'
    assert row0['X'] == pytest.approx(1.0)
    assert row0['Y'] == pytest.approx(1.0)
    assert row0['Z'] == pytest.approx(0.0)
    assert xyz.iloc[-1]['at_symb'] == 'H0C'


def test_get_XYZ_his(his):
    xyz = his.reader.get_XYZ(his.slice('XYZinit', 'XYZfinal'))
    assert len(xyz) == 29
    assert xyz.iloc[0]['at_symb'] == 'H1'
    # Non-"C00"-style names and real coordinates for a ring atom (NE2).
    ne2 = xyz.iloc[14]
    assert ne2['at_symb'] == 'NE2'
    assert ne2['at_num'] == 7
    assert ne2['X'] == pytest.approx(9.2613)
    assert ne2['Y'] == pytest.approx(4.22852)
    assert ne2['Z'] == pytest.approx(-1.01808)


# ---------------------------------------------------------------------
# get_pairs
# ---------------------------------------------------------------------

def test_get_pairs_phenol(phenol):
    pairs = phenol.reader.get_pairs(phenol.slice('PAIRinit', 'PAIRfinal'))
    assert len(pairs) == 46
    assert pairs[0] == '     1     4     1\n'
    assert pairs[-1] == '    12    13     1\n'


def test_get_pairs_his(his):
    pairs = his.reader.get_pairs(his.slice('PAIRinit', 'PAIRfinal'))
    assert len(pairs) == 328
    assert pairs[100] == '     6    17     1\n'
    assert pairs[200] == '    12    17     1\n'
    assert pairs[300] == '    21    25     1\n'


# ---------------------------------------------------------------------
# get_charge
# ---------------------------------------------------------------------

def test_get_charge_phenol(phenol):
    charge = phenol.reader.get_charge(phenol.charge_slice)
    assert charge == {
        'Reference-Solute': 0.0,
        '1st-Perturbed-Solute': 0.0,
        '2nd-Perturbed-Solute': 0.0,
    }


def test_get_charge_his(his):
    charge = his.reader.get_charge(his.charge_slice)
    assert charge == {
        'Reference-Solute': 0.0,
        '1st-Perturbed-Solute': 0.0,
        '2nd-Perturbed-Solute': 0.0,
    }


# ---------------------------------------------------------------------
# get_addihed (not in issue #23's required list, but it's the method
# directly responsible for MolData['ADD_DIHED'] -- the declared-dihedral
# sequence pair_declared_torsions() zips against TORSIONS_BY_DECL_IDX --
# and the his fixture's ring case is exactly what exercises it most.
# ---------------------------------------------------------------------

def test_get_addihed_his(his):
    addihed = his.reader.get_addihed(his.add_dihed_slice)
    assert len(addihed) == 39
    assert addihed[0] == ['10', '6', '5', '4']
    assert addihed[-1] == ['32', '28', '26', '29']
