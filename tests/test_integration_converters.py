"""
Corpus-based integration test + cross-converter parity check (issue #24).

Runs a small corpus of real BOSS reference Zmats (tests/fixtures/*.z --
copied from BOSS's own molecules/peptide/ corpus, see tests/fixtures/README.md)
through the full BOSSReader -> all 8 BOSS2*.py converter pipeline, end to
end, against a real BOSS binary. This is deliberately the SLOW, real-BOSS
layer -- it needs Docker/a real BOSS install (like the existing
tests/fixtures/phenol fixture), and is separate from any faster
BOSSReader-level unit test layer working off captured text fixtures
(issue #23) that may live alongside this file.

Two things this test suite checks that a plain "does it crash" smoke test
would not:

1. Cross-converter parity: the same Zmat run through all 8 converters
   should agree on total real-atom count, per-atom element list, and total
   declared-torsion count, because every converter is built on the same
   shared LigParGen.boss_common.bossData()/pair_declared_torsions() helpers.
   This is exactly the kind of check that would have immediately flagged
   "7 of 8 converters still have the bug" the way this session's earlier
   manual audit did -- a converter with a stale/un-deduplicated copy of the
   old buggy logic would disagree with the rest on these counts (or crash
   outright) rather than silently agreeing.
2. Ring-bond/-angle completion: peptide/his.z (imidazole side chain) is a
   real BOSS reference Zmat whose ring-closing bonds/angles are declared in
   each ring atom's own NA/NB tree columns but were never given "Additional
   Bonds"/"Additional Bond Angles" entries -- so BOSS's own output never
   tabulates OPLS parameters for them. tests/fixtures/his_ring_complete.z
   is the same Zmat with those entries injected (see
   tests/fixtures/README.md for exactly how). Running both through BOSS
   confirms the completed fixture makes BOSS tabulate a strictly larger
   BONDS/ANGLES set -- the "missing bonded terms" bug this session's
   earlier work fixed, now pinned down by a real BOSS run instead of a
   one-off manual check.
"""
from __future__ import print_function

import glob
import os
import pickle
import re
import shutil

import pytest

# rdkit must be imported before openbabel: with the pinned rdkit/openbabel
# wheels this repo uses, importing openbabel's native extension first can
# segfault when rdkit.Chem is imported afterward (see Converter.py's own
# import order and its comment for the full story). BOSSReader -> mol_boss
# imports openbabel, so rdkit has to load first, here, before any other
# LigParGen submodule.
from rdkit import Chem  # noqa: F401
from rdkit.Chem import GetPeriodicTable


def _boss_available():
    bossdir = os.environ.get("BOSSdir")
    if not bossdir:
        return False
    return os.path.isfile(os.path.join(bossdir, "scripts", "xZCM1A"))


# Deliberately re-checked here (not imported from tests/conftest.py) so this
# file's collection/skip behavior doesn't depend on how pytest's rootdir
# import machinery resolves a `tests` package -- see tests/conftest.py for
# the same check shared with any sibling test module that wants it as a
# fixture instead.
BOSS_AVAILABLE = _boss_available()
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")

pytestmark = [
    pytest.mark.requires_boss,
    pytest.mark.skipif(
        not BOSS_AVAILABLE,
        reason="Requires a real BOSS binary: set $BOSSdir to a BOSS install "
        "with scripts/xZCM1A (e.g. run inside the ligpargen:dev Docker "
        "image built from this repo's Dockerfile -- see "
        "tests/fixtures/phenol/README.md for the general approach).",
    ),
]

# Imported after the rdkit-first guard above, and after the module-level
# skip is registered, so a no-BOSS collection run never needs a working
# LigParGen/mol_boss (openbabel) import to just report "skipped".
if BOSS_AVAILABLE:
    from LigParGen.BOSSReader import BOSSReader
    from LigParGen.boss_common import bossData, pair_declared_torsions
    from LigParGen.BOSS2OPENMM import mainBOSS2OPM
    from LigParGen.BOSS2Q import mainBOSS2Q
    from LigParGen.BOSS2XPLOR import mainBOSS2XPLOR
    from LigParGen.BOSS2CHARMM import mainBOSS2CHARMM
    from LigParGen.BOSS2GMX import mainBOSS2GMX
    from LigParGen.BOSS2LAMMPS import mainBOSS2LAMMPS
    from LigParGen.BOSS2DESMOND import mainBOSS2DESMOND
    from LigParGen.BOSS2TINKER import mainBOSS2TINKER

# (fixture file, resname, human description)
FIXTURES = [
    ("his.z", "HRW", "his.z, raw BOSS reference Zmat (imidazole ring, ring bonds/angles NOT completed)"),
    ("his_ring_complete.z", "HRC", "his.z with ring-closing Additional Bonds/Bond Angles injected"),
    ("ala.z", "ALA", "ala.z, plain baseline peptide fragment (no ring)"),
]
FIXTURE_IDS = [resname for _, resname, _ in FIXTURES]


# --------------------------------------------------------------------------
# Helpers: parse each converter's own output file for atom/torsion counts.
# Deliberately independent, format-specific parsers (not just re-reading
# MolData) -- the whole point is to check what each converter actually
# WROTE, not what BOSSReader parsed.
# --------------------------------------------------------------------------

def _real_atom_count_from_zmat(path):
    """Ground truth atom count, parsed directly from the fixture file
    itself: real atoms are Zmat atom lines whose type column is > 1 (BOSS's
    own dummy-atom convention -- dummies are type -1), independent of any
    BOSS run at all."""
    lines = open(path).readlines()
    stop = next(i for i, l in enumerate(lines) if "Geometry Variations follow" in l)
    return sum(1 for l in lines[1:stop] if l.split() and float(l.split()[2]) > 1)


def _section_data_lines(text, header, stop_headers):
    """Count non-blank, non-comment lines between a `header` marker line and
    the next line containing any of `stop_headers` (or EOF)."""
    lines = text.splitlines()
    start = None
    for i, l in enumerate(lines):
        if header in l:
            start = i + 1
            break
    if start is None:
        return 0
    n = 0
    for l in lines[start:]:
        if any(h in l for h in stop_headers):
            break
        s = l.strip()
        if not s or s[0] in "#;!":
            continue
        n += 1
    return n


def atoms_from_gro(text):
    return int(text.splitlines()[1].strip())


def atoms_from_itp(text):
    return _section_data_lines(text, "[ atoms ]", ["[ bonds ]"])


def atoms_from_rtf(text):
    return len([l for l in text.splitlines() if l.startswith("ATOM ")])


def atoms_from_top(text):
    return len([l for l in text.splitlines() if l.startswith("ATOM ")])


def atoms_from_lib(text):
    return _section_data_lines(text, "[atoms]", ["[bonds]"])


def atoms_from_lmp(text):
    return int(re.search(r"(\d+)\s+atoms\b", text).group(1))


def atoms_from_cms(text):
    return int(re.search(r"m_atom\[(\d+)\]", text).group(1))


def atoms_from_key(text):
    return len([l for l in text.splitlines() if re.match(r"^atom\s+\d+", l.strip())])


def atoms_from_pdb(text):
    return len([l for l in text.splitlines() if l.startswith("ATOM")])


def torsions_from_lmp(text):
    dihed = int(re.search(r"(\d+)\s+dihedrals\b", text).group(1))
    imp = int(re.search(r"(\d+)\s+impropers\b", text).group(1))
    return dihed + imp


def torsions_from_cms(text):
    return int(re.search(r"ffio_dihedrals\[(\d+)\]", text).group(1))


def elements_from_pdb(text):
    # PDB columns 77-78 (0-indexed 76:78): element symbol, written by
    # BOSS2OPENMM.py's pdb_prep() from bossData()'s periodic-table lookup.
    return [l[76:78].strip() for l in text.splitlines() if l.startswith("ATOM")]


def elements_from_cms(text, n_atoms):
    """DESMOND's .cms m_atom[] block encodes each atom's atomic number
    (i_m_atomic_number, the 7th field of each data row) independently of
    the PDB's own element column -- convert back to a symbol via RDKit's
    periodic table for comparison."""
    pt = GetPeriodicTable()
    elems = []
    for l in text.splitlines():
        m = re.match(
            r"\s*\d+\s+\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+\d+\s+\"MOL\"\s+(\d+)\s+\"",
            l,
        )
        if m:
            elems.append(pt.GetElementSymbol(int(m.group(1))))
            if len(elems) == n_atoms:
                break
    return elems


# --------------------------------------------------------------------------
# Core pipeline runner: BOSSReader -> pickle -> all 8 converters, exactly
# mirroring LigParGen.Converter.convert()'s own zmat-input code path (see
# Converter.py: chdir to /tmp, BOSSReader(..., lbcc=False) for a supplied
# Zmat, pickle the mol, then run each mainBOSS2*() in the same order).
# --------------------------------------------------------------------------

@pytest.fixture(scope="module", params=FIXTURES, ids=FIXTURE_IDS)
def converted(request):
    fname, resname, _desc = request.param
    src = os.path.join(FIXTURES_DIR, fname)
    assert os.path.isfile(src), "missing fixture: %s" % src

    os.chdir("/tmp")
    for stale in glob.glob("/tmp/%s.*" % resname):
        os.remove(stale)
    shutil.copyfile(src, "/tmp/%s.z" % resname)

    mol = BOSSReader("%s.z" % resname, 0, 0, False)

    # Reference values computed directly from BOSSReader's own parsed
    # MolData / the shared boss_common helpers every converter is built on
    # -- NOT re-derived per converter, so these are a genuine independent
    # ground truth for the parity checks below.
    types, Qs, num2opls, st_no, num2typ2symb, num2pqrtype = bossData(mol)
    ref_elements = [num2typ2symb[i][3] for i in range(len(types))]

    ats = []
    for line in mol.MolData["ATOMS"][3:]:
        f = line.split()
        ats.append([int(f[0]), int(f[4]), int(f[6]), int(f[8])])
    for line in mol.MolData["ADD_DIHED"]:
        ats.append([int(x) for x in line])
    paired_ats, _paired_dhd = pair_declared_torsions(mol, ats)

    pickle.dump(mol, open(resname + ".pkl", "wb"))
    # Same order as LigParGen.Converter.convert() (TINKER last -- it reads
    # back the .pdb file OPENMM's converter just wrote).
    mainBOSS2OPM(resname, False)
    mainBOSS2Q(resname, False)
    mainBOSS2XPLOR(resname, False)
    mainBOSS2CHARMM(resname, False)
    mainBOSS2GMX(resname, False)
    mainBOSS2LAMMPS(resname, False)
    mainBOSS2DESMOND(resname, False)
    mainBOSS2TINKER(resname, False)

    files = {}
    for ext in ("gro", "itp", "rtf", "prm", "top", "param", "lib", "Q.prm", "lmp", "cms", "key", "pdb", "xml"):
        path = "/tmp/%s.%s" % (resname, ext)
        if os.path.isfile(path):
            files[ext] = open(path).read()

    yield {
        "resname": resname,
        "fname": fname,
        "mol": mol,
        "n_atoms_ref": len(types),
        "n_atoms_zmat": _real_atom_count_from_zmat(src),
        "n_bonds": len(mol.MolData["BONDS"]["cl1"]),
        "n_angles": len(mol.MolData["ANGLES"]["cl1"]),
        "n_torsions_ref": len(paired_ats),
        "ref_elements": ref_elements,
        "files": files,
    }

    mol.cleanup()


# --------------------------------------------------------------------------
# 1. Basic "runs clean, sane counts" integration test.
# --------------------------------------------------------------------------

def test_pipeline_runs_without_exceptions(converted):
    # If we get here, BOSSReader + all 8 mainBOSS2*() converters already ran
    # without raising (the `converted` fixture itself is where that
    # happens) -- this test just asserts the counts they all agreed on are
    # sane, not garbage/zero/absurdly large.
    n = converted["n_atoms_ref"]
    assert n > 0
    assert n == converted["n_atoms_zmat"], (
        "BOSSReader's own tabulated atom count (%d) disagrees with the real "
        "(non-dummy) atom count parsed directly from the fixture Zmat (%d)"
        % (n, converted["n_atoms_zmat"])
    )
    assert 0 < converted["n_bonds"] < n * 6
    assert 0 < converted["n_angles"] < n * 12
    assert 0 < converted["n_torsions_ref"] < n * 20
    # every expected output file exists and is non-empty
    expected_exts = {"gro", "itp", "rtf", "prm", "top", "param", "lib", "Q.prm", "lmp", "cms", "key", "pdb", "xml"}
    missing = expected_exts - set(converted["files"])
    assert not missing, "converters did not produce: %s" % sorted(missing)
    for ext, text in converted["files"].items():
        assert text.strip(), "%s.%s is empty" % (converted["resname"], ext)


# --------------------------------------------------------------------------
# 2. Cross-converter parity check: same Zmat, 8 converters -- atom count,
#    element list, and (where the format encodes it directly) torsion
#    count must all agree, since every converter shares the same
#    boss_common.bossData()/pair_declared_torsions() plumbing.
# --------------------------------------------------------------------------

def test_cross_converter_atom_count_parity(converted):
    files = converted["files"]
    counted = {
        "openmm(.pdb)": atoms_from_pdb(files["pdb"]),
        "gromacs(.gro)": atoms_from_gro(files["gro"]),
        "gromacs(.itp)": atoms_from_itp(files["itp"]),
        "charmm(.rtf)": atoms_from_rtf(files["rtf"]),
        "xplor(.top)": atoms_from_top(files["top"]),
        "q(.lib)": atoms_from_lib(files["lib"]),
        "lammps(.lmp)": atoms_from_lmp(files["lmp"]),
        "desmond(.cms)": atoms_from_cms(files["cms"]),
        "tinker(.key)": atoms_from_key(files["key"]),
    }
    ref = converted["n_atoms_ref"]
    bad = {k: v for k, v in counted.items() if v != ref}
    assert not bad, (
        "atom count disagreement for %s: reference=%d, mismatches=%s "
        "(this is exactly the class of bug a stale/un-deduplicated "
        "converter copy would produce)" % (converted["resname"], ref, bad)
    )


def test_cross_converter_torsion_count_parity(converted):
    # LAMMPS and DESMOND both encode the raw declared/paired-torsion count
    # (Proper + Improper together) as a literal header integer, computed
    # via the same shared pair_declared_torsions() every converter calls --
    # so these two must match the reference exactly.
    files = converted["files"]
    lmp_total = torsions_from_lmp(files["lmp"])
    cms_total = torsions_from_cms(files["cms"])
    ref = converted["n_torsions_ref"]
    assert lmp_total == ref, "LAMMPS declared torsion count %d != reference %d" % (lmp_total, ref)
    assert cms_total == ref, "DESMOND declared torsion count %d != reference %d" % (cms_total, ref)

    # The remaining formats deduplicate/collapse torsions by atom-type name
    # (not one row per declared quadruple), so an exact count match isn't
    # meaningful -- but if the reference has torsions at all, each format's
    # own torsion/dihedral section must be non-empty.
    if ref > 0:
        assert "<Proper " in files["xml"] or "<Improper " in files["xml"]
        assert "[ dihedrals ]" in files["itp"]
        assert _section_data_lines(files["prm"], "DIHEDRAL", ["IMPROPER"]) > 0
        assert "DIHEdral" in files["param"] or "IMPRoper" in files["param"]
        assert _section_data_lines(files["Q.prm"], "[torsions]", ["[impropers]"]) > 0
        assert "torsion " in files["key"] or "imptors " in files["key"]


def test_cross_converter_element_list_parity(converted):
    ref = converted["ref_elements"]
    pdb_elems = elements_from_pdb(converted["files"]["pdb"])
    cms_elems = elements_from_cms(converted["files"]["cms"], len(ref))
    assert pdb_elems == ref, "OpenMM PDB element column disagrees with reference element list"
    assert cms_elems == ref, "DESMOND atomic-number-derived element list disagrees with reference element list"


# --------------------------------------------------------------------------
# 3. Ring-bond/-angle completion regression check: his.z (raw) vs
#    his_ring_complete.z (Additional Bonds/Bond Angles injected) -- see
#    tests/fixtures/README.md for exactly how the completed fixture was
#    built.
# --------------------------------------------------------------------------

def test_ring_completion_adds_tabulated_bonded_terms():
    raw = [c for c in FIXTURES if c[1] == "HRW"][0]
    complete = [c for c in FIXTURES if c[1] == "HRC"][0]

    def run(fname, resname):
        src = os.path.join(FIXTURES_DIR, fname)
        os.chdir("/tmp")
        for stale in glob.glob("/tmp/%s.*" % resname):
            os.remove(stale)
        shutil.copyfile(src, "/tmp/%s.z" % resname)
        mol = BOSSReader("%s.z" % resname, 0, 0, False)
        result = (len(mol.MolData["BONDS"]["cl1"]), len(mol.MolData["ANGLES"]["cl1"]))
        mol.cleanup()
        return result

    raw_bonds, raw_angles = run(raw[0], raw[1])
    complete_bonds, complete_angles = run(complete[0], complete[1])

    assert complete_bonds > raw_bonds, (
        "his_ring_complete.z's injected 'Additional Bonds' did not make BOSS "
        "tabulate more bonds than the raw fixture (raw=%d, complete=%d) -- "
        "the ring-completion fixture may no longer be doing what it's for"
        % (raw_bonds, complete_bonds)
    )
    assert complete_angles > raw_angles, (
        "his_ring_complete.z's injected 'Additional Bond Angles' did not "
        "make BOSS tabulate more angles than the raw fixture (raw=%d, "
        "complete=%d)" % (raw_angles, complete_angles)
    )
