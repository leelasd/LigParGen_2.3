"""Shared fixture-loading helpers for the BOSSReader parsing-method test
layer (issue #23, parent map #19).

The whole point of this test layer is that it needs no Docker/BOSS/license
to run: every fixture under tests/fixtures/<name>/boss/{out,sum} is BOSS's
own captured, real output text, checked into the repo, and
BOSSReader.get_bonds()/get_angs()/get_tors()/get_tors_by_decl_idx()/
get_atinfo()/get_XYZ()/get_pairs()/get_charge() are pure functions of a
`data` argument (a list of lines) -- none of them touch `self`. So a test
can build a BOSSReader instance with BOSSReader.__new__(BOSSReader)
(skipping __init__, which would call Get_OPT() and actually try to run
BOSS) and call those methods directly against slices of the captured text.

Section boundaries are located via BOSSReader.find_boss_sections() --
the same banner-scanning logic get_ImpDat() itself uses -- rather than a
hand-copied duplicate of it, so this test layer can't silently desync
from what section boundaries BOSSReader actually computes. It's a plain
function of already-captured text (no Get_OPT()/real-BOSS call involved),
so it's safe to call directly here.
"""
import os

import pytest

from LigParGen.BOSSReader import BOSSReader, Refine_file, find_boss_sections

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')


class BossSections(object):
    """Captured BOSS out/sum text for one fixture molecule, pre-located
    into the same section boundaries BOSSReader.get_ImpDat() computes for
    its own parsing-method calls -- so tests can call those methods
    directly, on real captured text, without ever running BOSS."""

    def __init__(self, name):
        boss_dir = os.path.join(FIXTURES_DIR, name, 'boss')
        self.odat = Refine_file(os.path.join(boss_dir, 'out'))
        self.sdat = Refine_file(os.path.join(boss_dir, 'sum'))
        self.d = find_boss_sections(self.odat, self.sdat, zmat_name=name)
        self.reader = BOSSReader.__new__(BOSSReader)

    def slice(self, start_key, end_key):
        return self.odat[self.d[start_key]:self.d[end_key]]

    @property
    def charge_slice(self):
        # get_charge() expects a 4-line window: the 'Net Charge' banner
        # line itself, followed by the 3 solute charge lines.
        start = self.d['TotalQ']
        return self.odat[start:start + 4]

    @property
    def add_dihed_slice(self):
        return self.sdat[self.d['ADDinit']:self.d['ADDfinal']]


@pytest.fixture(scope='session')
def phenol():
    """LigParGen's own auto-generated Zmat baseline: 2 leading dummy
    atoms, 'C00'/'H0A'-style generated atom names, single ring but no
    Additional-declared ring bonds/angles needed (phenol's OH is the only
    rotatable bond). See tests/fixtures/phenol/README.md."""
    return BossSections('phenol')


@pytest.fixture(scope='session')
def his():
    """One of BOSS's own hand-crafted reference Zmats (Ac-His-NHMe): 3
    leading dummy atoms, real amino-acid atom names (not 'C00'-style), and
    an imidazole ring whose bonds/angles/torsions only got tabulated after
    injecting the ring's missing Additional Bonds/Angles -- see
    tests/fixtures/his/README.md for the full story."""
    return BossSections('his')


# NOTE: tests/test_integration_converters.py (issue #24) deliberately
# re-checks BOSS availability locally rather than importing it from here
# (see that file's own comment) -- so no boss_available()/BOSS_AVAILABLE
# helper is needed in this shared conftest; FIXTURES_DIR above already
# covers both this file's captured-text fixtures and that file's
# Docker/BOSS integration fixtures, since both live under tests/fixtures/.
