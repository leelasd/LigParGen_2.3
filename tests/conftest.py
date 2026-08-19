"""
Shared pytest configuration for LigParGen's test suite.

boss_available()/BOSS_AVAILABLE are used by the real-BOSS integration tests
(tests/test_integration_converters.py, issue #24) to skip cleanly when no
real BOSS binary is present (e.g. a plain local `pytest` run outside the
ligpargen:dev Docker image) rather than failing on a missing $BOSSdir.
This is deliberately independent of any BOSSReader-level *unit* test layer
(issue #23) that may exist alongside this file -- that layer works off
captured text fixtures and never touches $BOSSdir or a real BOSS install.
"""
import os


def boss_available():
    bossdir = os.environ.get("BOSSdir")
    if not bossdir:
        return False
    return os.path.isfile(os.path.join(bossdir, "scripts", "xZCM1A"))


BOSS_AVAILABLE = boss_available()
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
