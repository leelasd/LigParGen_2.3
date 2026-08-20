#!/bin/bash
# Single-point potential energy -- total AND decomposed by force type --
# for a LigParGen-generated <resid>.rtf + <resid>.prm + <resid>.pdb, using
# NAMD itself (via psfgen to build the PSF, then a genuine `run 0` -- not
# `minimize`, see ../README.md's "Single-point only" section for why that
# distinction matters).
#
# Unlike the OpenMM/GROMACS legs, this one is NOT Dockerized: NAMD is a
# licensed, proprietary binary supplied locally at runtime, never
# committed to this repo or baked into any image -- exactly the same
# constraint as BOSS itself (see docs/adr/0001). Point this script at your
# own NAMD install directory (the one containing the `namd3` and `psfgen`
# binaries) via NAMD_DIR.
#
# On macOS, both binaries will be Gatekeeper-quarantined the first time
# you run them (System Settings > Privacy & Security > Allow Anyway,
# separately for each binary) -- this only needs doing once per machine.
#
# Usage: NAMD_DIR=/path/to/namd/install ./eval_namd_energy.sh <resid>
# (expects <resid>.rtf, <resid>.prm, <resid>.pdb in the current directory,
# and this script's own directory to also contain psfgen_template.pgn and
# namd_sp_template.conf)
#
# Prints:
#   NAMD_ENERGY_KCAL_PER_MOL <total>
#   NAMD_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
set -euo pipefail
resid="$1"
here="$(cd "$(dirname "$0")" && pwd)"

if [ -z "${NAMD_DIR:-}" ]; then
    echo "ERROR: set NAMD_DIR to your local NAMD install directory (containing namd3, psfgen)." >&2
    exit 1
fi
if [ ! -x "$NAMD_DIR/psfgen" ] || [ ! -x "$NAMD_DIR/namd3" ]; then
    echo "ERROR: $NAMD_DIR does not contain executable psfgen/namd3." >&2
    exit 1
fi

sed "s/__RESID__/${resid}/g" "$here/psfgen_template.pgn" > "${resid}.pgn"
sed "s/__RESID__/${resid}/g" "$here/namd_sp_template.conf" > "${resid}_sp.conf"

"$NAMD_DIR/psfgen" "${resid}.pgn" > "${resid}_psfgen.log" 2>&1
if ! grep -q "psf file complete" "${resid}_psfgen.log"; then
    echo "psfgen failed:" >&2
    cat "${resid}_psfgen.log" >&2
    exit 1
fi

"$NAMD_DIR/namd3" +p1 "${resid}_sp.conf" > "${resid}_namd.log" 2>&1
if ! grep -q "^ENERGY:" "${resid}_namd.log"; then
    echo "NAMD run failed:" >&2
    cat "${resid}_namd.log" >&2
    exit 1
fi

# ETITLE:  TS  BOND  ANGLE  DIHED  IMPRP  ELECT  VDW  BOUNDARY  MISC  KINETIC  TOTAL  TEMP  POTENTIAL  TOTAL3  TEMPAVG
python3 -c "
line = [l for l in open('${resid}_namd.log') if l.startswith('ENERGY:')][0]
f = line.split()
bond, angle, dihed, imprp, elect, vdw = (float(f[i]) for i in (2, 3, 4, 5, 6, 7))
potential = float(f[13])
torsion = dihed + imprp
nonbonded = elect + vdw
print('NAMD_ENERGY_KCAL_PER_MOL %s' % potential)
print('NAMD_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (bond, angle, torsion, nonbonded))
"
