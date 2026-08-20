#!/bin/bash
# Run inside ligpargen-lammps:dev. Computes a single-point potential
# energy -- total AND decomposed by force type -- for a LigParGen-
# generated <resid>.lmp data file, using LAMMPS itself (`lmp`). No
# simulation, no minimization (`run 0`): evaluate the force field at the
# given (BOSS-optimized) geometry, matching what BOSS's own single-point
# "NEW E" and the OpenMM/GROMACS/NAMD legs all represent.
#
# Free (non-periodic) boundaries + a 100 A cutoff approximate BOSS/
# OpenMM/NAMD's vacuum evaluation directly -- unlike modern GROMACS,
# LAMMPS actually supports true non-periodic boundaries, so no enlarged-
# box workaround is needed for the CUTOFF -- but the .lmp file's own box
# (BOSS2LAMMPS.py writes xlo/ylo/zlo as the exact coordinate minimum, by
# construction placing at least one atom exactly ON the box's lower edge)
# can make LAMMPS's domain decomposition silently fail to assign that
# atom ("Did not assign all atoms correctly"), even with free boundaries.
# Confirmed directly on a real molecule (phenol) -- padding every bound
# by 1 A before running fixes it. Not a BOSS2LAMMPS.py bug (it's a
# perfectly valid LAMMPS data file); a property of how domain
# decomposition -- unrelated to the boundary style -- treats an atom
# sitting exactly at a box face.
#
# Usage: eval_lammps_energy.sh <resid>
# (expects <resid>.lmp in the current directory, and this script's own
# directory to also contain lammps_sp_template.in)
#
# Prints:
#   LAMMPS_ENERGY_KCAL_PER_MOL <total>
#   LAMMPS_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
set -euo pipefail
resid="$1"
here="$(cd "$(dirname "$0")" && pwd)"
cd /tmp

python3 -c "
import re
with open('${resid}.lmp') as f:
    lines = f.readlines()
for i, l in enumerate(lines):
    m = re.match(r'\s*([-\d.]+)\s+([-\d.]+)\s+(xlo xhi|ylo yhi|zlo zhi)', l)
    if m:
        lo, hi, tag = float(m.group(1)), float(m.group(2)), m.group(3)
        lines[i] = '  %.6f  %.6f %s\n' % (lo - 1.0, hi + 1.0, tag)
with open('${resid}.lmp', 'w') as f:
    f.writelines(lines)
"

sed "s/__RESID__/${resid}/g" "$here/lammps_sp_template.in" > "${resid}_sp.in"
lmp -in "${resid}_sp.in" > "${resid}_lammps.log" 2>&1

python3 -c "
import re

with open('${resid}_lammps.log') as f:
    lines = f.readlines()

header_idx = next(i for i, l in enumerate(lines) if l.startswith('Step '))
values = [float(x) for x in lines[header_idx + 1].split()]
cols = lines[header_idx].split()
d = dict(zip(cols, values))

bond = d['E_bond']
angle = d['E_angle']
torsion = d['E_dihed'] + d['E_impro']
nonbonded = d['E_vdwl'] + d['E_coul'] + d.get('E_long', 0.0)
total = d['PotEng']

print('LAMMPS_ENERGY_KCAL_PER_MOL %s' % total)
print('LAMMPS_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (bond, angle, torsion, nonbonded))
"
