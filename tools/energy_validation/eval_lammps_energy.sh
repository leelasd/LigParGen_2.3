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
# box workaround is needed here.
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
