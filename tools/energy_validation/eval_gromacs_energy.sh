#!/bin/bash
# Run inside ligpargen-gmx:dev. Computes a single-point potential energy --
# total AND decomposed by force type -- for a LigParGen-generated
# <resid>.itp + <resid>.gro pair, using GROMACS itself. No simulation, no
# minimization (nsteps=0): evaluate the force field at the given
# (BOSS-optimized) geometry, matching what BOSS's own single-point "NEW E"
# (see gen_and_boss_energy.py) and eval_openmm_energy.py's OpenMM
# evaluation both represent.
#
# Usage: eval_gromacs_energy.sh <resid>
# (expects <resid>.itp and <resid>.gro in the current directory, and this
# script's own directory to also contain vacuum_sp.mdp)
#
# Prints:
#   GROMACS_ENERGY_KCAL_PER_MOL <total>
#   GROMACS_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
#
# See ../README.md for the full methodology, including why this needs a
# hand-written [ defaults ] directive and an enlarged box (two real,
# non-obvious gotchas hit while first building this).
set -euo pipefail
resid="$1"
mdp="$(dirname "$0")/vacuum_sp.mdp"
cd /tmp

# GROMACS's itp files (as LigParGen writes them) are a self-contained
# per-molecule fragment (atomtypes + moleculetype + atoms/bonds/angles/
# dihedrals) -- there's no system-level .top wrapping it, since LigParGen
# never assumes what box/solvent/ion setup the itp will end up in. Build
# the minimal one needed for a single-molecule vacuum evaluation.
# [ defaults ] is REQUIRED before [ atomtypes ] can appear anywhere in the
# assembled topology (grompp: "Invalid order for directive atomtypes" if
# it's missing) -- it's conventionally supplied once, system-wide, by
# whatever includes a molecule's itp, not by the itp itself.
cat > system.top <<EOF
[ defaults ]
; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
1         3          yes        0.5      0.5

#include "${resid}.itp"

[ system ]
${resid} single molecule vacuum

[ molecules ]
${resid}    1
EOF

# Enlarge the box far past the molecule's own extent (LigParGen's own .gro
# output uses a nominal 1x1x1nm box, too small to avoid a molecule seeing
# its own periodic images within any reasonable cutoff) -- pairs with
# vacuum_sp.mdp's pbc=xyz + 2nm cutoff to approximate an isolated molecule.
python3 -c "
lines = open('${resid}.gro').readlines()
lines[-1] = '   5.00000   5.00000   5.00000\n'
open('${resid}.gro', 'w').writelines(lines)
"

gmx grompp -f "$mdp" -c "${resid}.gro" -p system.top -o sp.tpr -maxwarn 10 > grompp.log 2>&1
gmx mdrun -s sp.tpr -deffnm sp -nt 1 > mdrun.log 2>&1
echo -e "1\n2\n3\n4\n5\n6\n7\n8\n9\n0" | gmx energy -f sp.edr -o energy.xvg > energy.log 2>&1

python3 -c "
import re

terms = {}
with open('energy.log') as f:
    for line in f:
        m = re.match(
            r'^(Bond|Angle|Ryckaert-Bell\.|Per\. Imp\. Dih\.|LJ-14|Coulomb-14|LJ \(SR\)|Coulomb \(SR\)|Potential)\s+([-\d.]+)',
            line)
        if m:
            terms[m.group(1)] = float(m.group(2))

def kj2kcal(x):
    return x / 4.184

bond = kj2kcal(terms.get('Bond', 0.0))
angle = kj2kcal(terms.get('Angle', 0.0))
# GROMACS reports proper (Ryckaert-Bellemans form) and improper (periodic
# form) torsions as two separate energy categories -- BOSS/OpenMM's
# 'torsion' term is their sum.
torsion = kj2kcal(terms.get('Ryckaert-Bell.', 0.0) + terms.get('Per. Imp. Dih.', 0.0))
nonbonded = kj2kcal(
    terms.get('LJ-14', 0.0) + terms.get('Coulomb-14', 0.0) +
    terms.get('LJ (SR)', 0.0) + terms.get('Coulomb (SR)', 0.0))
total = kj2kcal(terms.get('Potential', 0.0))

print('GROMACS_ENERGY_KCAL_PER_MOL %s' % total)
print('GROMACS_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (bond, angle, torsion, nonbonded))
"
