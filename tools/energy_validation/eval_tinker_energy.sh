#!/bin/bash
# Run inside ligpargen-tinker:dev. Computes a single-point potential
# energy -- total AND decomposed by force type -- for a LigParGen-
# generated <resid>.new.xyz + <resid>.key pair, using TINKER's own
# `analyze` program (E option: total energy plus a per-component
# breakdown). No simulation, no minimization -- `analyze` itself is a
# pure single-point evaluator, matching what BOSS's own single-point
# "NEW E" and the OpenMM/GROMACS/NAMD/LAMMPS legs all represent.
#
# Usage: eval_tinker_energy.sh <resid>
# (expects <resid>.new.xyz and <resid>.key in the current directory --
# both written by gen_and_boss_energy.py via mainBOSS2TINKER)
#
# Prints:
#   TINKER_ENERGY_KCAL_PER_MOL <total>
#   TINKER_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
set -euo pipefail
resid="$1"
cd /tmp

# TINKER auto-associates a keyfile with an xyz file by matching basename
# (<base>.key next to <base>.xyz) -- mainBOSS2TINKER writes the keyfile
# as <resid>.key (matching the ORIGINAL <resid>.xyz OpenBabel produced,
# not the renamed <resid>.new.xyz it then edits in place), so it needs
# copying under the new.xyz basename for `analyze` to find it.
cp "${resid}.key" "${resid}.new.key"

# `echo ''` answers analyze's interactive "Enter Parameter File Name"
# prompt with a blank line -- the .key file already carries every force
# field parameter LigParGen generated (no external oplsaa.prm reference),
# so declining a separate parameter file is correct, not a workaround.
echo '' | analyze "${resid}.new.xyz" E > "${resid}_tinker.log" 2>&1
if ! grep -q "Total Potential Energy" "${resid}_tinker.log"; then
    echo "TINKER analyze failed:" >&2
    cat "${resid}_tinker.log" >&2
    exit 1
fi

python3 -c "
import re

with open('${resid}_tinker.log') as f:
    text = f.read()

# 'Bond Stretching'/'Angle Bending'/etc. as labels also appear earlier in
# the log, in an 'Additional ... Parameters' listing of every declared
# term -- restrict matching to the 'Energy Component Breakdown' block
# (and the 'Total Potential Energy' line just above it) so a label match
# there can't pick up an unrelated number from that earlier listing.
breakdown = text[text.index('Total Potential Energy'):]

def term(label):
    # Labels are followed by ' :' before the number in the 'Total
    # Potential Energy' line, but not in the per-component breakdown
    # rows -- skip any run of non-numeric characters rather than assuming
    # either form.
    m = re.search(re.escape(label) + r'[^\d.\-]*(-?[\d.]+)', breakdown)
    return float(m.group(1)) if m else 0.0

total = term('Total Potential Energy')
bond = term('Bond Stretching')
angle = term('Angle Bending')
# TINKER reports proper ('Torsional Angle') and improper ('Improper
# Torsion') dihedrals as separate categories when both are present --
# BOSS/OpenMM/GROMACS/LAMMPS/NAMD's 'torsion' term is their sum.
torsion = term('Torsional Angle') + term('Improper Torsion')
nonbonded = term('Van der Waals') + term('Charge-Charge')

print('TINKER_ENERGY_KCAL_PER_MOL %s' % total)
print('TINKER_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (bond, angle, torsion, nonbonded))
"
