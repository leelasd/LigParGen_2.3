#!/bin/bash
# Run inside ligpargen-q:dev. Computes a single-point potential energy --
# total AND decomposed by force type -- for a LigParGen-generated
# <resid>.lib + <resid>.Q.prm pair, using Q's own Qprep6 (topology build)
# and Qdyn6 (energy evaluation).
#
# Qdyn6 has no true steps=0 mode ("Need at least one step of dynamics")
# and refuses temperature=0 ("No dynamics at zero temperature!"), so this
# uses steps=1 at a near-zero temperature (0.001 K) and reads the "Energy
# summary at step      0" block Qdyn6 prints BEFORE that one integration
# step runs -- a genuine single-point evaluation of the input geometry,
# not an artifact of the one tiny step (confirmed identical to the "FINAL
# Energy summary" block in every run checked here).
#
# Usage: eval_q_energy.sh <resid>
# (expects <resid>.lib, <resid>.Q.prm, <resid>.pdb in the current
# directory, and this script's own directory to also contain
# q_prep_template.inp and q_sp_template.inp)
#
# Prints:
#   Q_ENERGY_KCAL_PER_MOL <total>
#   Q_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
set -euo pipefail
resid="$1"
here="$(cd "$(dirname "$0")" && pwd)"
cd /tmp

# Qprep6's PDB reader doesn't understand TER/CONECT/END/REMARK lines --
# confirmed directly, their presence makes it miscount "0 molecules"
# instead of 1, which corrupts topology assembly downstream
# ("Inconsistent molecule/residue start atoms", every bond/angle/torsion
# count coming out 0). Strip to ATOM lines only.
grep '^ATOM' "${resid}.pdb" > "${resid}_stripped.pdb"

sed "s/__RESID__/${resid}/g" "$here/q_prep_template.inp" > "${resid}_prep.inp"
Qprep6 < "${resid}_prep.inp" > "${resid}_qprep.log" 2>&1
if ! grep -q "Topology successfully generated" "${resid}_qprep.log"; then
    echo "Qprep6 failed:" >&2
    cat "${resid}_qprep.log" >&2
    exit 1
fi

sed "s/__RESID__/${resid}/g" "$here/q_sp_template.inp" > "${resid}_qsp.inp"
Qdyn6 "${resid}_qsp.inp" > "${resid}_qdyn.log" 2>&1
if ! grep -q "terminated normally" "${resid}_qdyn.log"; then
    echo "Qdyn6 failed:" >&2
    cat "${resid}_qdyn.log" >&2
    exit 1
fi

python3 -c "
import re

with open('${resid}_qdyn.log') as f:
    text = f.read()

m = re.search(r'Energy summary at step\s+0.*?\n(.*?)\nsolute\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', text, re.S)
el, vdw, bond, angle, torsion, improper = (float(x) for x in m.groups()[1:])

m2 = re.search(r'Energy summary at step\s+0.*?SUM\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', text, re.S)
total = float(m2.group(2))

nonbonded = el + vdw
torsion_total = torsion + improper

print('Q_ENERGY_KCAL_PER_MOL %s' % total)
print('Q_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (bond, angle, torsion_total, nonbonded))
"
