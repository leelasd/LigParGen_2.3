"""Run inside ligpargen-openmm:dev (or plain ligpargen:dev if you only need
the BOSS energy, not the OpenMM XML/PDB). Feeds a Zmatrix -- either one of
BOSS's own reference files (BOSSdir/molecules/*/*.z) or one pulled out of a
LigParGen-generated output zip -- through LigParGen's BOSSReader directly,
bypassing Converter.convert() so /tmp/out survives long enough to pull
BOSS's own single-point energy out of it before BOSSReader.cleanup() would
delete it (convert() calls cleanup() at the very end; nothing else in the
package does).

Usage (inside the container, BOSSdir already set by the image):
    python3 gen_and_boss_energy.py <resid> <path-to-zmatrix> [charge]

Prints (to stdout, one line each, machine-parseable):
    BOSS_ENERGY_KCAL_PER_MOL <total>
    BOSS_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>
    WROTE_FILES <resid>.xml <resid>.pdb

See ../README.md for the full methodology this is one half of.
"""
import sys
import os
import shutil
import pickle

sys.path.insert(0, '/app')
os.chdir('/tmp')

from LigParGen.BOSSReader import BOSSReader
from LigParGen.BOSS2OPENMM import mainBOSS2OPM
from LigParGen.BOSS2GMX import mainBOSS2GMX
from LigParGen.BOSS2CHARMM import mainBOSS2CHARMM
from LigParGen.BOSS2LAMMPS import mainBOSS2LAMMPS
from LigParGen.BOSS2TINKER import mainBOSS2TINKER
from LigParGen.BOSS2Q import mainBOSS2Q

resid = sys.argv[1]
zmat_src = sys.argv[2]
charge = int(sys.argv[3]) if len(sys.argv) > 3 else 0

shutil.copyfile(zmat_src, '%s.z' % resid)

# optim=0, lbcc=False: a genuine single-point evaluation of the geometry
# and parameters the Zmatrix already carries -- no re-optimization, no
# charge regeneration. See ../README.md's "Single-point only" section for
# why this matters.
mol = BOSSReader('%s.z' % resid, 0, charge, False)

boss_energy = None
boss_ebnd = boss_eang = boss_edih = boss_enb = None
with open('/tmp/out') as f:
    for line in f:
        if 'NEW E' in line:
            # e.g. "T(C) =   0.00  OLD E = 0.0D+00  NEW E = 0.66806231D+01"
            # BOSS prints Fortran double-precision exponent notation ('D'
            # instead of 'E') -- swap it before float() or this raises.
            val = line.split('NEW E =')[1].strip().split()[0].replace('D', 'E')
            boss_energy = float(val)
        if 'EBNDNE' in line:
            # " EBNDOL=   0.0000    EBNDNE=   0.2210    EANGOL=   0.0000    EANGNE=   0.0000"
            parts = line.replace('=', ' ').split()
            d = dict(zip(parts[0::2], parts[1::2]))
            boss_ebnd = float(d['EBNDNE'])
            boss_eang = float(d['EANGNE'])
        if 'EDIHNE' in line:
            parts = line.replace('=', ' ').split()
            d = dict(zip(parts[0::2], parts[1::2]))
            boss_edih = float(d['EDIHNE'])
            boss_enb = float(d['ENBNE'])

print('BOSS_ENERGY_KCAL_PER_MOL %s' % boss_energy)
print('BOSS_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (boss_ebnd, boss_eang, boss_edih, boss_enb))

# Also write the OpenMM XML+PDB, GROMACS itp+gro, CHARMM/NAMD rtf+prm,
# LAMMPS lmp, and TINKER new.xyz+key for this same BOSS-computed geometry,
# so eval_openmm_energy.py / eval_gromacs_energy.sh / eval_namd_energy.sh
# / eval_lammps_energy.sh / eval_tinker_energy.sh can evaluate them.
# Harmless (and cheap) to do even if you only wanted the BOSS numbers, or
# only some of the formats.
pickle.dump(mol, open('%s.pkl' % resid, 'wb'))
mainBOSS2OPM(resid, False)
mainBOSS2GMX(resid, False)
mainBOSS2CHARMM(resid)
mainBOSS2LAMMPS(resid)
mainBOSS2TINKER(resid)
mainBOSS2Q(resid)
print('WROTE_FILES %s.xml %s.pdb %s.itp %s.gro %s.rtf %s.prm %s.lmp %s.new.xyz %s.key %s.lib %s.Q.prm' % (resid, resid, resid, resid, resid, resid, resid, resid, resid, resid, resid))
