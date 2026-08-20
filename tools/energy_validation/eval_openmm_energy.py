"""Run inside ligpargen-openmm:dev. Computes a single-point potential
energy -- total AND decomposed by force type -- for a LigParGen-generated
<resid>.pdb + <resid>.xml pair, using OpenMM itself. No simulation, no
minimization: evaluate the force field at the given (BOSS-optimized)
geometry, matching what BOSS's own single-point "NEW E" (and its
EBNDNE/EANGNE/EDIHNE/ENBNE breakdown, see gen_and_boss_energy.py)
represents.

Total-energy agreement alone can hide compensating errors in different
terms (see https://github.com/openmm/openmm/issues/1463) -- assign each
Force a distinct force group and query each group's energy separately,
the same technique discussed there, rather than pull in ParmEd (which
does the same thing internally via its own energy_decomposition_system)
for a single System we already have direct access to build ourselves.

Usage: python3 eval_openmm_energy.py <resid>
(expects <resid>.pdb and <resid>.xml in the current directory)

Prints:
    OPENMM_ENERGY_KCAL_PER_MOL <total>
    OPENMM_TERMS bond=<..> angle=<..> torsion=<..> nonbonded=<..>

See ../README.md for the full methodology.
"""
import sys
import os
from openmm.app import PDBFile, ForceField, NoCutoff
from openmm import Context, VerletIntegrator, HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, NonbondedForce
from openmm.unit import kilocalories_per_mole

os.chdir('/tmp')
resid = sys.argv[1]

pdb = PDBFile('%s.pdb' % resid)
ff = ForceField('%s.xml' % resid)
# NoCutoff: an isolated-molecule (gas-phase) evaluation, matching BOSS's
# own NMOL=0 single-point setup -- no periodic images, no truncation.
system = ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff)

GROUP_NAMES = {
    HarmonicBondForce: ('bond', 0),
    HarmonicAngleForce: ('angle', 1),
    PeriodicTorsionForce: ('torsion', 2),
    NonbondedForce: ('nonbonded', 3),
}
term_of_group = {}
for force in system.getForces():
    for cls, (name, group) in GROUP_NAMES.items():
        if isinstance(force, cls):
            force.setForceGroup(group)
            term_of_group[group] = name

integrator = VerletIntegrator(0.001)
context = Context(system, integrator)
context.setPositions(pdb.positions)

total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
print('OPENMM_ENERGY_KCAL_PER_MOL %s' % total)

terms = {}
for group, name in term_of_group.items():
    e = context.getState(getEnergy=True, groups={group}).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    terms[name] = e

print('OPENMM_TERMS bond=%s angle=%s torsion=%s nonbonded=%s' % (
    terms.get('bond'), terms.get('angle'), terms.get('torsion'), terms.get('nonbonded')))
