from MMTK import *
from HarmonicOscillatorFF import HarmonicOscillatorForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest
from MMTK.Environment import PathIntegrals
from MMTK.Dynamics import VelocityVerletIntegrator

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector(0., 0., 1.05), name='C1', nbeads=2)
universe.atom2 = Atom('C', position=Vector(0., 1.05, 0.), name='C2', nbeads=2)
universe.path_integrals = PathIntegrals(50.*Units.K, include_spring_terms=True)

ff1 = HarmonicOscillatorForceField(universe.atom1, Vector(0., 0., 1.), 1000.)
ff2 = HarmonicOscillatorForceField(universe.atom2, Vector(0., 1., 0.), 1000.)
universe.setForceField(ff1+ff2)

# We need to move the beads apart because in the (unrealistic)
# initial configuration where all beads of an atom are at the same
# place, the numerical derivatives of the path integral spring terms
# are unstable.
universe.initializeVelocitiesToTemperature(50.*Units.K)
VelocityVerletIntegrator(universe)(steps=50)

e, g = universe.energyAndGradients()
print universe.energyTerms()
print e

gradientTest(universe)
forceConstantTest(universe)
