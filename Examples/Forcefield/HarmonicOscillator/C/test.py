from MMTK import *
from HarmonicOscillatorFF import HarmonicOscillatorForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest, virialTest

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector(0., 0., 1.05), name='C1')
universe.atom2 = Atom('C', position=Vector(0., 1.05, 0.), name='C2')

ff1 = HarmonicOscillatorForceField(universe.atom1, Vector(0., 0., 1.), 100.)
ff2 = HarmonicOscillatorForceField(universe.atom2, Vector(0., 1., 0.), 100.)
universe.setForceField(ff1+ff2)

e, g = universe.energyAndGradients()
print universe.energyTerms()
print e

gradientTest(universe)
forceConstantTest(universe)
virialTest(universe)
