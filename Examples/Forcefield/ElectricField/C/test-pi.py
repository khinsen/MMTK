from MMTK import *
from ElectricField import ElectricField
from MMTK.ForceFields.ForceFieldTest import gradientTest
from MMTK.Environment import PathIntegrals

universe = InfiniteUniverse()
universe.atom1 = Atom('C', nbeads=4, position=Vector(0., 0., 1.))
universe.atom1.test_charge = 1.
universe.atom2 = Atom('C', nbeads=4, position=Vector(0., 0., 0.))
universe.atom2.test_charge = -0.2

universe.pathintegrals = PathIntegrals(100*Units.K)


universe.setForceField(ElectricField(0.5*(Units.V/Units.m)
                                     * Vector(0., 0., 1.),
                                     'test_charge'))

print universe.energyTerms()
e, g = universe.energyAndGradients()
print g[universe.atom1.beads()[0]]
print g[universe.atom2.beads()[0]]

gradientTest(universe)
