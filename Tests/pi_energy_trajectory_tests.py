# Path integral energy tests
#
# Written by Chris Ing
#

import unittest
from MMTK import *
import MMTK.Utility
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from MMTK.PINormalModeIntegrator import PINormalModeIntegrator, PILangevinNormalModeIntegrator
from MMTK.ForceFields import Amber99ForceField
from Scientific import N
import os

def springEnergyCartesian(self):
    beta = self.beta
    e = 0.
    for a in self.universe.atomIterator():
        nb = a.numberOfBeads()
        x = N.array(a.beadPositions())
        sumsq = 0.
        for j in range(nb):
            for z in range(3):
                sumsq += (x[(j+1)%nb, z]-x[j, z])**2
        e += 0.5*nb*a.mass()*sumsq/(beta*beta*Units.hbar*Units.hbar)
    return e

def virialTerm(self):
    beadsH1 = self.universe.atomList()[0].beadPositions()
    beadsH2 = self.universe.atomList()[1].beadPositions()
    nb1 = len(beadsH1)
    nb2 = len(beadsH2)
    r_diff = N.array(beadsH1)-N.array(beadsH2)
    r_grad = self.universe.energyAndGradients()[1].array[nb1:nb1+nb2]
    return N.dot((r_diff).flat, (r_grad).flat)

class AmberInfiniteTest(unittest.TestCase):

    def test_variable_bead_atoms(self):
        temperature = 100*Units.K
        self.universe = InfiniteUniverse(Amber99ForceField())
        self.universe.addObject(Atom('H', nbeads=2, amber_atom_type="CT",amber_charge=0.))
        self.universe.addObject(Atom('H', nbeads=2, amber_atom_type="CT",amber_charge=0.))
        self.universe.atomList()[0].setBeadPositions([MMTK.Vector(-1.0,0.,0.),MMTK.Vector(-1.0,0.,0.)])
        self.universe.atomList()[1].setBeadPositions([MMTK.Vector(2.5,0.,0.),MMTK.Vector(1.5,0.,0.)])
        self.universe.addObject(PathIntegrals(temperature, False))
        e1 =  self.universe.energy()

        self.universe.atomList()[1].setNumberOfBeads(4)
        self.universe.atomList()[1].setBeadPositions([MMTK.Vector(2.5,0.,0.),MMTK.Vector(2.5,0.,0.),
                                                 MMTK.Vector(1.5,0.,0.),MMTK.Vector(1.5,0.,0.)])
        self.assert_(e1 - self.universe.energy() < 1e-5)

    def test_estimators(self):
        temperature = 100*Units.K
        self.universe = InfiniteUniverse(Amber99ForceField())
        self.universe.addObject(Atom('H', nbeads=4, position= Vector(-0.2, 0., 0.),
                                     amber_atom_type="CT",amber_charge=0.))
        self.universe.addObject(Atom('H', nbeads=4, position= Vector(0.2, 0., 0.),
                                     amber_atom_type="CT",amber_charge=0.))
        self.universe.addObject(PathIntegrals(temperature, False))
        self.beta=self.universe.environmentObjectList(Environment.PathIntegrals)[0].beta
        self.np = self.universe.numberOfPoints()

        self.universe.initializeVelocitiesToTemperature(temperature)

        integrator = PINormalModeIntegrator(self.universe, delta_t=1.*Units.fs)
        trajectory = Trajectory(self.universe, "test.nc", "w", "Test Trajectory")

        integrator(steps = 100,
           actions = [TrajectoryOutput(trajectory, ("time", "energy", 
                                                    "configuration", "auxiliary"), 0, None, 1)])
        skip = 0
        for step in trajectory:
            if skip > 0:
                self.universe.setConfiguration(step['configuration'])
                Vcl = self.universe.energy()
                Vqu = springEnergyCartesian(self)
                self.assertAlmostEqual(Vcl - Vqu + 1.5*self.np/self.beta, step['quantum_energy_primitive'], 5)
                self.assertAlmostEqual(Vcl - 0.5*virialTerm(self), step['quantum_energy_virial'], 5)
            skip += 1

        trajectory.close()
        try:
            os.remove('test.nc')
        except OSError:
            pass
        
def suite():
    return unittest.TestLoader().loadTestsFromTestCase(AmberInfiniteTest)

if __name__ == '__main__':
    unittest.main()
