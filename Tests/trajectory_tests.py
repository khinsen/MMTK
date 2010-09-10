# Trajectory tests
#
# Written by Konrad Hinsen
# last revision: 2010-9-9
#


import unittest
from MMTK import *
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from Scientific import N
import os

class InfiniteUniverseTest:

    def setUp(self):
        self.universe = InfiniteUniverse()
        self.universe.addObject(Molecule('water',
                                         position = Vector(-0.2, 0., 0.)))
        self.universe.addObject(Molecule('water',
                                         position = Vector(0.2, 0., 0.)))


class OrthorhombicUniverseTest:

    def setUp(self):
        self.universe = OrthorhombicPeriodicUniverse((1., 1., 0.5))
        self.universe.addObject(Molecule('water',
                                         position = Vector(-0.2, 0., 0.)))
        self.universe.addObject(Molecule('water',
                                         position = Vector(0.2, 0., 0.)))

class ParallelepipedicUniverseTest:

    def setUp(self):
        self.universe = ParallelepipedicPeriodicUniverse((Vector(1., 0., 0.),
                                                          Vector(0., 1., 0.),
                                                          Vector(0.5, 0., 0.5)))
        self.universe.addObject(Molecule('water',
                                         position = Vector(-0.2, 0., 0.)))
        self.universe.addObject(Molecule('water',
                                         position = Vector(0.2, 0., 0.)))

class InfiniteUniverseWithPITest:

    def setUp(self):
        self.universe = InfiniteUniverse()
        self.universe.addObject(Molecule('water',
                                         position = Vector(-0.2, 0., 0.)))
        for a in self.universe.atomList():
            a.setNBeads(4)

class SinglePrecisionTest:

    double_precision = False
    tolerance = 1.e-7


class DoublePrecisionTest:

    double_precision = True
    tolerance = 1.e-15


class TrajectoryTest:

    def tearDown(self):
        return
        try:
            os.remove('test.nc')
        except OSError:
            pass

    def test_snapshot(self):

        initial = self.universe.copyConfiguration()

        transformation = Translation(Vector(0.,0.,0.01)) \
                         * Rotation(Vector(0.,0.,1.), 1.*Units.deg)

        trajectory = Trajectory(self.universe, "test.nc", "w",
                                "trajectory test",
                                double_precision = self.double_precision)
        snapshot = SnapshotGenerator(self.universe,
                                     actions = [TrajectoryOutput(trajectory,
                                                                 ["all"],
                                                                 0, None, 1)])
        snapshot()
        for i in range(100):
            self.universe.setConfiguration(
                self.universe.contiguousObjectConfiguration())
            self.universe.applyTransformation(transformation)
            self.universe.foldCoordinatesIntoBox()
            snapshot()
        trajectory.close()

        self.universe.setConfiguration(initial)
        trajectory = Trajectory(None, "test.nc")
        t_universe = trajectory.universe
        for i in range(101):
            configuration = self.universe.configuration()
            t_configuration = trajectory[i]['configuration']
            max_diff = N.maximum.reduce(N.ravel(N.fabs(configuration.array
                                                     - t_configuration.array)))
            self.assert_(max_diff < self.tolerance)
            if configuration.cell_parameters is not None:
                max_diff = N.maximum.reduce(N.fabs(configuration.cell_parameters
                                            - t_configuration.cell_parameters))
                self.assert_(max_diff < self.tolerance)
            self.universe.setConfiguration(
                self.universe.contiguousObjectConfiguration())
            self.universe.applyTransformation(transformation)
            self.universe.foldCoordinatesIntoBox()
        trajectory.close()


class InfiniteUniverseTestSP(unittest.TestCase,
                             InfiniteUniverseTest,
                             TrajectoryTest,
                             SinglePrecisionTest):
    setUp = InfiniteUniverseTest.setUp

class InfiniteUniverseTestDP(unittest.TestCase,
                             InfiniteUniverseTest,
                             TrajectoryTest,
                             DoublePrecisionTest):
    setUp = InfiniteUniverseTest.setUp
    tearDown = TrajectoryTest.tearDown

class OrthorhombicUniverseTestSP(unittest.TestCase,
                                 OrthorhombicUniverseTest,
                                 TrajectoryTest,
                                 SinglePrecisionTest):
    setUp = OrthorhombicUniverseTest.setUp
    tearDown = TrajectoryTest.tearDown

class OrthorhombicUniverseTestDP(unittest.TestCase,
                                 OrthorhombicUniverseTest,
                                 TrajectoryTest,
                                 DoublePrecisionTest):
    setUp = OrthorhombicUniverseTest.setUp
    tearDown = TrajectoryTest.tearDown

class ParallelepipedicUniverseTestSP(unittest.TestCase,
                                     ParallelepipedicUniverseTest,
                                     TrajectoryTest,
                                     SinglePrecisionTest):
    setUp = ParallelepipedicUniverseTest.setUp
    tearDown = TrajectoryTest.tearDown

class ParallelepipedicUniverseTestDP(unittest.TestCase,
                                     ParallelepipedicUniverseTest,
                                     TrajectoryTest,
                                     DoublePrecisionTest):
    setUp = ParallelepipedicUniverseTest.setUp
    tearDown = TrajectoryTest.tearDown

class InfiniteUniverseWithPITestDP(unittest.TestCase,
                                   InfiniteUniverseWithPITest,
                                   TrajectoryTest,
                                   DoublePrecisionTest):
    setUp = InfiniteUniverseWithPITest.setUp
    tearDown = TrajectoryTest.tearDown

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(InfiniteUniverseTestSP))
    s.addTest(loader.loadTestsFromTestCase(InfiniteUniverseTestDP))
    s.addTest(loader.loadTestsFromTestCase(OrthorhombicUniverseTestSP))
    s.addTest(loader.loadTestsFromTestCase(OrthorhombicUniverseTestDP))
    s.addTest(loader.loadTestsFromTestCase(ParallelepipedicUniverseTestSP))
    s.addTest(loader.loadTestsFromTestCase(ParallelepipedicUniverseTestDP))
    s.addTest(loader.loadTestsFromTestCase(InfiniteUniverseWithPITestDP))
    return s


if __name__ == '__main__':
    unittest.main()
