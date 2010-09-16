# Universe tests
#
# Written by Konrad Hinsen
#

import unittest
from MMTK import *
from MMTK.Random import randomPointInBox
from Scientific import N

class PeriodicUniverseTest:

    def test_basisVectors(self):
        e1, e2, e3 = self.universe.basisVectors()
        self.assertEqual(e1, self.a)
        self.assertEqual(e2, self.b)
        self.assertEqual(e3, self.c)

    def test_reciprocalBasisVectors(self):
        r1, r2, r3 = self.universe.reciprocalBasisVectors()
        self.assertAlmostEqual(r1*self.a, 1., 14)
        self.assertAlmostEqual(r2*self.b, 1., 14)
        self.assertAlmostEqual(r3*self.c, 1., 14)
        self.assert_(abs(r1*self.b) < 1.e-15)
        self.assert_(abs(r1*self.c) < 1.e-15)
        self.assert_(abs(r2*self.a) < 1.e-15)
        self.assert_(abs(r2*self.c) < 1.e-15)
        self.assert_(abs(r3*self.a) < 1.e-15)
        self.assert_(abs(r3*self.b) < 1.e-15)

    def test_volume(self):
        volume = self.universe.cellVolume()
        self.assertAlmostEqual(volume, self.a*(self.b.cross(self.c)), 14)

    def test_boxCoordinates1(self):
        for i in range(500):
            p = randomPointInBox(5.)
            pt = self.universe.realToBoxCoordinates(p)
            self.assert_((self.universe.boxToRealCoordinates(pt)-p).length()
                         < 1.e-14)
            pt = self.universe.boxToRealCoordinates(p)
            self.assert_((self.universe.realToBoxCoordinates(pt)-p).length()
                         < 1.e-14)

    def test_boxCoordinates2(self):
        configuration = self.universe.copyConfiguration()
        configuration.convertToBoxCoordinates()
        for atom in self.universe.atomList():
            rb = self.universe.realToBoxCoordinates(atom.position())
            self.assert_((configuration[atom]-rb).length() < 1.e-15)
        configuration.convertFromBoxCoordinates()
        self.assert_(N.maximum.reduce(N.ravel(N.fabs(
            configuration.array-self.universe.configuration().array))) < 1.e-15)

    def test_coordinateFold(self):
        for i in range(500):
            self.universe.translateBy(randomPointInBox(0.5))
            self.universe.foldCoordinatesIntoBox()
            for atom in self.universe.atomList():
                rb = self.universe.realToBoxCoordinates(atom.position())
                self.assert_(rb[0] >= -0.5)
                self.assert_(rb[1] >= -0.5)
                self.assert_(rb[2] >= -0.5)
                self.assert_(rb[0] <= 0.5)
                self.assert_(rb[1] <= 0.5)
                self.assert_(rb[2] <= 0.5)
            cconf = self.universe.contiguousObjectConfiguration(
                                                    self.universe.objectList())
            for o in self.universe:
                for bu in o.bondedUnits():
                    for bond in bu.bonds:
                        l = (bond.a1.position(cconf)
                             - bond.a2.position(cconf)).length()
                        self.assert_(l < 0.16)


class ParallelepipedicPeriodicUniverseTest(unittest.TestCase,
                                           PeriodicUniverseTest):

    def setUp(self):
        self.a = Vector(1.5, 0.3, 0.)
        self.b = Vector(-0.1, 0.8, 0.2)
        self.c = Vector(0., -0.2, 1.1)
        self.universe = ParallelepipedicPeriodicUniverse((self.a,
                                                          self.b,
                                                          self.c))
        self.universe.addObject(Molecule('water'))


class OrthorhombicPeriodicUniverseTest(unittest.TestCase,
                                       PeriodicUniverseTest):

    def setUp(self):
        self.a = Vector(0.7, 0., 0.)
        self.b = Vector(0., 1., 0.)
        self.c = Vector(0., 0., 0.8)
        self.universe = OrthorhombicPeriodicUniverse((self.a.length(),
                                                      self.b.length(),
                                                      self.c.length()))
        self.universe.addObject(Molecule('water'))


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(ParallelepipedicPeriodicUniverseTest))
    s.addTest(loader.loadTestsFromTestCase(OrthorhombicPeriodicUniverseTest))
    return s


if __name__ == '__main__':
    unittest.main()
