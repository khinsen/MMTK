# Restraint tests
#
# Written by Konrad Hinsen
#

import unittest
from subsets import SubsetTest
from MMTK import *
import MMTK.ForceFields.Restraints as Restraints
from Scientific import N

class TrapTest(unittest.TestCase):

    def setUp(self):
        self.universe = InfiniteUniverse()
        atom = Atom('C', position=Vector(0.5, 0., 0.))
        self.universe.addObject(atom)
        ff = Restraints.HarmonicTrapForceField(atom, Vector(-0.5, 0., 0.), 1.)
        self.universe.setForceField(ff)

    def test_energy(self):
        e, g, fc = self.universe.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.5, 7)
        self.assert_((g[0]-Vector(1., 0., 0.)).length() < 1.e-5)
        fc = fc[0, 0]
        for i in range(3):
            for j in range(3):
                if i == j:
                    self.assertAlmostEqual(fc[i, j], 1., 7)
                else:
                    self.assertEqual(fc[i,j], 0.)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TrapTest))
    return s

if __name__ == '__main__':
    unittest.main()
