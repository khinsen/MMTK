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

        self.universe1 = InfiniteUniverse()
        atom = Atom('C', position=Vector(0.5, 0., 0.))
        self.universe1.addObject(atom)
        ff = Restraints.HarmonicTrapForceField(atom,
                                               Vector(-0.5, 0., 0.),
                                               1.)
        self.universe1.setForceField(ff)

        self.universe2 = InfiniteUniverse()
        atom1 = Atom('C', position=Vector(0.25, 0., 0.))
        atom2 = Atom('C', position=Vector(0.75, 0., 0.))
        cluster = Collection([atom1, atom2])
        self.universe2.addObject(cluster)
        ff = Restraints.HarmonicTrapForceField(cluster,
                                               Vector(-0.5, 0., 0.),
                                               1.)
        self.universe2.setForceField(ff)

        self.universe3 = InfiniteUniverse()
        water = Molecule('water', position=Vector(0.5, 0., 0.))
        self.universe3.addObject(water)
        ff = Restraints.HarmonicTrapForceField(water,
                                               Vector(-0.5, 0., 0.),
                                               1.)
        self.universe3.setForceField(ff)

    def test_energy_one_atom(self):
        e, g, fc = self.universe1.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.5, 7)
        self.assert_((g[0]-Vector(1., 0., 0.)).length() < 1.e-5)
        fc = fc[0, 0]
        for i in range(3):
            for j in range(3):
                if i == j:
                    self.assertAlmostEqual(fc[i, j], 1., 7)
                else:
                    self.assertEqual(fc[i,j], 0.)

    def test_energy_two_atoms(self):
        e, g, fc = self.universe2.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.5, 7)
        self.assert_((g[0]-Vector(0.5, 0., 0.)).length() < 1.e-5)
        self.assert_((g[1]-Vector(0.5, 0., 0.)).length() < 1.e-5)
        for a1 in [0, 1]:
            for a2 in [0, 1]:
                fc12 = fc[a1, a2]
                for i in range(3):
                    for j in range(3):
                        if i == j:
                            self.assertAlmostEqual(fc12[i, j], 0.25, 7)
                        else:
                            self.assertEqual(fc12[i,j], 0.)

    def test_energy_water(self):

        H1 = self.universe3[0].H1
        H2 = self.universe3[0].H2
        O = self.universe3[0].O
        e, g, fc = self.universe3.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.5, 7)
        self.assert_((g[H1]-g[H2]).length() < 1.e-5)
        self.assert_((g[H1]-g[O]*H1.mass()/O.mass()).length() < 1.e-5)

        m = 2*H1.mass() + O.mass()
        for a1 in [H1, H2, O]:
            for a2 in [H1, H2, O]:
                fc12 = fc[a1, a2]
                f = a1.mass()*a2.mass()/(m*m)
                for i in range(3):
                    for j in range(3):
                        if i == j:
                            self.assertAlmostEqual(fc12[i, j], f, 7)
                        else:
                            self.assertEqual(fc12[i,j], 0.)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TrapTest))
    return s

if __name__ == '__main__':
    unittest.main()
