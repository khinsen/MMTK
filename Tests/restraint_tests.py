# Restraint tests
#
# Written by Konrad Hinsen
#

import unittest
from subsets import SubsetTest
from MMTK import *
from MMTK.ForceFields import Restraints, SPCEForceField
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
        self.assertAlmostEqual(e, 1., 7)
        self.assert_((g[0]-Vector(2., 0., 0.)).length() < 1.e-5)
        fc = fc[0, 0]
        for i in range(3):
            for j in range(3):
                if i == j:
                    self.assertAlmostEqual(fc[i, j], 2., 7)
                else:
                    self.assertEqual(fc[i,j], 0.)

    def test_energy_two_atoms(self):
        e, g, fc = self.universe2.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 1., 7)
        self.assert_((g[0]-Vector(1., 0., 0.)).length() < 1.e-5)
        self.assert_((g[1]-Vector(1., 0., 0.)).length() < 1.e-5)
        for a1 in [0, 1]:
            for a2 in [0, 1]:
                fc12 = fc[a1, a2]
                for i in range(3):
                    for j in range(3):
                        if i == j:
                            self.assertAlmostEqual(fc12[i, j], 0.5, 7)
                        else:
                            self.assertEqual(fc12[i,j], 0.)

    def test_energy_water(self):

        H1 = self.universe3[0].H1
        H2 = self.universe3[0].H2
        O = self.universe3[0].O
        e, g, fc = self.universe3.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 1., 7)
        self.assert_((g[H1]-g[H2]).length() < 1.e-5)
        self.assert_((g[H1]-g[O]*H1.mass()/O.mass()).length() < 1.e-5)

        m = 2*H1.mass() + O.mass()
        for a1 in [H1, H2, O]:
            for a2 in [H1, H2, O]:
                fc12 = fc[a1, a2]
                f = 2.*a1.mass()*a2.mass()/(m*m)
                for i in range(3):
                    for j in range(3):
                        if i == j:
                            self.assertAlmostEqual(fc12[i, j], f, 7)
                        else:
                            self.assertEqual(fc12[i,j], 0.)

    def test_subsets(self):
        H1 = self.universe3[0].H1
        H2 = self.universe3[0].H2
        O = self.universe3[0].O
        subset = Collection([H1, H2, O])
        e = self.universe3.energy(subset)
        self.assertAlmostEqual(e, 1., 7)
        subset = Collection()
        e = self.universe3.energy(subset)
        self.assertAlmostEqual(e, 0., 7)
        subset = Collection([H1, H2])
        self.assertRaises(ValueError, self.universe3.energy, subset)
        

class DistanceTest(unittest.TestCase):

    def setUp(self):

        self.universe1 = InfiniteUniverse()
        atom1 = Atom('C', position=Vector(0.5, 0., 0.))
        atom2 = Atom('C', position=Vector(-0.5, 0., 0.))
        self.universe1.addObject(atom1)
        self.universe1.addObject(atom2)
        ff = Restraints.HarmonicDistanceRestraint(atom1, atom2,
                                                   0.9, 1.)
        self.universe1.setForceField(ff)

        self.universe2 = InfiniteUniverse()
        atom1 = Atom('C', position=Vector(0.25, 0., 0.))
        atom2 = Atom('C', position=Vector(0.75, 0., 0.))
        cluster1 = Collection([atom1, atom2])
        atom3 = Atom('C', position=Vector(-0.75, 0., 0.))
        atom4 = Atom('C', position=Vector(-0.25, 0., 0.))
        cluster2 = Collection([atom3, atom4])
        self.universe2.addObject(cluster1)
        self.universe2.addObject(cluster2)
        ff = Restraints.HarmonicDistanceRestraint(cluster1, cluster2,
                                                  0.9, 1.)
        self.universe2.setForceField(ff)

        self.universe3 = InfiniteUniverse()
        water1 = Molecule('water', position=Vector(0.5, 0., 0.))
        water2 = Molecule('water', position=Vector(-0.5, 0., 0.))
        self.universe3.addObject(water1)
        self.universe3.addObject(water2)
        ff = Restraints.HarmonicDistanceRestraint(water1, water2,
                                                  0.9, 1.)
        self.universe3.setForceField(ff)

    def test_energy_single_atoms(self):
        e, g, fc = self.universe1.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.01, 7)
        self.assert_((g[0]-Vector(0.2, 0., 0.)).length() < 1.e-5)
        fc = fc[0, 0]
        fc_ref = [2., 0.2, 0.2]
        for i in range(3):
            self.assertAlmostEqual(fc[i, i], fc_ref[i], 7)
            for j in range(3):
                if i != j:
                    self.assertEqual(fc[i,j], 0.)

    def test_energy_atom_pairs(self):
        e, g, fc = self.universe2.energyGradientsAndForceConstants()
        self.assertAlmostEqual(e, 0.01, 7)
        self.assert_((g[0]-Vector(0.1, 0., 0.)).length() < 1.e-5)
        self.assert_((g[1]-Vector(0.1, 0., 0.)).length() < 1.e-5)
        self.assert_((g[2]-Vector(-0.1, 0., 0.)).length() < 1.e-5)
        self.assert_((g[3]-Vector(-0.1, 0., 0.)).length() < 1.e-5)
        p1 = [0, 1]
        p2 = [2, 3]
        fc_ref = [0.5, 0.05, 0.05]
        for a1 in range(4):
            for a2 in range(4):
                fcp = fc[a1, a2]
                for i in range(3):
                    if (a1 in p1 and a2 in p1) or (a1 in p2 and a2 in p2):
                        self.assertAlmostEqual(fcp[i, i], fc_ref[i], 7)
                    else:
                        self.assertAlmostEqual(fcp[i, i], -fc_ref[i], 7)
                    for j in range(3):
                        if i != j:
                            self.assertEqual(fcp[i,j], 0.)

    def test_energy_water(self):
        e, g, fc = self.universe3.energyGradientsAndForceConstants()
        fc_ref = [2., 0.2, 0.2]
        self.assertAlmostEqual(e, 0.01, 7)
        for m, s in [(self.universe3[0], 1), (self.universe3[1], -1),]:
            self.assert_((g[m.H1]-g[m.H2]).length() < 1.e-5)
            self.assert_((g[m.H1]-g[m.O]*m.H1.mass()/m.O.mass()).length()
                         < 1.e-5)
            mass = 2*m.H1.mass() + m.O.mass()
            for a1 in [m.H1, m.H2, m.O]:
                for a2 in [m.H1, m.H2, m.O]:
                    fc12 = fc[a1, a2]
                    f = a1.mass()*a2.mass()/(mass*mass)
                    for i in range(3):
                        for j in range(3):
                            if i == j:
                                self.assertAlmostEqual(fc12[i, j], f*fc_ref[i], 7)
                            else:
                                self.assertEqual(fc12[i,j], 0.)

    def test_subsets(self):
        m1 = self.universe3[0]
        m2 = self.universe3[1]
        subset = Collection([m1.H1, m1.H2, m1.O, m2.H1, m2.H2, m2.O])
        e = self.universe3.energy(subset)
        self.assertAlmostEqual(e, 0.01, 7)
        subset = Collection()
        e = self.universe3.energy(subset)
        self.assertAlmostEqual(e, 0., 7)
        subset = Collection([m1.H1, m2.H2])
        self.assertRaises(ValueError, self.universe3.energy, subset)

class NonbondedExclusionTest(unittest.TestCase):

    def test_exclusion(self):
        universe = InfiniteUniverse()
        w1 = Molecule('water', position=(Vector(-0.2, 0., 0.)))
        w2 = Molecule('water', position=(Vector(+0.2, 0., 0.)))
        universe.addObject(w1)
        universe.addObject(w2)
        ff = SPCEForceField() + \
               Restraints.HarmonicDistanceRestraint(w1.O, w2.O, 0.35, 1., False)
        universe.setForceField(ff)
        self.assert_(universe.energyTerms()['Lennard-Jones'] < 0.)
        ff = SPCEForceField() + \
               Restraints.HarmonicDistanceRestraint(w1.O, w2.O, 0.35, 1., True)
        universe.setForceField(ff)
        self.assert_(universe.energyTerms()['Lennard-Jones'] == 0.)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(TrapTest))
    s.addTest(loader.loadTestsFromTestCase(DistanceTest))
    s.addTest(loader.loadTestsFromTestCase(NonbondedExclusionTest))
    return s

if __name__ == '__main__':
    unittest.main()
