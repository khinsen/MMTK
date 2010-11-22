# Subspace tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
from MMTK.Proteins import Protein
from MMTK.ForceFields import HarmonicForceField
from MMTK.Subspace import RigidMotionSubspace, PairDistanceSubspace
from MMTK import NormalModes
from Scientific.Geometry.Transformation import Rotation, Translation

class RigidBodyTest(unittest.TestCase):

    """
    Test rigid-body motion subspace
    """
    
    def test_projections(self):
        universe = MMTK.InfiniteUniverse()
        universe.m1 = MMTK.Molecule('water',
                                         position=MMTK.Vector(0., 0., 0.))
        universe.m2 = MMTK.Molecule('water',
                                         position=MMTK.Vector(1., 0., 0.))
        rbs = [universe.m1,
               MMTK.Collection([universe.m2.H1,
                                universe.m2.H2])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s), 11)
        basis = s.getBasis()
        self.assertEqual(len(basis), 11)
        complement = s.complement()
        complement_basis = complement.getBasis()
        self.assertEqual(len(complement_basis),
                         universe.degreesOfFreedom()-11)
        for rb in rbs:
            for t in [Translation(MMTK.Vector(0.1, -0.2, 1.5)),
                      Rotation(MMTK.Vector(0.3, 1.2, -2.3), 0.001)]:
                d = rb.displacementUnderTransformation(t)
                self.assert_((s.projectionOf(d)-d).norm() < 1.e-7)
                self.assert_(s.projectionComplementOf(d).norm() < 1.e-7)

    def test_lrb(self):
        universe = MMTK.InfiniteUniverse()
        for p in [MMTK.Vector(0., 0., 0.),
                  MMTK.Vector(1., 0., 0.),
                  MMTK.Vector(0., 1., 1.),
                  MMTK.Vector(0., 1., 0.),
                  MMTK.Vector(1., 0., 1.),
                  MMTK.Vector(0., 0., 1.),
                  MMTK.Vector(1., 1., 0.),
                  MMTK.Vector(1., 1., 1.)]:
            universe.addObject(MMTK.Atom('C', position=p))
        atoms = universe.atomList()
        # all atoms independent
        s = RigidMotionSubspace(universe, atoms)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 3*len(atoms))
        # 1 rb, four free atoms
        rbs = [MMTK.Collection(atoms[:4])] + atoms[4:]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 6+4*3)
        # 2 independent rbs with > 2 atoms
        rbs = [MMTK.Collection(atoms[:4]),
               MMTK.Collection(atoms[4:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 12)
        # 2 independent rbs, one 2 atoms, one 6 atoms
        rbs = [MMTK.Collection(atoms[:6]),
               MMTK.Collection(atoms[6:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 11)
        # 2 rbs > 2 atoms each, one atom in common
        rbs = [MMTK.Collection(atoms[:5]),
               MMTK.Collection(atoms[4:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 9)
        # 2 rbs > 2 atoms each, two atoms in common
        rbs = [MMTK.Collection(atoms[:5]),
               MMTK.Collection(atoms[3:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 7)
        # 2 rbs > 2 atoms each, three atoms in common
        rbs = [MMTK.Collection(atoms[:5]),
               MMTK.Collection(atoms[2:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 6)
        # 2 rbs with 2 and 7 atoms, one atom in common
        rbs = [MMTK.Collection(atoms[:2]),
               MMTK.Collection(atoms[1:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 8)
        # 3 rbs > 2 atoms each, one atom in common between rb[n] and rb[n+1]
        rbs = [MMTK.Collection(atoms[:3]),
               MMTK.Collection(atoms[2:6]),
               MMTK.Collection(atoms[5:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 12)
        # 3 rbs > 2 atoms each, one atom in common for each pair
        # (chain of 3 rbs)
        rbs = [MMTK.Collection(atoms[:3]+atoms[7:]),
               MMTK.Collection(atoms[2:6]),
               MMTK.Collection(atoms[5:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 10)
        # 2 rbs > 2 atoms each, linked by a rigid bond
        rbs = [MMTK.Collection(atoms[:4]),
               MMTK.Collection(atoms[3:5]),
               MMTK.Collection(atoms[4:])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 11)

    def test_projections_pi(self):
        universe = MMTK.InfiniteUniverse()
        universe.m1 = MMTK.Molecule('water',
                                         position=MMTK.Vector(0., 0., 0.))
        universe.m2 = MMTK.Molecule('water',
                                         position=MMTK.Vector(1., 0., 0.))
        nb = 4
        for atom in universe.atomIterator():
            atom.setNumberOfBeads(nb)
        rbs = [universe.m1,
               MMTK.Collection([universe.m2.H1,
                                universe.m2.H2])]
        s = RigidMotionSubspace(universe, rbs)
        self.checkOrthonormality(s)
        self.assertEqual(len(s), nb*11)
        basis = s.getBasis()
        self.assertEqual(len(basis), nb*11)
        complement = s.complement()
        complement_basis = complement.getBasis()
        self.assertEqual(len(complement_basis),
                         universe.degreesOfFreedom()-nb*11)
        for rb in rbs:
            for t in [Translation(MMTK.Vector(0.1, -0.2, 1.5)),
                      Rotation(MMTK.Vector(0.3, 1.2, -2.3), 0.0001)]:
                d = rb.displacementUnderTransformation(t)
                self.assert_((s.projectionOf(d)-d).norm() < 1.e-7)
                self.assert_(s.projectionComplementOf(d).norm() < 1.e-7)

    def test_nbeads_1(self):
        universe = MMTK.InfiniteUniverse()
        universe.m = MMTK.Molecule('water',
                                   position=MMTK.Vector(0., 0., 0.))
        for a in universe.m.atomIterator():
            a.setNumberOfBeads(2)
        s = RigidMotionSubspace(universe, [universe.m])
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 12)

    def test_nbeads_2(self):
        universe = MMTK.InfiniteUniverse()
        universe.m = MMTK.Molecule('water',
                                   position=MMTK.Vector(0., 0., 0.))
        for a in [universe.m.H1, universe.m.H2]:
            a.setNumberOfBeads(2)
        s = RigidMotionSubspace(universe, [universe.m])
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 9)

    def test_nbeads_3(self):
        universe = MMTK.InfiniteUniverse()
        universe.m = MMTK.Molecule('water',
                                   position=MMTK.Vector(0., 0., 0.))
        universe.m.O.setNumberOfBeads(2)
        s = RigidMotionSubspace(universe, [universe.m])
        self.checkOrthonormality(s)
        self.assertEqual(len(s.getBasis()), 7)

    def checkOrthonormality(self, subspace):
        basis = subspace.getBasis()
        ndim = len(basis)
        for i in range(ndim):
            v = basis[i]
            self.assertAlmostEqual(v.dotProduct(v), 1., 10)
            for j in range(i+1, ndim):
                self.assert_(abs(v.dotProduct(basis[j])) < 1.e-10)


class PeptideNormalModeTest(unittest.TestCase):

    """
    Test mode projections on a subspace
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(HarmonicForceField())
        self.universe.peptide = Protein('bala1')
        self.emodes = NormalModes.EnergeticModes(self.universe)
        self.rm = RigidMotionSubspace(self.universe,
                                      self.universe.peptide.residues())
        self.pd = PairDistanceSubspace(self.universe,
                                       [(self.universe.peptide[0][0].O,
                                         self.universe.peptide[0][-1].CH3)])

    def test_rmProjections(self):
        p = self.rm.projectionOf(self.emodes.rawMode(6)).norm()
        self.assertAlmostEqual(p, 0.932906130224)
        p = self.rm.projectionOf(self.emodes.rawMode(7)).norm()
        self.assertAlmostEqual(p, 0.946620894164)
        p = self.rm.projectionOf(self.emodes.rawMode(8)).norm()
        self.assertAlmostEqual(p, 0.918711763134)
        p = self.rm.projectionOf(self.emodes.rawMode(9)).norm()
        self.assertAlmostEqual(p, 0.914888590566)
        p = self.rm.projectionOf(self.emodes.rawMode(10)).norm()
        self.assertAlmostEqual(p, 0.954987614764)
        p = self.rm.projectionOf(self.emodes.rawMode(11)).norm()
        self.assertAlmostEqual(p, 0.490550926574)
        p = self.rm.projectionOf(self.emodes.rawMode(12)).norm()
        self.assertAlmostEqual(p, 0.974543149457)
        p = self.rm.projectionOf(self.emodes.rawMode(13)).norm()
        self.assertAlmostEqual(p, 0.902492405825)
        p = self.rm.projectionOf(self.emodes.rawMode(14)).norm()
        self.assertAlmostEqual(p, 0.873380594746)
        p = self.rm.projectionOf(self.emodes.rawMode(15)).norm()
        self.assertAlmostEqual(p, 0.772133237884)
        p = self.rm.projectionOf(self.emodes.rawMode(16)).norm()
        self.assertAlmostEqual(p, 0.883756751312)
        p = self.rm.projectionOf(self.emodes.rawMode(17)).norm()
        self.assertAlmostEqual(p, 0.577569822435)
        p = self.rm.projectionOf(self.emodes.rawMode(18)).norm()
        self.assertAlmostEqual(p, 0.703601140498)
        p = self.rm.projectionOf(self.emodes.rawMode(19)).norm()
        self.assertAlmostEqual(p, 0.341866735125)
        p = self.rm.projectionOf(self.emodes.rawMode(20)).norm()
        self.assertAlmostEqual(p, 0.590049489902)
        p = self.rm.projectionOf(self.emodes.rawMode(21)).norm()
        self.assertAlmostEqual(p, 0.432251617574)
        p = self.rm.projectionOf(self.emodes.rawMode(22)).norm()
        self.assertAlmostEqual(p, 0.245304052688)
        p = self.rm.projectionOf(self.emodes.rawMode(23)).norm()
        self.assertAlmostEqual(p, 0.211381043933)
        p = self.rm.projectionOf(self.emodes.rawMode(24)).norm()
        self.assertAlmostEqual(p, 0.401154809528)
        p = self.rm.projectionOf(self.emodes.rawMode(25)).norm()
        self.assertAlmostEqual(p, 0.161820108479)
        p = self.rm.projectionOf(self.emodes.rawMode(26)).norm()
        self.assertAlmostEqual(p, 0.208594682213)
        p = self.rm.projectionOf(self.emodes.rawMode(27)).norm()
        self.assertAlmostEqual(p, 0.329444853782)
        p = self.rm.projectionOf(self.emodes.rawMode(28)).norm()
        self.assertAlmostEqual(p, 0.335842626795)
        p = self.rm.projectionOf(self.emodes.rawMode(29)).norm()
        self.assertAlmostEqual(p, 0.246955423968)
        p = self.rm.projectionOf(self.emodes.rawMode(30)).norm()
        self.assertAlmostEqual(p, 0.102801852182)
        p = self.rm.projectionOf(self.emodes.rawMode(31)).norm()
        self.assertAlmostEqual(p, 0.0445276342964)
        p = self.rm.projectionOf(self.emodes.rawMode(32)).norm()
        self.assertAlmostEqual(p, 0.0810509154227)
        p = self.rm.projectionOf(self.emodes.rawMode(33)).norm()
        self.assertAlmostEqual(p, 0.166723733735)
        p = self.rm.projectionOf(self.emodes.rawMode(34)).norm()
        self.assertAlmostEqual(p, 0.254081370404)
        p = self.rm.projectionOf(self.emodes.rawMode(35)).norm()
        self.assertAlmostEqual(p, 0.221723386921)
        p = self.rm.projectionOf(self.emodes.rawMode(36)).norm()
        self.assertAlmostEqual(p, 0.133163531911)
        p = self.rm.projectionOf(self.emodes.rawMode(37)).norm()
        self.assertAlmostEqual(p, 0.0662977755488)
        p = self.rm.projectionOf(self.emodes.rawMode(38)).norm()
        self.assertAlmostEqual(p, 0.156961937663)
        p = self.rm.projectionOf(self.emodes.rawMode(39)).norm()
        self.assertAlmostEqual(p, 0.310640400206)
        p = self.rm.projectionOf(self.emodes.rawMode(40)).norm()
        self.assertAlmostEqual(p, 0.424586976986)
        p = self.rm.projectionOf(self.emodes.rawMode(41)).norm()
        self.assertAlmostEqual(p, 0.439565895727)
        p = self.rm.projectionOf(self.emodes.rawMode(42)).norm()
        self.assertAlmostEqual(p, 0.220029858675)
        p = self.rm.projectionOf(self.emodes.rawMode(43)).norm()
        self.assertAlmostEqual(p, 0.26737834485)
        p = self.rm.projectionOf(self.emodes.rawMode(44)).norm()
        self.assertAlmostEqual(p, 0.367708993752)
        p = self.rm.projectionOf(self.emodes.rawMode(45)).norm()
        self.assertAlmostEqual(p, 0.40324227143)
        p = self.rm.projectionOf(self.emodes.rawMode(46)).norm()
        self.assertAlmostEqual(p, 0.192735085639)
        p = self.rm.projectionOf(self.emodes.rawMode(47)).norm()
        self.assertAlmostEqual(p, 0.0840421758091)
        p = self.rm.projectionOf(self.emodes.rawMode(48)).norm()
        self.assertAlmostEqual(p, 0.15035662682)
        p = self.rm.projectionOf(self.emodes.rawMode(49)).norm()
        self.assertAlmostEqual(p, 0.154260578723)
        p = self.rm.projectionOf(self.emodes.rawMode(50)).norm()
        self.assertAlmostEqual(p, 0.125542420544)
        p = self.rm.projectionOf(self.emodes.rawMode(51)).norm()
        self.assertAlmostEqual(p, 0.0451196639208)
        p = self.rm.projectionOf(self.emodes.rawMode(52)).norm()
        self.assertAlmostEqual(p, 0.0458762025086)
        p = self.rm.projectionOf(self.emodes.rawMode(53)).norm()
        self.assertAlmostEqual(p, 0.00887400641221)
        p = self.rm.projectionOf(self.emodes.rawMode(54)).norm()
        self.assertAlmostEqual(p, 0.017086967116)
        p = self.rm.projectionOf(self.emodes.rawMode(55)).norm()
        self.assertAlmostEqual(p, 0.0609002158356)
        p = self.rm.projectionOf(self.emodes.rawMode(56)).norm()
        self.assertAlmostEqual(p, 0.0343388857234)
        p = self.rm.projectionOf(self.emodes.rawMode(57)).norm()
        self.assertAlmostEqual(p, 0.0559503540031)
        p = self.rm.projectionOf(self.emodes.rawMode(58)).norm()
        self.assertAlmostEqual(p, 0.0332620617319)
        p = self.rm.projectionOf(self.emodes.rawMode(59)).norm()
        self.assertAlmostEqual(p, 0.0583023511826)
        p = self.rm.projectionOf(self.emodes.rawMode(60)).norm()
        self.assertAlmostEqual(p, 0.0663342593349)
        p = self.rm.projectionOf(self.emodes.rawMode(61)).norm()
        self.assertAlmostEqual(p, 0.0568833304731)
        p = self.rm.projectionOf(self.emodes.rawMode(62)).norm()
        self.assertAlmostEqual(p, 0.161627264195)
        p = self.rm.projectionOf(self.emodes.rawMode(63)).norm()
        self.assertAlmostEqual(p, 0.1722951381)
        p = self.rm.projectionOf(self.emodes.rawMode(64)).norm()
        self.assertAlmostEqual(p, 0.239075225836)
        p = self.rm.projectionOf(self.emodes.rawMode(65)).norm()
        self.assertAlmostEqual(p, 0.241338995808)

    def test_pdProjections(self):
        p = self.pd.projectionOf(self.emodes.rawMode(6)).norm()
        self.assertAlmostEqual(p, 0.0283570029711)
        p = self.pd.projectionOf(self.emodes.rawMode(7)).norm()
        self.assertAlmostEqual(p, 0.00570349238657)
        p = self.pd.projectionOf(self.emodes.rawMode(8)).norm()
        self.assertAlmostEqual(p, 0.132182279933)
        p = self.pd.projectionOf(self.emodes.rawMode(9)).norm()
        self.assertAlmostEqual(p, 0.372076699964)
        p = self.pd.projectionOf(self.emodes.rawMode(10)).norm()
        self.assertAlmostEqual(p, 0.294163776112)
        p = self.pd.projectionOf(self.emodes.rawMode(11)).norm()
        self.assertAlmostEqual(p, 0.183147528471)
        p = self.pd.projectionOf(self.emodes.rawMode(12)).norm()
        self.assertAlmostEqual(p, 0.0271240228062)
        p = self.pd.projectionOf(self.emodes.rawMode(13)).norm()
        self.assertAlmostEqual(p, 0.178665480086)
        p = self.pd.projectionOf(self.emodes.rawMode(14)).norm()
        self.assertAlmostEqual(p, 0.128834250691)
        p = self.pd.projectionOf(self.emodes.rawMode(15)).norm()
        self.assertAlmostEqual(p, 0.0591397692555)
        p = self.pd.projectionOf(self.emodes.rawMode(16)).norm()
        self.assertAlmostEqual(p, 0.198830357382)
        p = self.pd.projectionOf(self.emodes.rawMode(17)).norm()
        self.assertAlmostEqual(p, 0.0421877341135)
        p = self.pd.projectionOf(self.emodes.rawMode(18)).norm()
        self.assertAlmostEqual(p, 0.00971216016261)
        p = self.pd.projectionOf(self.emodes.rawMode(19)).norm()
        self.assertAlmostEqual(p, 0.0736856356204)
        p = self.pd.projectionOf(self.emodes.rawMode(20)).norm()
        self.assertAlmostEqual(p, 0.114536706027)
        p = self.pd.projectionOf(self.emodes.rawMode(21)).norm()
        self.assertAlmostEqual(p, 0.274928559325)
        p = self.pd.projectionOf(self.emodes.rawMode(22)).norm()
        self.assertAlmostEqual(p, 0.0686378784385)
        p = self.pd.projectionOf(self.emodes.rawMode(23)).norm()
        self.assertAlmostEqual(p, 0.0130940520067)
        p = self.pd.projectionOf(self.emodes.rawMode(24)).norm()
        self.assertAlmostEqual(p, 0.0281201750134)
        p = self.pd.projectionOf(self.emodes.rawMode(25)).norm()
        self.assertAlmostEqual(p, 0.196135377387)
        p = self.pd.projectionOf(self.emodes.rawMode(26)).norm()
        self.assertAlmostEqual(p, 0.0194536676741)
        p = self.pd.projectionOf(self.emodes.rawMode(27)).norm()
        self.assertAlmostEqual(p, 0.0120311581935)
        p = self.pd.projectionOf(self.emodes.rawMode(28)).norm()
        self.assertAlmostEqual(p, 0.0479444837967)
        p = self.pd.projectionOf(self.emodes.rawMode(29)).norm()
        self.assertAlmostEqual(p, 0.113826710408)
        p = self.pd.projectionOf(self.emodes.rawMode(30)).norm()
        self.assertAlmostEqual(p, 0.0346732822148)
        p = self.pd.projectionOf(self.emodes.rawMode(31)).norm()
        self.assertAlmostEqual(p, 0.0693022722289)
        p = self.pd.projectionOf(self.emodes.rawMode(32)).norm()
        self.assertAlmostEqual(p, 0.0427314320652)
        p = self.pd.projectionOf(self.emodes.rawMode(33)).norm()
        self.assertAlmostEqual(p, 0.0145340306202)
        p = self.pd.projectionOf(self.emodes.rawMode(34)).norm()
        self.assertAlmostEqual(p, 0.0655196353298)
        p = self.pd.projectionOf(self.emodes.rawMode(35)).norm()
        self.assertAlmostEqual(p, 0.0446686455027)
        p = self.pd.projectionOf(self.emodes.rawMode(36)).norm()
        self.assertAlmostEqual(p, 0.0447414683339)
        p = self.pd.projectionOf(self.emodes.rawMode(37)).norm()
        self.assertAlmostEqual(p, 0.105694607933)
        p = self.pd.projectionOf(self.emodes.rawMode(38)).norm()
        self.assertAlmostEqual(p, 0.0369811368042)
        p = self.pd.projectionOf(self.emodes.rawMode(39)).norm()
        self.assertAlmostEqual(p, 0.099422515813)
        p = self.pd.projectionOf(self.emodes.rawMode(40)).norm()
        self.assertAlmostEqual(p, 0.214824447562)
        p = self.pd.projectionOf(self.emodes.rawMode(41)).norm()
        self.assertAlmostEqual(p, 0.0113247783174)
        p = self.pd.projectionOf(self.emodes.rawMode(42)).norm()
        self.assertAlmostEqual(p, 0.0428190611544)
        p = self.pd.projectionOf(self.emodes.rawMode(43)).norm()
        self.assertAlmostEqual(p, 0.0395942495882)
        p = self.pd.projectionOf(self.emodes.rawMode(44)).norm()
        self.assertAlmostEqual(p, 0.184574456157)
        p = self.pd.projectionOf(self.emodes.rawMode(45)).norm()
        self.assertAlmostEqual(p, 0.0389662995269)
        p = self.pd.projectionOf(self.emodes.rawMode(46)).norm()
        self.assertAlmostEqual(p, 0.0443754728854)
        p = self.pd.projectionOf(self.emodes.rawMode(47)).norm()
        self.assertAlmostEqual(p, 0.0189823224475)
        p = self.pd.projectionOf(self.emodes.rawMode(48)).norm()
        self.assertAlmostEqual(p, 0.106332137903)
        p = self.pd.projectionOf(self.emodes.rawMode(49)).norm()
        self.assertAlmostEqual(p, 0.182336029271)
        p = self.pd.projectionOf(self.emodes.rawMode(50)).norm()
        self.assertAlmostEqual(p, 0.222418516151)
        p = self.pd.projectionOf(self.emodes.rawMode(51)).norm()
        self.assertAlmostEqual(p, 0.208379764186)
        p = self.pd.projectionOf(self.emodes.rawMode(52)).norm()
        self.assertAlmostEqual(p, 0.161707260103)
        p = self.pd.projectionOf(self.emodes.rawMode(53)).norm()
        self.assertAlmostEqual(p, 0.29786008418)
        p = self.pd.projectionOf(self.emodes.rawMode(54)).norm()
        self.assertAlmostEqual(p, 0.00783831845487)
        p = self.pd.projectionOf(self.emodes.rawMode(55)).norm()
        self.assertAlmostEqual(p, 0.0265930190584)
        p = self.pd.projectionOf(self.emodes.rawMode(56)).norm()
        self.assertAlmostEqual(p, 0.00543966569196)
        p = self.pd.projectionOf(self.emodes.rawMode(57)).norm()
        self.assertAlmostEqual(p, 0.00405705977917)
        p = self.pd.projectionOf(self.emodes.rawMode(58)).norm()
        self.assertAlmostEqual(p, 0.0622064557348)
        p = self.pd.projectionOf(self.emodes.rawMode(59)).norm()
        self.assertAlmostEqual(p, 0.0489324788794)
        p = self.pd.projectionOf(self.emodes.rawMode(60)).norm()
        self.assertAlmostEqual(p, 0.199110460545)
        p = self.pd.projectionOf(self.emodes.rawMode(61)).norm()
        self.assertAlmostEqual(p, 0.148263694589)
        p = self.pd.projectionOf(self.emodes.rawMode(62)).norm()
        self.assertAlmostEqual(p, 0.14306616768)
        p = self.pd.projectionOf(self.emodes.rawMode(63)).norm()
        self.assertAlmostEqual(p, 0.0798824699709)
        p = self.pd.projectionOf(self.emodes.rawMode(64)).norm()
        self.assertAlmostEqual(p, 0.021465578636)
        p = self.pd.projectionOf(self.emodes.rawMode(65)).norm()
        self.assertAlmostEqual(p, 0.0166782530763)


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(RigidBodyTest))
    s.addTest(loader.loadTestsFromTestCase(PeptideNormalModeTest))
    return s


if __name__ == '__main__':
    unittest.main()
