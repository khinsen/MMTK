# Energy tests
#
# Written by Konrad Hinsen
#

import unittest
from subsets import SubsetTest
from MMTK import *
from MMTK.MoleculeFactory import MoleculeFactory
from MMTK.ForceFields import Amber99ForceField, LennardJonesForceField
from MMTK_forcefield import NonbondedList
from MMTK.Random import randomPointInBox
from MMTK.Geometry import SCLattice
from MMTK.Utility import pairs
from Scientific.Geometry import ex, ey, ez
from Scientific import N
from cStringIO import StringIO
import itertools

factory = MoleculeFactory()
factory.createGroup('dihedral_test')

factory.addAtom('dihedral_test', 'C1', 'C')
factory.addAtom('dihedral_test', 'C2', 'C')
factory.addAtom('dihedral_test', 'C3', 'C')
factory.addAtom('dihedral_test', 'C4', 'C')

factory.addBond('dihedral_test', 'C1', 'C2')
factory.addBond('dihedral_test', 'C2', 'C3')
factory.addBond('dihedral_test', 'C3', 'C4')

factory.setAttribute('dihedral_test', 'C1.amber_atom_type', 'D1')
factory.setAttribute('dihedral_test', 'C2.amber_atom_type', 'D2')
factory.setAttribute('dihedral_test', 'C3.amber_atom_type', 'D3')
factory.setAttribute('dihedral_test', 'C4.amber_atom_type', 'D4')
factory.setAttribute('dihedral_test', 'C1.amber_charge', 0.)
factory.setAttribute('dihedral_test', 'C2.amber_charge', 0.)
factory.setAttribute('dihedral_test', 'C3.amber_charge', 0.)
factory.setAttribute('dihedral_test', 'C4.amber_charge', 0.)

factory.setPosition('dihedral_test', 'C1', Vector(1., 0., 0.))
factory.setPosition('dihedral_test', 'C2', Vector(0., 0., 0.))
factory.setPosition('dihedral_test', 'C3', Vector(0., 0., 1.))
factory.setPosition('dihedral_test', 'C4', Vector(1., 0., 1.))


def sorted_tuple(pair):
    if pair[0] < pair[1]:
        return (pair[0], pair[1])
    else:
        return (pair[1], pair[0])


class DihedralTest(unittest.TestCase):

    def setUp(self):
        self.universe = InfiniteUniverse()
        self.universe.addObject(factory.retrieveMolecule('dihedral_test'))
        self.mod_template = """Amber parameters

MASS
D1 12.0
D2 12.0
D3 12.0
D4 12.0

BOND
D1-D2    0.0    1.000
D2-D3    0.0    1.000
D3-D4    0.0    1.000

ANGL
D1-D2-D3      0.0      90.0
D2-D3-D4      0.0      90.0

DIHEDRAL
D1-D2-D3-D4   1%15.6f%15.6f%15.6f

NONB
  D1          1.0000    0.0000   0.00000
  D2          1.0000    0.0000   0.00000
  D3          1.0000    0.0000   0.00000
  D4          1.0000    0.0000   0.00000
"""

    def _gradientTest(self, delta = 0.0001):
        e0, grad = self.universe.energyAndGradients()
        atoms = self.universe.atomList()
        for a in atoms:
            num_grad = []
            for v in [ex, ey, ez]:
                x = a.position()
                a.setPosition(x+delta*v)
                eplus = self.universe.energy()
                a.setPosition(x-delta*v)
                eminus = self.universe.energy()
                a.setPosition(x)
                num_grad.append(0.5*(eplus-eminus)/delta)
            self.assert_((Vector(num_grad)-grad[a]).length() < 1.e-2)

    def _forceConstantTest(self, delta = 0.0001):
        e0, grad0, fc = self.universe.energyGradientsAndForceConstants()
        atoms = self.universe.atomList()
        for a1, a2 in itertools.chain(itertools.izip(atoms, atoms),
                                      pairs(atoms)):
            num_fc = []
            for v in [ex, ey, ez]:
                x = a1.position()
                a1.setPosition(x+delta*v)
                e_plus, grad_plus = self.universe.energyAndGradients()
                a1.setPosition(x-delta*v)
                e_minus, grad_minus = self.universe.energyAndGradients()
                a1.setPosition(x)
                num_fc.append(0.5*(grad_plus[a2]-grad_minus[a2])/delta)
            num_fc = N.array([a.array for a in num_fc])
            diff = N.fabs(N.ravel(num_fc-fc[a1, a2].array))
            error = N.maximum.reduce(diff)
            self.assert_(error < 5.e-2)
   
    def _dihedralTerm(self, n, phase, V):

        mod_file = self.mod_template % \
                   (V/(Units.kcal/Units.mol), phase/Units.deg, n)
        ff = Amber99ForceField(mod_files=[StringIO(mod_file)])
        self.universe.setForceField(ff)

        param = self.universe.energyEvaluatorParameters()
        i1, i2, i3, i4, n_test, phase_test, V_test = \
            param['cosine_dihedral_term'][0]
        self.assertEqual(n_test, n)
        # The accuracy is no better than five digits because the
        # parameters pass through a text representation.
        self.assertAlmostEqual(phase_test, phase, 5)
        self.assertAlmostEqual(V_test, V, 5)

        two_pi = 2.*N.pi
        m = self.universe[0]
        for angle in N.arange(0., two_pi, 0.1):
            m.C4.setPosition(Vector(N.cos(angle), N.sin(angle), 1.))
            e = self.universe.energyTerms()['cosine dihedral angle']
            da = self.universe.dihedral(m.C1, m.C2, m.C3, m.C4)
            e_ref = V*(1.+N.cos(n*angle-phase))
            self.assertAlmostEqual(angle % two_pi, da % two_pi, 14)
            self.assertAlmostEqual(e, e_ref, 5)
            self._gradientTest()
            self._forceConstantTest()

    def test_dihedral(self):
        for n in [1, 2, 3, 4]:
            for phase in N.arange(0., 6., 0.5):
                for V in [-10., 2.]:
                    self._dihedralTerm(n, phase, V)

class LennardJonesTest(unittest.TestCase):
    """
    Lennard Jones Tests
    """
    def test_periodic(self):
       universe = OrthorhombicPeriodicUniverse((1.,1.,1.),
                                        LennardJonesForceField(0.4))
       universe.addObject(Atom('Ar', position=Vector(0.0, -0.5, -0.1)))
       universe.addObject(Atom('Ar', position=Vector(0.5, 0.0, 0.4)))
       self.assertAlmostEqual(universe.energy(), 0.0, 6)

       universe.atomList()[0].setPosition(Vector(-0.45, 0., -0.))
       universe.atomList()[1].setPosition(Vector(0.3, 0., 0.))
       eperiodic = universe.energy()
       r = universe.distance(universe.atomList()[0], universe.atomList()[1])
       numerical_v = 4*120*Units.K*Units.k_B*(((0.34)/r)**12-((0.34)/r)**6)
       self.assert_(eperiodic > 1.e-15)
       self.assertAlmostEqual(eperiodic, numerical_v, 6)

       universe.atomList()[0].setPosition(Vector(0., 0., 0.))
       universe.atomList()[1].setPosition(Vector(0.21, 0.21, 0.21))
       self.assert_(universe.energy() - eperiodic < 1.e-15)

    def test_pi_periodic(self):
       universe = OrthorhombicPeriodicUniverse((1.,1.,1.),
                                        LennardJonesForceField(0.40))

       universe.addObject(Environment.PathIntegrals(100.*Units.K))
       for x in range(2):
            universe.addObject(Atom('Ar',nbeads=2))
       universe.atomList()[0].setBeadPositions([Vector(0.0,0.,0.),Vector(0.1,0.,0.)])
       universe.atomList()[1].setBeadPositions([Vector(0.5,0.,0.),Vector(-0.4,0.,0.)])
       self.assertAlmostEqual(universe.pathIntegralEnergies()[0], 0.0, 6)

       universe.atomList()[0].setBeadPositions([Vector(-0.45,0.,0.),Vector(0.4,0.,0.)])
       universe.atomList()[1].setBeadPositions([Vector(0.3,0.,0.),Vector(0.45,0.,0.)])
       eperiodic = universe.pathIntegralEnergies()[0]
       r1 = universe.distance(universe.atomList()[0].beadPositions()[0], universe.atomList()[1].beadPositions()[0])
       r2 = universe.distance(universe.atomList()[0].beadPositions()[1], universe.atomList()[1].beadPositions()[1])
       numerical_v = 2*120*Units.K*Units.k_B*(((0.34)/r1)**12-((0.34)/r1)**6)
       numerical_v = numerical_v + 2*120*Units.K*Units.k_B*(((0.34)/r2)**12-((0.34)/r2)**6)
       self.assert_(eperiodic > 1.e-15)
       self.assertAlmostEqual(eperiodic, numerical_v, 4)

       universe.atomList()[0].setBeadPositions([Vector(0.,0.,0.),Vector(0.4,0.,0.)])
       universe.atomList()[1].setBeadPositions([Vector(0.25,0.,0.),Vector(0.45,0.,0.)])
       self.assert_(universe.pathIntegralEnergies()[0] - eperiodic < 1.e-15)

class AmberPathIntegralTest(unittest.TestCase):
    """
    Test bonded energies with the Amber99 Forcefield
    with no parameters
    """

    def setUp(self):
        self.universe = InfiniteUniverse(Amber99ForceField())
        self.temperature = 100.*Units.K
        self.universe.addObject(Environment.PathIntegrals(self.temperature))

    def test_pi_quantum_spring_potential(self):
        self.universe.h = Atom('H',amber_atom_type="CT",amber_charge=0.)
        nb=4
        self.universe.h.setNumberOfBeads(nb)
        kT = Units.k_B*self.temperature
        k = kT*kT*float(nb*nb)*self.universe.h.mass() / (Units.hbar*Units.hbar*2.)
        beadpos=[Vector(-2.,-3.,-4.),Vector(0.,0.,0.),Vector(5.,1.,3.),Vector(1,-2,3)]
        self.universe.h.setBeadPositions(beadpos)
        numerical_vqu = 0.
        for bead in range(nb): 
            dr = (beadpos[bead]-beadpos[(bead+1)%nb]).length()
            numerical_vqu += k*dr*dr

        piterms = self.universe.energyTerms()
        self.assertAlmostEqual(piterms['path integral spring'], numerical_vqu, 10)
        piterms = self.universe.energyTerms(self.universe.h)
        self.assertAlmostEqual(piterms['path integral spring'], numerical_vqu, 10)
        self.universe.h2 = Atom('H', amber_atom_type='CT', amber_charge = 0.,
                                position = Vector(1., 0., 0.))
        piterms = self.universe.energyTerms(self.universe.h)
        self.assertAlmostEqual(piterms['path integral spring'], numerical_vqu, 10)

    def test_pi_water_energy(self):
        self.universe.water = Molecule('water')
        atoms = self.universe.atomList()
        e0, grad = self.universe.energyAndGradients()
        g0 = [grad[a] for a in atoms]
        for a in self.universe.atomList():
            a.setNumberOfBeads(4)
        piterms = self.universe.energyTerms()
        self.assertEqual(piterms['path integral spring'], 0.)
        self.assertAlmostEqual(piterms['harmonic bond'], e0, 10)
        e1, grad = self.universe.energyAndGradients()
        self.assertAlmostEqual(e0, e1, 7)
        g1 = [4.*grad[a.beads()[0]] for a in atoms]
        for i in range(len(atoms)):
            self.assert_((g1[i]-g0[i]).length() == 0.)
        e2, grad = self.universe.energyAndGradients(self.universe.water)
        self.assertAlmostEqual(e0, e2, 7)

    def test_pi_two_water_energy(self):
        self.universe.addObject(Molecule('water', position=Vector(0., 0., 0.)))
        self.universe.addObject(Molecule('water', position=Vector(1., 0., 0.)))
        atoms = self.universe.atomList()
        piterms = self.universe.energyTerms()
        e0, grad = self.universe.energyAndGradients()
        g0 = [grad[a] for a in atoms]
        for a in self.universe.atomList():
            a.setNumberOfBeads(2)
        piterms = self.universe.energyTerms()
        self.assertEqual(piterms['path integral spring'], 0.)
        e1, grad = self.universe.energyAndGradients()
        self.assertAlmostEqual(e0, e1, 7)
        g1 = [2.*grad[a.beads()[0]] for a in atoms]
        for i in range(len(atoms)):
            self.assert_((g1[i]-g0[i]).length() < 1.e-7)

class NonbondedListTest:

    def test_nonbondedList(self):

        self.universe.configuration()
        atoms = self.universe.atomList()
        atom_indices = N.array([a.index for a in atoms])
        empty = N.zeros((0, 2), N.Int)

        for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]:

            nblist = NonbondedList(empty, empty, atom_indices,
                                   1, N.array(len(atoms)*[[0, 1]], N.Int16),
                                   self.universe._spec, cutoff)
            nblist.update(self.universe.configuration().array)
            distances = nblist.pairDistances()
            pairs1 = nblist.pairIndices()
            self.assert_((pairs1[:, 2] == 1).all())
            pairs1 = [sorted_tuple(pairs1[i, :2]) for i in range(len(pairs1))
                      if distances[i] < cutoff]
            pairs1.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            pairs2 = []
            for i in range(len(atoms)):
                for j in range(i+1, len(atoms)):
                    d = self.universe.distance(atoms[i], atoms[j])
                    if d < cutoff:
                        pairs2.append(sorted_tuple((atoms[i].index,
                                                    atoms[j].index)))
            pairs2.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            self.assertEqual(pairs1, pairs2)

    def test_nonbondedList_pi(self):

        for a in self.universe.atomIterator():
            a.setNumberOfBeads(2)
        self.universe.configuration()
        atoms = self.universe.atomList()
        beads = sum((a.beads() for a in atoms), [])
        bead_indices = N.array([b.index for b in beads])
        empty = N.zeros((0, 2), N.Int)

        for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]:

            nblist = NonbondedList(empty, empty, bead_indices,
                                   4, N.array(len(atoms)*[[0, 2], [1, 2]], N.Int16),
                                   self.universe._spec, cutoff)
            nblist.update(self.universe.configuration().array)
            distances = nblist.pairDistances()
            pairs1 = nblist.pairIndices()
            self.assert_((pairs1[:, 2] == 2).all())
            pairs1 = [sorted_tuple(pairs1[i, :2]) for i in range(len(pairs1))
                      if distances[i] < cutoff]
            pairs1.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            pairs2 = []
            for i in range(len(atoms)):
                for j in range(i+1, len(atoms)):
                    d = self.universe.distance(atoms[i], atoms[j])
                    if d < cutoff:
                        pairs2.append(sorted_tuple((atoms[i].index,
                                                    atoms[j].index)))
                        pairs2.append(sorted_tuple((atoms[i].index+1,
                                                    atoms[j].index+1)))
            pairs2.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            self.assertEqual(pairs1, pairs2)

class InfiniteUniverseNonbondedListTest(unittest.TestCase,
                                        NonbondedListTest):
    def setUp(self):
        self.universe = InfiniteUniverse()
        for i in range(100):
            p = randomPointInBox(1.)
            self.universe.addObject(Atom('C', position = p))

class OrthorhombicUniverseNonbondedListTest(unittest.TestCase,
                                            NonbondedListTest):
    def setUp(self):
        self.universe = OrthorhombicPeriodicUniverse((1.5, 0.8, 1.1))
        for i in range(100):
            p = self.universe.boxToRealCoordinates(randomPointInBox(1.))
            self.universe.addObject(Atom('C', position = p))

class ParallelepipedicUniverseNonbondedListTest(unittest.TestCase,
                                                NonbondedListTest):
    def setUp(self):
        a = Vector(1.5, 0.3, 0.)
        b = Vector(-0.1, 0.8, 0.2)
        c = Vector(0., -0.2, 1.1)
        self.universe = ParallelepipedicPeriodicUniverse((a, b, c))
        for i in range(100):
            p = self.universe.boxToRealCoordinates(randomPointInBox(1.))
            self.universe.addObject(Atom('C', position = p))

class LennardJonesSubsetTest(unittest.TestCase,
                             SubsetTest):

    def setUp(self):
        self.universe = OrthorhombicPeriodicUniverse((2., 2., 2.),
                                                     LennardJonesForceField())
        for point in SCLattice(0.5, 4):
            self.universe.addObject(Atom('Ar', position=point))
        self.subset1 = Collection(self.universe.atomList()[:10])
        self.subset2 = Collection(self.universe.atomList()[20:30])

class PathIntegralConsistencyTest(unittest.TestCase):

    def _consistent(self, i, nbi, j, nbj):
        if nbi > nbj:
            i, nbi, j, nbj = j, nbj, i, nbi
        # integer division - needs to be changed for Python 3
        f = nbj/nbi
        return i == j/f

    def test_consistency(self):
        universe = InfiniteUniverse(Amber99ForceField())
        universe.addObject(Molecule('water'))
        universe.addObject(Environment.PathIntegrals(50.*Units.K))
        atoms = universe.atomList()
        atoms[0].setNumberOfBeads(4)
        atoms[1].setNumberOfBeads(2)
        universe.configuration()
        index_to_atom = dict(sum(([(b.index, b.atom) for b in a.beads()]
                                  for a in atoms), []))
        index_to_bead = dict(sum(([(b.index, b.bead_number) for b in a.beads()]
                                  for a in atoms), []))
        try:
            params = universe.energyEvaluatorParameters()
        except ValueError:
            self.fail()
        try:
            ev = universe.energyEvaluator()
        except ValueError:
            self.fail()
        # first, test the test function
        self.assertTrue(self._consistent(0, 2, 0, 4))
        self.assertTrue(self._consistent(0, 2, 1, 4))
        self.assertFalse(self._consistent(0, 2, 2, 4))
        self.assertFalse(self._consistent(0, 2, 3, 4))
        # next, do the real test
        for i, j, d0, k in params['harmonic_distance_term']:
            ai = index_to_atom[i]
            aj = index_to_atom[j]
            if ai == aj:
                self.assert_((i-j+1) % ai.numberOfBeads() == 0)
            else:
                self.assert_(self._consistent(index_to_bead[i], ai.numberOfBeads(),
                                              index_to_bead[j], aj.numberOfBeads()))

    def test_inconsistency1(self):
        universe = InfiniteUniverse(Amber99ForceField())
        universe.addObject(Molecule('water'))
        universe.addObject(Environment.PathIntegrals(50.*Units.K))
        atoms = universe.atomList()
        atoms[0].setNumberOfBeads(4)
        atoms[1].setNumberOfBeads(3)
        self.assertRaises(ValueError, universe.energyEvaluatorParameters)
        self.assertRaises(ValueError, universe.energyEvaluator)

    def test_inconsistency2(self):
        universe = InfiniteUniverse(Amber99ForceField())
        universe.addObject(Molecule('water'))
        universe.addObject(Environment.PathIntegrals(50.*Units.K))
        atoms = universe.atomList()
        atoms[0].setNumberOfBeads(2)
        atoms[1].setNumberOfBeads(3)
        atoms[2].setNumberOfBeads(6)
        self.assertRaises(ValueError, universe.energyEvaluatorParameters)
        self.assertRaises(ValueError, universe.energyEvaluator)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(DihedralTest))
    s.addTest(loader.loadTestsFromTestCase(InfiniteUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(OrthorhombicUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(ParallelepipedicUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(LennardJonesSubsetTest))
    s.addTest(loader.loadTestsFromTestCase(AmberPathIntegralTest))
    s.addTest(loader.loadTestsFromTestCase(PathIntegralConsistencyTest))
    s.addTest(loader.loadTestsFromTestCase(LennardJonesTest))
    return s


if __name__ == '__main__':
    unittest.main()
