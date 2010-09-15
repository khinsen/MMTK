# Energy tests
#
# Written by Konrad Hinsen
# last revision: 2010-9-7
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

class AmberPathIntegralTest(unittest.TestCase):
    """
    Test bonded energies with the Amber99 Forcefield
    with no parameters
    """

    def setUp(self):
        self.universe = InfiniteUniverse(Amber99ForceField())

    def test_pi_quantum_spring_potential(self):
        self.universe.h = Atom('H',amber_atom_type="CT",amber_charge=0.)
        temperature = 100*Units.K
        self.universe.addObject(Environment.PathIntegrals(temperature))
        nb=4
        self.universe.h.setNumberOfBeads(nb) 
        k = (Units.k_B*temperature)*(Units.k_B*temperature)*float(nb*nb)*self.universe.h.mass() / (Units.hbar*Units.hbar*2.)
        beadpos=[Vector(-2.,-3.,-4.),Vector(0.,0.,0.),Vector(5.,1.,3.),Vector(1,-2,3)]
        self.universe.h.setBeadPositions(beadpos)
        numerical_vqu = 0.
        for bead in range(nb): 
            dr = (beadpos[bead]-beadpos[(bead+1)%nb]).length()
            numerical_vqu += k*dr*dr

        piterms = self.universe.energyTerms()
        self.assertAlmostEqual(piterms['path integral spring'], numerical_vqu, 10)

    def test_pi_water_energy(self):
        self.universe.water = Molecule('water')
        e0, grad = self.universe.energyAndGradients()
        for a in self.universe.atomList():
            a.setNumberOfBeads(4)
        self.universe.addObject(Environment.PathIntegrals(100*Units.K))
        piterms = self.universe.energyTerms()
        self.assertEqual(piterms['path integral spring'], 0.)
        self.assertAlmostEqual(piterms['harmonic bond'], e0, 10)
        #e1, grad1 = self.universe.energyAndGradients()
        #self.assertEqual(e0, e1)
        #for a in self.universe.atomList():
        #    self.assert_((grad1[a]-grad[a]).length() == 0.)


class NonbondedListTest:

    def test_nonbondedList(self):

        self.universe.configuration()
        atoms = self.universe.atomList()
        atom_indices = N.array([a.index for a in self.universe.atomList()])
        empty = N.zeros((0, 2), N.Int)

        for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]:

            nblist = NonbondedList(empty, empty, atom_indices,
                                   self.universe._spec, cutoff)
            nblist.update(self.universe.configuration().array)
            distances = nblist.pairDistances()
            pairs1 = nblist.pairIndices()
            pairs1 = [sorted_tuple(pairs1[i]) for i in range(len(pairs1))
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

            
def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(DihedralTest))
    s.addTest(loader.loadTestsFromTestCase(InfiniteUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(OrthorhombicUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(ParallelepipedicUniverseNonbondedListTest))
    s.addTest(loader.loadTestsFromTestCase(LennardJonesSubsetTest))
    s.addTest(loader.loadTestsFromTestCase(AmberPathIntegralTest))
    return s


if __name__ == '__main__':
    unittest.main()
