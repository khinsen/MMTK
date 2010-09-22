# Basic tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
import MMTK.Utility

class PositionTest(unittest.TestCase):

    def test_atom_position(self):
        atom = MMTK.Atom('C')
        atom.index = 42
        universe = MMTK.InfiniteUniverse()
        universe.addObject(atom)
        conf = universe.configuration()
        self.assert_(atom.index == 0)
        self.assert_((conf.array > MMTK.Utility.undefined_limit).all())
        self.assert_(atom.position() is None)
        p = MMTK.Vector(0., 0., 1.)
        atom.setPosition(p)
        self.assert_(atom.position() == p)
        
    def test_pi_atom_position(self):
        atom = MMTK.Atom('C')
        atom.index = 42
        atom.setNumberOfBeads(3)
        universe = MMTK.InfiniteUniverse()
        universe.addObject(atom)
        self.assert_(universe.numberOfPoints() == 3)
        conf = universe.configuration()
        self.assert_(atom.index == 0)
        self.assert_((conf.array > MMTK.Utility.undefined_limit).all())
        self.assert_(atom.position() is None)
        p = MMTK.Vector(0., 0., 1.)
        atom.setPosition(p)
        self.assert_(atom.position() == p)
        self.assert_(atom.beadPositions() == 3*[p])
        p = [MMTK.Vector(0., 1., 0.), MMTK.Vector(1., 0., 0.), MMTK.Vector(0., 0., 1.)]
        atom.setBeadPositions(p)
        self.assert_(atom.beadPositions() == p)

class GroupOfAtomTest(unittest.TestCase):

    """
    Test the methods of Collections.GroupOfAtoms
    """

    def setUp(self):
        self.molecule = MMTK.Molecule('water')
        self.results = {}
        self.results['numberOfAtoms'] = 3

    def test_numberOfAtoms(self):
        self.assertEqual(self.molecule.numberOfAtoms(),
                         self.results['numberOfAtoms'])


def suite():
    return unittest.TestLoader().loadTestsFromTestCase(GroupOfAtomTest)

if __name__ == '__main__':
    unittest.main()
