# Basic tests
#
# Written by Konrad Hinsen
# last revision: 2008-10-29
#

import unittest
import MMTK

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
