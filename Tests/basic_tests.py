# Basic tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
from MMTK.Proteins import Protein
from Scientific import N

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


class SuperpositionTest(unittest.TestCase):

    """
    Test the findTransformation method
    """

    def test_rotation_translation(self):
        for m in [MMTK.Molecule('water'), Protein('bala1')]:
            universe = MMTK.InfiniteUniverse()
            universe.addObject(m)
            ref_conf = universe.copyConfiguration()
            universe.translateBy(MMTK.Vector(0.1, -1.3, 1.2))
            universe.rotateAroundOrigin(MMTK.Vector(0., 1., 1.), 0.7)
            tr, rms = universe.findTransformation(ref_conf)
            self.assert_(abs(rms) < 1.e-5)
            universe.applyTransformation(tr)
            self.assert_(universe.rmsDifference(ref_conf) < 1.e-5)
            axis, angle = tr.rotation().axisAndAngle()
            self.assert_(abs(angle - 0.7) < 1.e-5)
            self.assert_(abs(abs(N.cos(axis.angle(MMTK.Vector(0., 1., 1.))))-1)
                         < 1.e-5)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(GroupOfAtomTest))
    s.addTest(loader.loadTestsFromTestCase(SuperpositionTest))
    return s

if __name__ == '__main__':
    unittest.main()
