# Basic tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
import MMTK.Utility
from Scientific import N

class TranslationTest(unittest.TestCase):

    def test_atom_translation_and_rotation(self):
        atom = MMTK.Atom('C')
        p = MMTK.Vector(0., 0., 1.)
        atom.setPosition(p)
        universe = MMTK.InfiniteUniverse()
        universe.addObject(atom)
        opos = atom.position()
        displace = MMTK.Vector(1.,2.,3.)
        atom.translateBy(displace)
        npos = atom.position()
        self.assertAlmostEqual(opos[0],npos[0]-1.,5)
        self.assertAlmostEqual(opos[1],npos[1]-2.,5)
        self.assertAlmostEqual(opos[2],npos[2]-3.,5)
        for i in range(3):
            atom.rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(0.,1.,0.),N.pi/3.)
        for i in range(3):
            atom.rotateAroundOrigin(MMTK.Vector(0.,1.,0.),N.pi/3.)
        dpos = atom.position()
        self.assertAlmostEqual(dpos[0],npos[0],5)
        self.assertAlmostEqual(dpos[1],npos[1],5)
        self.assertAlmostEqual(dpos[2],npos[2],5)


    def test_pi_atom_translation_and_rotation(self):
        atom = MMTK.Atom('C')
        np = 3
        atom.setNumberOfBeads(np)
        beadpos=[MMTK.Vector(-2.,-3.,-4.),MMTK.Vector(5.,1.,3.),MMTK.Vector(1,-2,3)]
        atom.setBeadPositions(beadpos)
        universe = MMTK.InfiniteUniverse()
        universe.addObject(atom)
        opos = atom.beadPositions()
        displace = MMTK.Vector(1.,2.,3.)
        atom.translateBy(displace)
        npos = atom.beadPositions()
        for x in range(len(opos)):
            self.assertAlmostEqual(opos[x][0]+1.,npos[x][0],5) 
            self.assertAlmostEqual(opos[x][1]+2.,npos[x][1],5) 
            self.assertAlmostEqual(opos[x][2]+3.,npos[x][2],5) 

        conf = universe.configuration()
        for i in range(6):
            atom.rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(0.,0.,1.),N.pi/3.)
        for i in range(6):
            atom.rotateAroundOrigin(MMTK.Vector(0.,1.,0.),N.pi/3.)
        for i in range(7):
            atom.rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(1.,0.,0.),2*N.pi/7.)
        dpos = atom.beadPositions()
        for x in range(len(opos)):
            self.assertAlmostEqual(dpos[x][0],npos[x][0],5) 
            self.assertAlmostEqual(dpos[x][1],npos[x][1],5) 
            self.assertAlmostEqual(dpos[x][2],npos[x][2],5) 

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
