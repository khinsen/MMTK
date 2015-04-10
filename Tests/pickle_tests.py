# Pickle tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from MMTK import NormalModes
from Scientific import N
import os

class PeptideTest(unittest.TestCase):

    """
    Test pickling a small protein
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(Amber99ForceField())
        self.universe.peptide = Protein('bala1')

    def tearDown(self):
        try:
            os.remove('test.pickle')
        except OSError:
            pass

    def test_energy(self):
        # Check that the energy terms after pickling and reloading
        # are the same.
        energy_terms = self.universe.energyTerms()
        MMTK.save(self.universe, 'test.pickle')
        restored_universe = MMTK.load('test.pickle')
        r_energy_terms = restored_universe.energyTerms()
        self.assertEqual(energy_terms.keys(), r_energy_terms.keys())

    def test_normal_modes(self):
        emodes = NormalModes.EnergeticModes(self.universe)
        MMTK.save(emodes, 'test.pickle')
        restored_modes = MMTK.load('test.pickle')
        for i in range(len(emodes)):
            self.assertEqual(emodes[i].force_constant,
                             restored_modes[i].force_constant)
            err = N.minimum.reduce(N.fabs(N.ravel(emodes[i].array
                                                  -restored_modes[i].array)))
            self.assert_(err < 1.e-15)


class WaterTest(unittest.TestCase):

    """
    Test pickling a water molecule
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(Amber99ForceField())
        self.universe.water1 = MMTK.Molecule('water',
                                             position=MMTK.Vector(0.,0.,0.))
        self.universe.water2 = MMTK.Molecule('water',
                                             position=MMTK.Vector(0.5,0.,0.))

    def tearDown(self):
        try:
            os.remove('test.pickle')
        except OSError:
            pass

    def test_energy(self):
        # Check that the energy terms after pickling and reloading
        # are the same.
        energy_terms = self.universe.energyTerms()
        MMTK.save(self.universe, 'test.pickle')
        restored_universe = MMTK.load('test.pickle')
        self.assertEqual(energy_terms, restored_universe.energyTerms())

    def test_parent(self):
        # Check that pickling an object from a universe sets that
        # objects's parent to None, rather than pickling the whole universe.
        MMTK.save(self.universe.water1, 'test.pickle')
        restored_molecule = MMTK.load('test.pickle')
        self.assertEqual(restored_molecule.parent, None)

    def test_type(self):
        # Check that the MoleculeType object is not pickled, but
        # always taken from the database.
        MMTK.save(self.universe.water1, 'test.pickle')
        restored_molecule = MMTK.load('test.pickle')
        self.assertEqual(restored_molecule.type, self.universe.water1.type)


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(PeptideTest))
    s.addTest(loader.loadTestsFromTestCase(WaterTest))
    return s


if __name__ == '__main__':
    unittest.main()
