# Elastic network model tests
#
# Written by Konrad Hinsen
# last revision: 2010-9-7
#

import unittest
from MMTK import *
from MMTK.ForceFields import CalphaForceField, DeformationForceField, \
                             AnisotropicNetworkForceField
from MMTK.Proteins import Protein

class SubsetTest(object):

    def test_singleSubset(self):
        outside = set(self.universe.atomList())-set(self.subset1.atomList())
        e, fc = self.universe.energyAndForceConstants(self.subset1)
        for a1 in outside:
            for a2 in outside:
                self.assert_((fc[a1,a2].array == 0).all())
        for a1 in self.subset1.atomList():
            for a2 in self.subset1.atomList():
                self.assert_((fc[a1,a2].array != 0).any())

    def test_twoSubsets(self):
        inside = set(self.subset1.atomList()) | set(self.subset2.atomList())
        outside = set(self.universe.atomList())-inside
        e, fc = self.universe.energyAndForceConstants(self.subset1, self.subset2)
        for a1 in self.subset1.atomList():
            for a2 in self.subset1.atomList():
                if a1 != a2:
                    self.assert_((fc[a1,a2].array == 0).all())
        for a1 in self.subset2.atomList():
            for a2 in self.subset2.atomList():
                if a1 != a2:
                    self.assert_((fc[a1,a2].array == 0).all())
        for a1 in self.subset1.atomList():
            for a2 in self.subset2.atomList():
                self.assert_((fc[a1,a2].array != 0).any())
        for a1 in outside:
            for a2 in outside:
                self.assert_((fc[a1,a2].array == 0).all())


class CalphaFFSubsetTest(unittest.TestCase,
                         SubsetTest):

    def setUp(self):
        self.universe = InfiniteUniverse(CalphaForceField())
        protein = Protein('insulin_calpha')
        self.universe.addObject(protein)        
        self.subset1 = protein[0]
        self.subset2 = protein[1]


class DeformationFFSubsetTest(unittest.TestCase,
                              SubsetTest):

    def setUp(self):
        self.universe = InfiniteUniverse(DeformationForceField(cutoff=5.))
        protein = Protein('insulin_calpha')
        self.universe.addObject(protein)        
        self.subset1 = protein[0]
        self.subset2 = protein[1]


class ANMFFSubsetTest(unittest.TestCase,
                      SubsetTest):

    def setUp(self):
        self.universe = InfiniteUniverse(AnisotropicNetworkForceField(cutoff=5.))
        protein = Protein('insulin_calpha')
        self.universe.addObject(protein)        
        self.subset1 = protein[0]
        self.subset2 = protein[1]


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(CalphaSubsetTest))
    s.addTest(loader.loadTestsFromTestCase(DeformationSubsetTest))
    s.addTest(loader.loadTestsFromTestCase(ANMSubsetTest))
    return s


if __name__ == '__main__':
    unittest.main()
