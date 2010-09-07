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
from subsets import SubsetTest

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
