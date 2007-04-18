# Normal mode tests
#
# Written by Konrad Hinsen
# last revision: 2007-3-22
#

import unittest
import MMTK
from MMTK.Proteins import Protein
from MMTK.ForceFields import HarmonicForceField
from MMTK.Subspace import RigidMotionSubspace, PairDistanceSubspace
from MMTK import NormalModes

class WaterTest(unittest.TestCase):

    """
    Test NormalModes.EnergeticModes
    and NormalModes.VibrationalModes
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(HarmonicForceField())
        self.universe.water = MMTK.Molecule('water')

    def test_energeticModes(self):
        emodes = NormalModes.EnergeticModes(self.universe)
        fc = emodes.force_constants[6:]
        self.assertAlmostEqual(fc[0], 757849.3957485439, 5)
        self.assertAlmostEqual(fc[1], 1041520.6252049773, 5)
        self.assertAlmostEqual(fc[2], 1388251.2, 5)
        for i in range(len(emodes)):
            mi = emodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(emodes)):
                overlap = mi.dotProduct(emodes.rawMode(j))
                self.assert_(overlap < 1.e-15)
        self.assertAlmostEqual(emodes[6].norm(), 0.0025656740328985168)
        self.assertAlmostEqual(emodes[7].norm(), 0.0021885627126988836)
        self.assertAlmostEqual(emodes[8].norm(), 0.0018956532696964624)
        f = emodes.fluctuations()
        self.assertAlmostEqual(f[self.universe.water.O],  3.359215370401356e-06)
        self.assertAlmostEqual(f[self.universe.water.H1], 2.061890142153454e-06)
        self.assertAlmostEqual(f[self.universe.water.H2], 2.061890142153454e-06)
        af = emodes.anisotropicFluctuations()
        self.assertAlmostEqual(f[self.universe.water.O],
                               af[self.universe.water.O].trace())
        self.assertAlmostEqual(f[self.universe.water.H1],
                               af[self.universe.water.H1].trace())
        self.assertAlmostEqual(f[self.universe.water.H2],
                               af[self.universe.water.H2].trace())

    def test_vibrationalModes(self):
        vmodes = NormalModes.VibrationalModes(self.universe)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0], 87.220841117866954)
        self.assertAlmostEqual(freq[1], 112.00521970255059)
        self.assertAlmostEqual(freq[2], 181.05242207954848)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-15)
        self.assertAlmostEqual(vmodes[6].norm(), 0.0038454773367577063)
        self.assertAlmostEqual(vmodes[7].norm(), 0.0030510913286990616)
        self.assertAlmostEqual(vmodes[8].norm(), 0.001953454891033823)
        f = vmodes.fluctuations()
        self.assertAlmostEqual(f[self.universe.water.O], 8.0141737611250569e-08)
        self.assertAlmostEqual(f[self.universe.water.H1], 6.9381391949153065e-06)
        self.assertAlmostEqual(f[self.universe.water.H2], 6.9381391949153023e-06)
        af = vmodes.anisotropicFluctuations()
        self.assertAlmostEqual(f[self.universe.water.O],
                               af[self.universe.water.O].trace())
        self.assertAlmostEqual(f[self.universe.water.H1],
                               af[self.universe.water.H1].trace())
        self.assertAlmostEqual(f[self.universe.water.H2],
                               af[self.universe.water.H2].trace())


class PeptideTest(unittest.TestCase):

    """
    Test VibrationalModes with a RigidMotionSubspace
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(HarmonicForceField())
        self.universe.peptide = Protein('bala1')
        self.subspace = RigidMotionSubspace(self.universe,
                                            self.universe.peptide.residues())

    def test_subspaceModes(self):
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=self.subspace)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0], 2.9043011743224980)
        self.assertAlmostEqual(freq[1], 3.6128808005364310)
        self.assertAlmostEqual(freq[2], 4.3000233773335133)
        self.assertAlmostEqual(freq[3], 5.2090714221161276)
        self.assertAlmostEqual(freq[4], 7.0952710075055210)
        self.assertAlmostEqual(freq[5], 7.9176238459485777)
        self.assertAlmostEqual(freq[6], 9.6411463817591212)
        self.assertAlmostEqual(freq[7], 12.0152992300718431)
        self.assertAlmostEqual(freq[8], 14.9048078046098240)
        self.assertAlmostEqual(freq[9], 15.2183689507424500)
        self.assertAlmostEqual(freq[10], 24.5076529744237135)
        self.assertAlmostEqual(freq[11], 28.8319451574310399)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(), 0.0571232913318908)
        self.assertAlmostEqual(vmodes[7].norm(), 0.0529117541861086)
        self.assertAlmostEqual(vmodes[8].norm(), 0.0400676921877648)
        self.assertAlmostEqual(vmodes[9].norm(), 0.0362476194356468)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0268401314942991)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0320133176798322)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0300364605341464)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0187776117285175)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0110276987544110)
        self.assertAlmostEqual(vmodes[15].norm(), 0.0107808471624157)
        self.assertAlmostEqual(vmodes[16].norm(), 0.0057739025201874)
        self.assertAlmostEqual(vmodes[17].norm(), 0.0054749382544073)

    def test_subspaceModesWithExclusion(self):
        excluded = PairDistanceSubspace(self.universe,
                                        [(self.universe.peptide[0][0].O,
                                          self.universe.peptide[0][-1].CH3)])
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=(excluded,
                                                        self.subspace))
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0], 3.6017105438829620)
        self.assertAlmostEqual(freq[1], 4.0003081312982314)
        self.assertAlmostEqual(freq[2], 5.1866580036742684)
        self.assertAlmostEqual(freq[3], 7.0881826743326561)
        self.assertAlmostEqual(freq[4], 7.8556351790494565)
        self.assertAlmostEqual(freq[5], 9.6339267095045233)
        self.assertAlmostEqual(freq[6], 11.7536434294379255)
        self.assertAlmostEqual(freq[7], 14.6787264425208157)
        self.assertAlmostEqual(freq[8], 15.1577783773569710)
        self.assertAlmostEqual(freq[9], 23.3702950486827703)
        self.assertAlmostEqual(freq[10], 28.6687153116058759)
        self.assertAlmostEqual(freq[11], 40.6387222002291679)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-15)
        self.assertAlmostEqual(vmodes[6].norm(), 0.0539485999212520)
        self.assertAlmostEqual(vmodes[7].norm(), 0.0457018173643255)
        self.assertAlmostEqual(vmodes[8].norm(), 0.0366795025558611)
        self.assertAlmostEqual(vmodes[9].norm(), 0.0268655512968592)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0324186752148948)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0302214878165687)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0196586315981532)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0112857136900571)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0106366256304052)
        self.assertAlmostEqual(vmodes[15].norm(), 0.0062204565671206)
        self.assertAlmostEqual(vmodes[16].norm(), 0.0054224348308254)
        self.assertAlmostEqual(vmodes[17].norm(), 0.0048003670430221)

    def test_subspaceModesWithNumericalDifferentiation(self):
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=self.subspace,
                                              delta=0.01)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0], 1.83856873, 5)
        self.assertAlmostEqual(freq[1], 2.97757135, 5)
        self.assertAlmostEqual(freq[2], 4.02648368, 5)
        self.assertAlmostEqual(freq[3], 4.87021586, 5)
        self.assertAlmostEqual(freq[4], 6.64485664, 5)
        self.assertAlmostEqual(freq[5], 7.74286200, 5)
        self.assertAlmostEqual(freq[6], 9.61511625, 5)
        self.assertAlmostEqual(freq[7], 11.6184651, 5)
        self.assertAlmostEqual(freq[8], 14.5713711, 5)
        self.assertAlmostEqual(freq[9], 14.8473398, 5)
        self.assertAlmostEqual(freq[10], 24.415408, 5)
        self.assertAlmostEqual(freq[11], 28.721024, 5)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.0946584, 5)
        self.assertAlmostEqual(vmodes[7].norm(),  0.0626193, 5)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0436544, 5)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0383367, 5)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0284992, 5)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0319058, 5)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0308663, 5)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0196713, 5)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0113195, 5)
        self.assertAlmostEqual(vmodes[15].norm(), 0.0109957, 5)
        self.assertAlmostEqual(vmodes[16].norm(), 0.0057879, 5)
        self.assertAlmostEqual(vmodes[17].norm(), 0.0054876, 5)

if __name__ == '__main__':
    unittest.main()
