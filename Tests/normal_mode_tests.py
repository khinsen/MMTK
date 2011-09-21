# Normal mode tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from MMTK.Subspace import RigidMotionSubspace, PairDistanceSubspace
from MMTK import NormalModes

class WaterTest(unittest.TestCase):

    """
    Test NormalModes.EnergeticModes
    and NormalModes.VibrationalModes
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(Amber99ForceField())
        self.universe.water = MMTK.Molecule('water')

    def test_energeticModes(self):
        emodes = NormalModes.EnergeticModes(self.universe)
        fc = emodes.force_constants[6:]
        self.assertAlmostEqual(fc[0], 757849.3957485439, 5)
        self.assertAlmostEqual(fc[1], 1041551.2240706938, 5)
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
        self.assertAlmostEqual(freq[1], 112.01238171677888)
        self.assertAlmostEqual(freq[2], 181.05242207954848)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-15)
        self.assertAlmostEqual(vmodes[6].norm(), 0.0038454773367577063)
        self.assertAlmostEqual(vmodes[7].norm(), 0.0030509249071616175)
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
        self.universe = MMTK.InfiniteUniverse(Amber99ForceField())
        self.universe.peptide = Protein('bala1')
        self.subspace = RigidMotionSubspace(self.universe,
                                            self.universe.peptide.residues())

    def test_subspaceModes(self):
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=self.subspace)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0],  3.29430105)
        self.assertAlmostEqual(freq[1],  3.5216842)
        self.assertAlmostEqual(freq[2],  6.36658335)
        self.assertAlmostEqual(freq[3],  6.81335036)
        self.assertAlmostEqual(freq[4],  8.96368475)
        self.assertAlmostEqual(freq[5],  9.80151378)
        self.assertAlmostEqual(freq[6],  12.26528168)
        self.assertAlmostEqual(freq[7],  12.64098176)
        self.assertAlmostEqual(freq[8],  16.13558466)
        self.assertAlmostEqual(freq[9],  17.15343073)
        self.assertAlmostEqual(freq[10], 24.86662786)
        self.assertAlmostEqual(freq[11], 29.0265789)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-13)
        self.assertAlmostEqual(vmodes[6].norm(),  0.0577716669202)
        self.assertAlmostEqual(vmodes[7].norm(),  0.0571769348227)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0276758947516)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0254182530029)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0213488665278)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0195579088896)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0194026958464)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0251256023118)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0100970642246)
        self.assertAlmostEqual(vmodes[15].norm(), 0.00935723869926)
        self.assertAlmostEqual(vmodes[16].norm(), 0.00561663849164)
        self.assertAlmostEqual(vmodes[17].norm(), 0.00541615037839)

    def test_subspaceModesWithExclusion(self):
        excluded = PairDistanceSubspace(self.universe,
                                        [(self.universe.peptide[0][0].O,
                                          self.universe.peptide[0][-1].CH3)])
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=(excluded,
                                                        self.subspace))
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0], 3.48436894)
        self.assertAlmostEqual(freq[1], 5.36671021)
        self.assertAlmostEqual(freq[2], 6.57184292)
        self.assertAlmostEqual(freq[3], 8.37798098)
        self.assertAlmostEqual(freq[4], 8.95620773)
        self.assertAlmostEqual(freq[5], 12.06606868)
        self.assertAlmostEqual(freq[6], 12.62157723)
        self.assertAlmostEqual(freq[7], 16.02985338)
        self.assertAlmostEqual(freq[8], 16.90566338)
        self.assertAlmostEqual(freq[9], 23.55890937)
        self.assertAlmostEqual(freq[10],28.92170134)
        self.assertAlmostEqual(freq[11],41.55110572)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 2.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.0578155246705)
        self.assertAlmostEqual(vmodes[7].norm(),  0.030148251724)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0275287609872)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0225332103616)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0216566218055)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0215173538578)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0249366816925)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0102016850103)
        self.assertAlmostEqual(vmodes[14].norm(), 0.00939822099626)
        self.assertAlmostEqual(vmodes[15].norm(), 0.00611049609175)
        self.assertAlmostEqual(vmodes[16].norm(), 0.0053756335944)
        self.assertAlmostEqual(vmodes[17].norm(), 0.00475443155674)

    def test_subspaceModesWithNumericalDifferentiation(self):
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=self.subspace,
                                              delta=0.01)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0],  3.31849983,  5)
        self.assertAlmostEqual(freq[1],  3.54902184,  5)
        self.assertAlmostEqual(freq[2],  6.39993263,  5)
        self.assertAlmostEqual(freq[3],  6.82609026,  5)
        self.assertAlmostEqual(freq[4],  8.97880448,  5)
        self.assertAlmostEqual(freq[5],  9.85460764,  5)
        self.assertAlmostEqual(freq[6],  12.34554287, 5)
        self.assertAlmostEqual(freq[7],  12.89029636, 5)
        self.assertAlmostEqual(freq[8],  16.13726961, 5)
        self.assertAlmostEqual(freq[9],  17.15617066, 5)
        self.assertAlmostEqual(freq[10], 24.86718661, 5)
        self.assertAlmostEqual(freq[11], 29.02770494, 5)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.057265850728,   5)
        self.assertAlmostEqual(vmodes[7].norm(),  0.0568157841206,  5)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0272242339256,  5)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0253529888364,  5)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0212445420919,  5)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0188740025114,  5)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0190325360167,  5)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0252070482368,  5)
        self.assertAlmostEqual(vmodes[14].norm(), 0.010099061532,   5)
        self.assertAlmostEqual(vmodes[15].norm(), 0.00935799475998, 5)
        self.assertAlmostEqual(vmodes[16].norm(), 0.00561654885265, 5)
        self.assertAlmostEqual(vmodes[17].norm(), 0.00541642076462, 5)


def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(WaterTest))
    s.addTest(loader.loadTestsFromTestCase(PeptideTest))
    return s


if __name__ == '__main__':
    unittest.main()
