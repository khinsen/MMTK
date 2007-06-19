# Normal mode tests
#
# Written by Konrad Hinsen
# last revision: 2007-6-19
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
        self.assertAlmostEqual(freq[0],  2.85497617804)
        self.assertAlmostEqual(freq[1],  3.59284024328)
        self.assertAlmostEqual(freq[2],  4.26855357171)
        self.assertAlmostEqual(freq[3],  5.16159136625)
        self.assertAlmostEqual(freq[4],  7.09274180791)
        self.assertAlmostEqual(freq[5],  7.83097776215)
        self.assertAlmostEqual(freq[6],  9.59758521911)
        self.assertAlmostEqual(freq[7],  11.9531585682)
        self.assertAlmostEqual(freq[8],  14.8601127282)
        self.assertAlmostEqual(freq[9],  15.2065248058)
        self.assertAlmostEqual(freq[10], 24.5040873321)
        self.assertAlmostEqual(freq[11], 28.8319251443)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.057820673306)
        self.assertAlmostEqual(vmodes[7].norm(),  0.0531479847048)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0406909341987)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0363663505062)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0268048215741)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0316000296808)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0306960641089)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0189994693992)
        self.assertAlmostEqual(vmodes[14].norm(), 0.011071709684)
        self.assertAlmostEqual(vmodes[15].norm(), 0.0107411277505)
        self.assertAlmostEqual(vmodes[16].norm(), 0.00577428662572)
        self.assertAlmostEqual(vmodes[17].norm(), 0.00547489239119)

    def test_subspaceModesWithExclusion(self):
        excluded = PairDistanceSubspace(self.universe,
                                        [(self.universe.peptide[0][0].O,
                                          self.universe.peptide[0][-1].CH3)])
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=(excluded,
                                                        self.subspace))
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0],  3.5782639181)
        self.assertAlmostEqual(freq[1],  3.97271215735)
        self.assertAlmostEqual(freq[2],  5.13166264054)
        self.assertAlmostEqual(freq[3],  7.08556627871)
        self.assertAlmostEqual(freq[4],  7.76872007158)
        self.assertAlmostEqual(freq[5],  9.59209084127)
        self.assertAlmostEqual(freq[6],  11.6952492359)
        self.assertAlmostEqual(freq[7],  14.6259188161)
        self.assertAlmostEqual(freq[8],  15.1551463297)
        self.assertAlmostEqual(freq[9],  23.3710255654)
        self.assertAlmostEqual(freq[10], 28.6686533921)
        self.assertAlmostEqual(freq[11], 40.6368580968)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.0543573983684)
        self.assertAlmostEqual(vmodes[7].norm(),  0.0461441234471)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0369274232204)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0268331435223)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0320502610566)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0308393261861)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0198837573395)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0112980200309)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0106269306549)
        self.assertAlmostEqual(vmodes[15].norm(), 0.00621915629047)
        self.assertAlmostEqual(vmodes[16].norm(), 0.00542239330934)
        self.assertAlmostEqual(vmodes[17].norm(), 0.00480126559494)

    def test_subspaceModesWithNumericalDifferentiation(self):
        vmodes = NormalModes.VibrationalModes(self.universe,
                                              subspace=self.subspace,
                                              delta=0.01)
        freq = vmodes.frequencies[6:]
        self.assertAlmostEqual(freq[0],  1.72645040264, 5)
        self.assertAlmostEqual(freq[1],  2.9671174865,  5)
        self.assertAlmostEqual(freq[2],  3.99735617583, 5)
        self.assertAlmostEqual(freq[3],  4.8180612886,  5)
        self.assertAlmostEqual(freq[4],  6.64151129102, 5)
        self.assertAlmostEqual(freq[5],  7.64626271892, 5)
        self.assertAlmostEqual(freq[6],  9.57677736649, 5)
        self.assertAlmostEqual(freq[7],  11.5541230127, 5)
        self.assertAlmostEqual(freq[8],  14.5230764485, 5)
        self.assertAlmostEqual(freq[9],  14.8405362256, 5)
        self.assertAlmostEqual(freq[10], 24.4117593029, 5)
        self.assertAlmostEqual(freq[11], 28.7209926889, 5)
        for i in range(len(vmodes)):
            mi = vmodes.rawMode(i)
            norm_sq = mi.dotProduct(mi)
            self.assertAlmostEqual(norm_sq, 1.)
            for j in range(i+1, len(vmodes)):
                overlap = mi.dotProduct(vmodes.rawMode(j))
                self.assert_(overlap < 1.e-14)
        self.assertAlmostEqual(vmodes[6].norm(),  0.100814855236,  5)
        self.assertAlmostEqual(vmodes[7].norm(),  0.062859055838,  5)
        self.assertAlmostEqual(vmodes[8].norm(),  0.0442747883645, 5)
        self.assertAlmostEqual(vmodes[9].norm(),  0.0385239298851, 5)
        self.assertAlmostEqual(vmodes[10].norm(), 0.0284589465079, 5)
        self.assertAlmostEqual(vmodes[11].norm(), 0.0315884822204, 5)
        self.assertAlmostEqual(vmodes[12].norm(), 0.0314460194228, 5)
        self.assertAlmostEqual(vmodes[13].norm(), 0.0199090672857, 5)
        self.assertAlmostEqual(vmodes[14].norm(), 0.0113607276252, 5)
        self.assertAlmostEqual(vmodes[15].norm(), 0.0109610042085, 5)
        self.assertAlmostEqual(vmodes[16].norm(), 0.0057883026904, 5)
        self.assertAlmostEqual(vmodes[17].norm(), 0.0054875425125, 5)

if __name__ == '__main__':
    unittest.main()
