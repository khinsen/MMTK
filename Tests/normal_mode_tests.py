# Normal mode tests
#
# Written by Konrad Hinsen
# last revision: 2007-3-19
#

import unittest
import MMTK
from MMTK.ForceFields import HarmonicForceField
import MMTK.NormalModes

class WaterTest(unittest.TestCase):

    """
    Test MMTK.NormalModes.EnergeticModes
    """

    def setUp(self):
        self.universe = MMTK.InfiniteUniverse(HarmonicForceField())
        self.universe.water = MMTK.Molecule('water')

    def test_energeticModes(self):
        emodes = MMTK.NormalModes.EnergeticModes(self.universe)
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
        vmodes = MMTK.NormalModes.VibrationalModes(self.universe)
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

if __name__ == '__main__':
    unittest.main()
