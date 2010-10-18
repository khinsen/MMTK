# Ewald summation tests
#
# Written by Konrad Hinsen
#

import unittest
from MMTK import *
from MMTK.ForceFields import Amber99ForceField
from MMTK.Environment import PathIntegrals
from Scientific import N

class EwaldTest(unittest.TestCase):

    def test_selfTerm(self):
        beta = 1.
        real_cutoff = 1.
        reciprocal_cutoff = 1.
        ff = Amber99ForceField(1.,
                               {'method': 'ewald',
                                'real_cutoff': real_cutoff,
                                'beta': beta,
                                'reciprocal_cutoff': reciprocal_cutoff})
        for nb in range(1, 10):
            universe = CubicPeriodicUniverse(2., ff)
            universe.addObject(Molecule('water', position=Vector(-0.5, 0., 0.)))
            universe.addObject(Molecule('water', position=Vector( 0.5, 0., 0.)))
            for a in universe.atomIterator():
                a.setNumberOfBeads(nb)
            if nb > 1:
                universe.addObject(PathIntegrals(10.))
            charge_sum = sum(a.charge()**2 for a in universe.atomIterator())
            self_term = -charge_sum*Units.electrostatic_energy *  (beta/N.sqrt(N.pi))
            self.assertAlmostEqual(self_term,
                                   universe.energyTerms()['electrostatic/ewald self term'],
                                   6)

    def test_allTerms(self):
        beta = 1.
        real_cutoff = 1.
        reciprocal_cutoff = 1.

        ff = Amber99ForceField(1.,
                               {'method': 'ewald',
                                'real_cutoff': real_cutoff,
                                'beta': beta,
                                'reciprocal_cutoff': reciprocal_cutoff})
        es = []
        gs = []
        for nb in range(1, 10):
            universe = CubicPeriodicUniverse(2., ff)
            universe.addObject(Molecule('water', position=Vector(-0.5, 0., 0.)))
            universe.addObject(Molecule('water', position=Vector( 0.5, 0., 0.)))
            for a in universe[0].atomIterator():
                a.setNumberOfBeads(nb)
            if nb > 1:
                universe.addObject(PathIntegrals(10.))
            e, g = universe.energyAndGradients()
            es.append(universe.energyTerms()['electrostatic'])
            gs.append(N.array([a.numberOfBeads()*g[a.beads()[0]].array
                               for a in universe.atomIterator()]))
        for ev in es[1:]:
            self.assertAlmostEqual(ev, es[0], 6)
        for gv in gs[1:]:
            self.assert_(N.fabs(N.ravel(gv-gs[0])).max() < 1.e-10)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(EwaldTest))
    return s

if __name__ == '__main__':
    unittest.main()
