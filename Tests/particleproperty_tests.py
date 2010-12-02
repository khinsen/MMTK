#
# Written by Chris Ing
#

import unittest
from MMTK import *
from Scientific import N

class ParticlePropertiesTest(unittest.TestCase):

    def setUp(self):
        self.universe = InfiniteUniverse()
        self.universe.addObject(Atom('H', position=Vector([1.,1.,1.])))
        self.universe.addObject(Atom('C', position=Vector([-1.,-1.,-1.])))
        self.universe.addObject(Atom('C', position=Vector([2.,2.,2.])))
        self.universe.addObject(Atom('C', position=Vector([-2.,-2.,-2.])))

    def test_configuration(self):
        x = 3.
        y = 1.
        z = 2.
        self.universe_o = OrthorhombicPeriodicUniverse((x, y, z))
        self.universe_o.addObject(Atom('C', position=Vector([1.,1.,1.])))
        self.universe_o.addObject(Atom('C', position=Vector([-1.,-1.,-1.])))
        conf = self.universe_o.configuration()

        # Check copy and conversion to box coordinates work as expected
        confcopy = copy(conf)
        confcopy.convertToBoxCoordinates()
        self.assertAlmostEqual(confcopy.array[0, 0], 1./x, 5)
        self.assertAlmostEqual(confcopy.array[0, 1], 1./y, 5)
        self.assertAlmostEqual(confcopy.array[0, 2], 1./z, 5)
        
        confcopy.convertFromBoxCoordinates()
        self.assertAlmostEqual(confcopy.array[0, 0], conf.array[0, 0], 5)
        self.assertAlmostEqual(confcopy.array[0, 1], conf.array[0, 1], 5)
        self.assertAlmostEqual(confcopy.array[0, 2], conf.array[0, 2], 5)

        # Check that an atom with undefined position returns false for this check
        self.universe_o.addObject(Atom('C'))
        self.assert_(self.universe_o.configuration().hasValidPositions() == False)

    def test_basic(self):
        # Check that an error is thrown if data_array has more 
        # or less vectors than universe has objects
        vdata = N.array([(1.,1.,1.),(-1.,-1.,-1.)])
        self.assertRaises(ValueError, ParticleProperties.ParticleVector, \
                          universe=self.universe, data_array=vdata)
        vdata = N.array([(1.,1.,1.),(-1.,-1.,-1.),(1.,1.,1.),(-1.,-1.,-1.),(1.,1.,1.)])
        self.assertRaises(ValueError, ParticleProperties.ParticleVector, \
                          universe=self.universe, data_array=vdata)

        # This isn't a very useful test
        # self.assert_(self.universe.masses().isParticleProperty == True)

    def test_particle_scalar(self):
        mps = self.universe.masses()
        mps.array[3] = 0.5

        # Check that maximum and minimum values work
        self.assertAlmostEqual(mps.maximum(), self.universe.atomList()[2].mass(), 7)
        self.assertAlmostEqual(mps.minimum(), 0.5, 6)
        
        # Check that summing the scalar returns the correct value
        ms = 0.
        msroot = 0.
        for m in mps.array:
            ms += m
            msroot += N.sqrt(m)
        self.assertAlmostEqual(mps.sumOverParticles(), ms, 7)
        self.assertAlmostEqual(N.sum(mps.applyFunction(N.sqrt)), msroot, 7)

    def test_particle_vector(self):
        self.universe.atomList()[0].setMass(2.)
        self.universe.atomList()[1].setMass(1.)
        self.universe.atomList()[2].setMass(1.)
        self.universe.atomList()[3].setMass(1.)
        vdata = N.array([(1.,1.,1.),(-1.,-1.,-1.),(1.,-1.,1.),(-1.,1.,-1.)])
        vdata2 = N.array([(1.,1.,1.),(-2.,-2.,-2.),(2.,-2.,2.),(-2.,2.,-2.)])
        vpv = ParticleProperties.ParticleVector(self.universe, data_array=vdata)
        vpv2 = ParticleProperties.ParticleVector(self.universe, data_array=vdata2)


        # Check that mass * dot product of two vectors is a known value
        self.assertAlmostEqual(vpv.massWeightedDotProduct(vpv2), 24., 7)

        # Check that the sum of velocities is zero
        self.assert_(N.fabs(N.ravel(vpv.sumOverParticles())).max() < 1.e-10)

        # Check that the length of each vector (all the same) is a known value
        self.assert_(N.fabs(N.ravel(vpv.length()-N.sqrt(3))).max() < 1.e-10)

        # Check that the norm of the array is a known value
        self.assertAlmostEqual(vpv.norm(), N.sqrt(12.), 7)

        # Check that norm scaling works for known norms
        norm_ratio = vpv2.norm()/vpv.norm()
        self.assert_(N.fabs(N.ravel(vpv*norm_ratio - vpv.scaledToNorm(vpv2.norm()))).max() < 1.e-10)

        # Check that mass weighting works for known values
        self.assert_(vpv.massWeightedNorm() - N.sqrt(3) < 1e-10)
        self.assertAlmostEqual(vpv.massWeightedDotProduct(vpv), 5.*(vpv.massWeightedNorm())**2, 5)
        self.assertAlmostEqual(vpv.massWeightedNorm(), vpv.norm()/2, 7)

        # Check the dyadic product for a simple case of 1 atom
        self.universe.removeObject(self.universe.atomList()[3])
        self.universe.removeObject(self.universe.atomList()[2])
        self.universe.removeObject(self.universe.atomList()[1])
        vdata = N.array([(1.,2.,3.)])
        vdata2 = N.array([(-1.,2.,-3.)])
        vpv = ParticleProperties.ParticleVector(self.universe, data_array=vdata)
        vpv2 = ParticleProperties.ParticleVector(self.universe, data_array=vdata2)
        vtensor = vpv.dyadicProduct(vpv2).array
        for x in range(3):
            for y in range(3):
                self.assertAlmostEqual(vtensor[0, x, y], vpv.array[0, x]*vpv2.array[0,y], 7)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(ParticlePropertiesTest))
    return s

if __name__ == '__main__':
    unittest.main()
