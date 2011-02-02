# Basic tests
#
# Written by Konrad Hinsen
#

import unittest
import MMTK
import MMTK.Utility
from MMTK.Proteins import Protein
from Scientific import N

class GroupOfAtomsTest(unittest.TestCase):

    def test_moment_of_inertia(self):
        col1 = MMTK.Collection()
        pos1 = MMTK.Vector(1.,0.,0.)
        pos2 = MMTK.Vector(-1.,0.,0.)
        col1.addObject(MMTK.Atom('C',position=pos1))
        col1.addObject(MMTK.Atom('C',position=pos2))
        mass1 = col1[0].mass()
        mass2 = col1[1].mass()

        #Test if the trace of the moment of inertia tensor is a known value
        #Linear molecule on x-axis so only I_yy and I_zz are non-zero equal to:
        i_num = (mass1*(pos1[0])**2) + (mass2*(pos2[0])**2)
        t = col1.centerAndMomentOfInertia()[1]
        self.assertAlmostEqual(N.sum(t.array.flat), i_num*2., 7)

    def test_center_of_mass(self):
        universe = MMTK.InfiniteUniverse()
        for x in range(6):
            universe.addObject(MMTK.Atom('C'))

        universe.atomList()[0].setMass(10.0*MMTK.Units.amu)
        universe.atomList()[1].setMass(1.0*MMTK.Units.amu)
        universe.atomList()[0].setPosition(MMTK.Vector(1.,0.,0.))
        universe.atomList()[1].setPosition(MMTK.Vector(-10.,0.,0.))
        universe.atomList()[2].setPosition(MMTK.Vector(0.,1.,0.))
        universe.atomList()[3].setPosition(MMTK.Vector(0.,-1.,0.))
        universe.atomList()[4].setPosition(MMTK.Vector(0.,0.,1.))
        universe.atomList()[5].setPosition(MMTK.Vector(0.,0.,-1.))

        #Test the center of mass for the collection is zero, also for the universe
        self.assert_((MMTK.Collection(universe.atomList()).centerOfMass().length() < 1.e-7))
        self.assert_(universe.centerOfMass() == MMTK.Vector(0.0,0.0,0.0))

    def test_pi_center_of_mass(self):
        universe = MMTK.InfiniteUniverse()
        for x in range(4):
            universe.addObject(MMTK.Atom('C',nbeads=2))

        #Test the center of mass for the collection is zero, also for the universe
        universe.atomList()[0].setBeadPositions([MMTK.Vector(-1.0,0.,0.),MMTK.Vector(-1.0,0.,0.)])
        universe.atomList()[1].setBeadPositions([MMTK.Vector(0.5,0.,0.),MMTK.Vector(1.5,0.,0.)])
        universe.atomList()[2].setBeadPositions([MMTK.Vector(0.,0.5,0.),MMTK.Vector(0.,-1.5,0.)])
        universe.atomList()[3].setBeadPositions([MMTK.Vector(0.,-0.5,0.),MMTK.Vector(0.,1.5,0.)])
        self.assert_((MMTK.Collection(universe.atomList()).centerOfMass().length() < 1.e-10))
        self.assert_((universe.centerOfMass().length() < 1.e-10))
        
        #Test that a non-symmetrical collection of atoms has a center of mass that is zero
        universe.atomList()[0].setMass(10.0*MMTK.Units.amu)
        universe.atomList()[1].setMass(1.0*MMTK.Units.amu)
        universe.atomList()[0].setBeadPositions([MMTK.Vector(-40.,0.,0.),MMTK.Vector(41.,0.,0.)])
        universe.atomList()[1].setBeadPositions([MMTK.Vector(-21.,0.,0.),MMTK.Vector(11.,0.,0.)])
        universe.atomList()[2].setBeadPositions([MMTK.Vector(0.,0.5,0.),MMTK.Vector(0.,-1.5,0.)])
        universe.atomList()[3].setBeadPositions([MMTK.Vector(0.,-0.5,0.),MMTK.Vector(0.,1.5,0.)])
        self.assert_((MMTK.Collection(universe.atomList()).centerOfMass().length() < 1.e-10))

    def test_find_transformation(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C'))
        universe.addObject(MMTK.Atom('C'))
        universe.atomList()[0].setPosition(MMTK.Vector(1.,0.,0.))
	universe.atomList()[1].setPosition(MMTK.Vector(0.,0.,0.))
        conf1 = universe.copyConfiguration()
        
        dispvec = MMTK.Vector(5.,5.,5.)
        universe[0].translateBy(dispvec)
        universe[1].translateBy(dispvec)
        universe[0].rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(0.,1.,0.),N.pi)
        universe[1].rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(0.,1.,0.),N.pi)
        trans = universe.findTransformation(conf1)

        #Test that we can recover the rotation by asserting the trace -1 +1 -1
        self.assert_(trans[0].tensor.array.trace() + 1.0 < 1e-5)
        #Test that we can recover the translation by asserting the sum of the displacement vector is -15
        self.assert_(15.0 + trans[0].vector.array.sum() < 1.e-5)

# not supported
#    def test_pi_find_transformation(self):
#        universe = MMTK.InfiniteUniverse()
#        universe.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector(0.,0.,0.)))
#        universe.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector(1.,1.,1.)))
#        conf1 = universe.copyConfiguration()
#        
#        universe[0].setBeadPositions([MMTK.Vector(10.,10.,10.),MMTK.Vector(-10.,-10.,-10.)])
#        universe[1].setBeadPositions([MMTK.Vector(11.,11.,11.),MMTK.Vector(-9.,-9.,-9.)])
#        trans = universe.findTransformation(conf1)
#        print trans[0].vector.array
#        print trans[0].tensor.array

    def test_rms_diff(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C',position=MMTK.Vector(0.,0.,0.)))
        universe.addObject(MMTK.Atom('H',position=MMTK.Vector(1.,1.,1.)))
        conf1 = universe.copyConfiguration()

        #Test the rms is zero for identical conformations
        self.assert_(universe.rmsDifference(conf1) < 1e-5)

        #Test the rms numerical value for a known translation
        universe[0].translateBy(MMTK.Vector(-1.,3.,-2.))
        universe[1].translateBy(MMTK.Vector(-3.,1.,2.))
        self.assert_(universe.rmsDifference(conf1) - N.sqrt(14.) < 1e-5)

    def test_pi_rms_diff(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C', nbeads=2))
        universe.addObject(MMTK.Atom('H', nbeads=2))
        universe[0].setBeadPositions([MMTK.Vector(10.,10.,10.),MMTK.Vector(-10.,-10.,-10.)])
        universe[1].setBeadPositions([MMTK.Vector(11.,11.,11.),MMTK.Vector(-9.,-9.,-9.)])
        conf1 = universe.copyConfiguration()
        #Test the rms is zero for identical conformations
        self.assert_(universe.rmsDifference(conf1) < 1e-5)

        #Test that the average of the bead positions still results in zero rms
        universe[0].setBeadPositions([MMTK.Vector(20.,40.,5.5),MMTK.Vector(-20.,-40.,-5.5)])
        universe[1].setBeadPositions([MMTK.Vector(16.,16.,16.),MMTK.Vector(-14.,-14.,-14.)])
        self.assert_(universe.rmsDifference(conf1) < 1e-5)

        #Test the rms numerical value for a known translation
        universe[0].translateBy(MMTK.Vector(-1.,3.,-2.))
        universe[1].translateBy(MMTK.Vector(-3.,1.,2.))
        self.assert_(universe.rmsDifference(conf1) - N.sqrt(14) < 1e-5)

    def test_bounding(self):
        universe = MMTK.InfiniteUniverse()
        for x in range(3):
            universe.addObject(MMTK.Atom('C'))
        universe.atomList()[0].setPosition(MMTK.Vector(1.,0.,0.))
	universe.atomList()[1].setPosition(MMTK.Vector(0.,0.,1.))
        universe.atomList()[2].setPosition(MMTK.Vector(0.,1.,0.))
        bb = universe.boundingBox()

        #Test the dimensions of the box that contains the three atoms
        self.assert_((bb[0].length() < 1.e-5)) 
        self.assert_((bb[1].length() - N.sqrt(3) < 1.e-5)) 
 
        universe.atomList()[0].setPosition(MMTK.Vector(0.,0.,0.))
        universe.atomList()[1].setPosition(MMTK.Vector(-1.,0.,0.))
        universe.atomList()[2].setPosition(MMTK.Vector(1.,0.,0.))
        bs = universe.boundingSphere()

        #Test the center of the sphere is the origin and that the radius is 1
        self.assertAlmostEquals(bs.center.length(), 0., 10) 
        self.assertAlmostEquals(bs.radius, 1., 10) 

    def test_pi_bounding(self):
        universe = MMTK.InfiniteUniverse()
        for x in range(3):
            universe.addObject(MMTK.Atom('C',nbeads=2))
        universe.atomList()[0].setBeadPositions([MMTK.Vector(0.,0.,0.5),MMTK.Vector(1.,0.,0.)])
        universe.atomList()[1].setBeadPositions([MMTK.Vector(0.5,0.,0.),MMTK.Vector(0.,1.,0.)])
        universe.atomList()[2].setBeadPositions([MMTK.Vector(0.,0.5,0.),MMTK.Vector(0.,0.,1.)])
        bb = universe.boundingBox()

        #Test the dimensions of the box that contains the three atoms
        self.assert_((bb[0].length() < 1.e-5)) 
        self.assert_((bb[1].length() - N.sqrt(3) < 1.e-5)) 
 
        universe.atomList()[0].setBeadPositions([MMTK.Vector(0.,0.,-0.1),MMTK.Vector(0.,0.,0.1)])
        universe.atomList()[1].setBeadPositions([MMTK.Vector(-1.,0.,0.),MMTK.Vector(0.,1.,0.)])
        universe.atomList()[2].setBeadPositions([MMTK.Vector(1.,0.,0.),MMTK.Vector(0.,-1.,0.)])
        bs = universe.boundingSphere()

        #Test the center of the sphere is the origin and that the radius is 1
        self.assertAlmostEquals(bs.center.length(), 0., 10) 
        self.assertAlmostEquals(bs.radius, 1., 10) 

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

        #Test translations move an atom by a known magnitude
        self.assertAlmostEqual(opos[0],npos[0]-1.,5)
        self.assertAlmostEqual(opos[1],npos[1]-2.,5)
        self.assertAlmostEqual(opos[2],npos[2]-3.,5)
        for i in range(3):
            atom.rotateAroundAxis(MMTK.Vector(0.,0.,0.),MMTK.Vector(0.,1.,0.),N.pi/3.)
        for i in range(3):
            atom.rotateAroundOrigin(MMTK.Vector(0.,1.,0.),N.pi/3.)
        dpos = atom.position()

        #Test rotations rotate an atom back to its original position 
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

        #Test translations move an atom by a known magnitude
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

        #Test rotations rotate an atom back to its original position 
        for x in range(len(opos)):
            self.assertAlmostEqual(dpos[x][0],npos[x][0],5) 
            self.assertAlmostEqual(dpos[x][1],npos[x][1],5) 
            self.assertAlmostEqual(dpos[x][2],npos[x][2],5) 

    def test_boolean_mask(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C'))
        col1 = MMTK.Collection() 
        for x in range(3):
            col1.addObject(MMTK.Atom('C'))
        universe.addObject(col1)
        universe.addObject(MMTK.Atom('C'))
        bmask = col1.booleanMask()

        #Test that the last carbon is not part of the boolean mask, and the other elements sum to zero
        self.assert_(bmask[4] == 0)
        self.assert_(N.sum(bmask) == 3)

#not supported
#
#    def test_pi_boolean_mask(self):
#        universe = MMTK.InfiniteUniverse()
#        universe.addObject(MMTK.Atom('C'))
#        col1 = MMTK.Collection() 
#        for x in range(3):
#            col1.addObject(MMTK.Atom('C', nbeads=2))
#        universe.addObject(col1)
#        universe.addObject(MMTK.Atom('C'))
#        print col1.booleanMask()
#        self.assert_(col1.booleanMask()[4] == 0)

    def test_kinetic(self):
        universe = MMTK.InfiniteUniverse()
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', position=MMTK.Vector([0.,0.,0.])))
        col1.addObject(MMTK.Atom('H', position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        vbeads = N.array([(-1.,1.,-1.),(1.,-1.,1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        ke1 = universe.kineticEnergy()
        ke2 = col1.kineticEnergy()
        ke3 = 3*(universe[0]._mass+universe[1]._mass)/2

        #Test numerically calculated kinetic energy against universe and collection methods
        self.assertAlmostEqual(ke3, ke2, 7)
        self.assertAlmostEqual(ke1, ke2, 7)

    def test_pi_kinetic(self):
        universe = MMTK.InfiniteUniverse()
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector([0.,0.,0.])))
        col1.addObject(MMTK.Atom('H', nbeads=2, position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        vbeads = N.array([(1.,1.,-1.),(-1.,1.,1.),(1.,-1.,1.),(-1.,1.,1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        ke1 = universe.kineticEnergy()
        ke2 = col1.kineticEnergy()
        ke3 = 3*(universe[0]._mass+universe[1]._mass)/2

        #Test numerically calculated kinetic energy against universe and collection methods
        self.assertAlmostEqual(ke3, ke2, 7)
        self.assertAlmostEqual(ke1, ke2, 7)

    def test_temperature(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C', position=MMTK.Vector([-1.,-1.,-1.])))
        universe.addObject(MMTK.Atom('C', position=MMTK.Vector([-2.,-2.,-2.])))
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', position=MMTK.Vector([0.,0.,0.])))
        col1.addObject(MMTK.Atom('C', position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        #Test numerical temperature of a known system
        vdata = N.array([(1.,1.,1.),(-1.,-1.,-1.), (1.,-1.,1.),(1.,1.,-1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vdata)
        universe.setVelocities(v)
        temp1 = col1.temperature()
        self.assertAlmostEqual((universe[0].mass()+universe[1].mass())/(2.*MMTK.Units.k_B), temp1, 6)
        self.assertAlmostEqual(universe.temperature(), temp1, 6)

        #Test if different velocity atoms result in an averaged temperature
        vdata = N.array([(1.,1.,1.),(-1.,-1.,-1.), (2.,-2.,2.),(2.,2.,-2.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vdata)
        universe.setVelocities(v)
        temp2 = col1.temperature()
        self.assertAlmostEqual(universe.temperature(), (temp1+temp2)/2, 6)
 

    def test_pi_temperature(self):
        universe = MMTK.InfiniteUniverse()
        universe.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector([-2.,-2.,-2.])))
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        #Test numerical temperature of a known system
        vbeads = N.array([(1.,1.,1.),(-1.,-1.,-1.), (1.,-1.,1.),(1.,1.,-1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        temp1 = col1.temperature()
        #Twice as many degrees of freedom
        self.assertAlmostEqual(universe[0].mass()/(2.*MMTK.Units.k_B), temp1, 6)
        self.assertAlmostEqual(universe.temperature(), temp1, 6)

        #Test if different velocity atoms result in an averaged temperature
        vbeads = N.array([(1.,1.,1.),(-1.,-1.,-1.), (2.,-2.,2.),(2.,2.,-2.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        temp2 = col1.temperature()
        temp3 = universe.temperature()
        self.assertAlmostEqual(temp3, (temp1+temp2)/2, 6)

        #Test if swapping bead velocities affects the temperature
        vbeads = N.array([(1.,1.,1.),(-2.,-2.,-2.), (2.,-2.,2.),(1.,1.,-1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        self.assertAlmostEqual(temp3, universe.temperature(), 6)

    def test_momentum(self):
        universe = MMTK.InfiniteUniverse()
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', position=MMTK.Vector([0.,0.,0.])))
        col1.addObject(MMTK.Atom('C', position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        #Test if momentum adds up to zero
        vbeads = N.array([(-2.,1.,-1.),(1.,-1.,2.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        self.assert_(N.sum(col1.momentum().array) < 1e-6)

        #Test if the angular momentum adds up to a known value
        am1 = col1.angularMomentum()
        num_am1 = universe[1].mass()*MMTK.Vector([1.,1.,1.]).cross(MMTK.Vector([1.,-1.,2.]))
        self.assertAlmostEqual(am1[0],num_am1[0],5) 
        self.assertAlmostEqual(am1[1],num_am1[1],5) 
        self.assertAlmostEqual(am1[2],num_am1[2],5) 

        #Test is momentum adds up to a calculated value
        universe.atomList()[0].setMass(10.0*MMTK.Units.amu)
        universe.atomList()[1].setMass(1.0*MMTK.Units.amu)
        vbeads = N.array([(2.,2.,2.),(1.,1.,1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        num_p = universe[0].mass()*6. + universe[1].mass()*3.
        self.assert_(N.sum(col1.momentum().array) - num_p < 1e-6)
 
    def test_pi_momentum(self):
        universe = MMTK.InfiniteUniverse()
        col1 = MMTK.Collection() 
        col1.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector([0.,0.,0.])))
        col1.addObject(MMTK.Atom('C', nbeads=2, position=MMTK.Vector([1.,1.,1.])))
        universe.addObject(col1)

        #Test if momentum adds up to zero
        vbeads = N.array([(-1.,1.,-1.), (1.,-1.,1.), (1.,2.,1.), (-2.,-1.,-1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        self.assert_(N.sum(col1.momentum().array) < 1e-6)

        #Test if the angular momentum adds up to a known value
        #Angular momentum for a path integral atom is the average of all bead angular momentums
        am1 = col1.angularMomentum()
        num_am1 = universe[1].mass()*MMTK.Vector([1.,1.,1.]).cross(MMTK.Vector([1.,2.,1.]))/2.
        num_am1 += universe[1].mass()*MMTK.Vector([1.,1.,1.]).cross(MMTK.Vector([-2.,-1.,-1.]))/2.
        self.assertAlmostEqual(am1[0],num_am1[0],5) 
        self.assertAlmostEqual(am1[1],num_am1[1],5) 
        self.assertAlmostEqual(am1[2],num_am1[2],5) 

        #Test is momentum adds up to a calculated value
        universe.atomList()[0].setMass(10.0*MMTK.Units.amu)
        universe.atomList()[1].setMass(1.0*MMTK.Units.amu)
        vbeads = N.array([(2.,2.,2.),(1.,1.,1.),(1.,1.,2.),(2.,2.,1.)])
        v = MMTK.ParticleProperties.ParticleVector(universe, data_array=vbeads)
        universe.setVelocities(v)
        num_p = universe[0].mass()*4.5 + universe[1].mass()*4.5
        self.assert_(N.sum(col1.momentum().array) - num_p < 1e-6)
        
 
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

        print universe.temperature(), col1.temperature()
        self.assertAlmostEqual(num_temp, temp, 6)
 
 
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


class SuperpositionTest(unittest.TestCase):

    """
    Test the findTransformation method
    """

    def test_rotation_translation(self):
        for m in [MMTK.Molecule('water'), Protein('bala1')]:
            universe = MMTK.InfiniteUniverse()
            universe.addObject(m)
            ref_conf = universe.copyConfiguration()
            universe.translateBy(MMTK.Vector(0.1, -1.3, 1.2))
            universe.rotateAroundOrigin(MMTK.Vector(0., 1., 1.), 0.7)
            tr, rms = universe.findTransformation(ref_conf)
            self.assert_(abs(rms) < 1.e-5)
            universe.applyTransformation(tr)
            self.assert_(universe.rmsDifference(ref_conf) < 1.e-5)
            axis, angle = tr.rotation().axisAndAngle()
            self.assert_(abs(angle - 0.7) < 1.e-5)
            self.assert_(abs(abs(N.cos(axis.angle(MMTK.Vector(0., 1., 1.))))-1)
                         < 1.e-5)

def suite():
    loader = unittest.TestLoader()
    s = unittest.TestSuite()
    s.addTest(loader.loadTestsFromTestCase(GroupOfAtomTest))
    s.addTest(loader.loadTestsFromTestCase(SuperpositionTest))
    return s

if __name__ == '__main__':
    unittest.main()
