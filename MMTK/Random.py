# Functions for finding random points and orientations.
#
# Written by: Konrad Hinsen
# Last revision: 2009-4-29
# 

"""
Random quantities for use in molecular simulations
"""

__docformat__ = 'epytext'

from Scientific.Geometry import Vector
from Scientific.Geometry.Transformation import Rotation
from Scientific import N
from MMTK import ParticleProperties, Units

try:
    numeric = N.package
except AttributeError:
    numeric = "Numeric"

RNG = None
try:
    if numeric == "Numeric":
        import RNG
    elif numeric == "NumPy":
        import numpy.oldnumeric.rng as RNG
except ImportError:
    pass


if RNG is None:

    if numeric == "Numeric":
        from RandomArray import uniform, seed
    elif numeric == "NumPy":
        from numpy.oldnumeric.random_array import uniform, seed
    random = __import__('random')
    seed(1, 1)
    random.seed(1)

    def initializeRandomNumbersFromTime():
        random.seed()
        seed(0, 0)

    def gaussian(mean, std, shape=None):
        if shape is None:
            x = random.normalvariate(0., 1.)
        else:
            x = N.zeros(shape, N.Float)
            xflat = N.ravel(x)
            for i in range(len(xflat)):
                xflat[i] = random.normalvariate(0., 1.)
        return mean + std*x

else:

    _uniform_generator = \
                    RNG.CreateGenerator(-1, RNG.UniformDistribution(0., 1.))
    _gaussian_generator = \
                    RNG.CreateGenerator(-1, RNG.NormalDistribution(0., 1.))

    def initializeRandomNumbersFromTime():
        global _uniform_generator, _gaussian_generator
        _uniform_generator = \
                       RNG.CreateGenerator(0, RNG.UniformDistribution(0., 1.))
        _gaussian_generator = \
                       RNG.CreateGenerator(0, RNG.NormalDistribution(0., 1.))

    def uniform(x1, x2, shape=None):
        if shape is None:
            x = _uniform_generator.ranf()
        else:
            n = N.multiply.reduce(shape)
            x = _uniform_generator.sample(n)
            x.shape = shape
        return x1+(x2-x1)*x

    def gaussian(mean, std, shape=None):
        if shape is None:
            x = _gaussian_generator.ranf()
        else:
            n = N.multiply.reduce(shape)
            x = _gaussian_generator.sample(n)
            x.shape = shape
        return mean+std*x

del numeric

#
# Random point in a rectangular box centered around the origin
#
def randomPointInBox(a, b = None, c = None):
    """
    @param a: the edge length of a box along the x axis
    @type a: C{float}
    @param b: the edge length of a box along the y axis (default: a)
    @type b: C{float}
    @param c: the edge length of a box along the z axis (default: a)
    @type c: C{float}
    @returns: a vector drawn from a uniform distribution within a
              rectangular box with edge lengths a, b, c.
    @rtype: C{Scientific.Geometry.Vector}
    """
    if b is None: b = a
    if c is None: c = a
    x = uniform(-0.5*a, 0.5*a)
    y = uniform(-0.5*b, 0.5*b)
    z = uniform(-0.5*c, 0.5*c)
    return Vector(x, y, z)

#
# Random point in a sphere around the origin.
#
def randomPointInSphere(r):
    """
    @param r: the radius of a sphere
    @type r: C{float}
    @returns: a vector drawn from a uniform distribution within
              a sphere of radius r.
    @rtype: C{Scientific.Geometry.Vector}
    """
    rsq = r*r
    while 1:
        x = N.array([uniform(-r, r), uniform(-r, r), uniform(-r, r)])
        if N.dot(x, x) < rsq: break
    return Vector(x)

#
# Random direction (unit vector).
#
def randomDirection():
    """
    @returns: a vector drawn from a uniform distribution on the surface
              of a unit sphere.
    @rtype: C{Scientific.Geometry.Vector}
    """
    r = randomPointInSphere(1.)
    return r.normal()

def randomDirections(n):
    """
    @param n: the number of direction vectors
    @returns: a list of n vectors drawn from a uniform distribution on
              the surface of a unit sphere. If n is negative, returns
              a deterministic list of not more than -n vectors of unit
              length (useful for testing purposes).
    @rtype: C{list}
    """
    if n < 0:
        vs = [Vector(1., 0., 0.), Vector(0., -1., 0.), Vector(0., 0., 1.),
              Vector(-1., 1., 0.).normal(), Vector(-1., 0., 1.).normal(),
              Vector(0., 1., -1.).normal(), Vector(1., -1., 1.).normal()]
        return vs[:-n]
    else:
        return [randomDirection() for i in range(n)]

#
# Random rotation.
#
def randomRotation(max_angle = N.pi):
    """
    @param max_angle: the upper limit for the rotation angle
    @type max_angle: C{float}
    @returns: a random rotation with a uniform axis distribution
              and angles drawn from a uniform distribution between
              -max_angle and max_angle.
    @rtype: C{Scientific.Geometry.Transformations.Rotation}
    """
    return Rotation(randomDirection(), uniform(-max_angle, max_angle))

#
# Random velocity (gaussian)
#
def randomVelocity(temperature, mass):
    """
    @param temperature: the temperature defining a Maxwell distribution
    @type temperature: C{float}
    @param mass: the mass of a particle
    @type mass: C{float}
    @returns: a random velocity vector for a particle of a given mass
              at a given temperature
    @rtype: C{Scientific.Geometry.Vector}
    """
    sigma = N.sqrt((temperature*Units.k_B)/(mass*Units.amu))
    return Vector(gaussian(0., sigma, (3,)))

#
# Random ParticleVector (gaussian)
#
def randomParticleVector(universe, width):
    """
    @param universe: a universe
    @type universe: L{MMTK.Universe.Universe}
    @param width: the width (standard deviation) of a Gaussian distribution
    @type width: C{float}
    @returns: a set of vectors drawn from a Gaussian distribution
              with a given width centered  around zero.
    @rtype: L{MMTK.ParticleProperties.ParticleVector}
    """
    data = gaussian(0., 0.577350269189*width, (universe.numberOfPoints(), 3))
    return ParticleProperties.ParticleVector(universe, data)
