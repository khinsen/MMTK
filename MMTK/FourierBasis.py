# Fourier basis for low-frequency normal mode calculations.
#
# Written by Konrad Hinsen
#

"""
Fourier basis for low-frequency normal mode calculations

This module provides a basis that is suitable for the
calculation of low-frequency normal modes. The basis is
derived from vector fields whose components are stationary
waves in a box surrounding the system. For a description
see
     K. Hinsen
     Analysis of domain motions by approximate normal mode calculations
     Proteins 33 (1998): 417-429
"""

from MMTK import ParticleProperties
from Scientific.Geometry import Vector
from Scientific import N

class FourierBasis(object):

    """
    Collective-motion basis for low-frequency normal mode calculations

    A FourierBasis behaves like a sequence of
    L{MMTK.ParticleProperties.ParticleVector} objects. The vectors are
    B{not} orthonormal, because orthonormalization is handled
    automatically by the normal mode class.
    """

    def __init__(self, universe, cutoff):
        """
        @param universe: the universe for which the basis will be used
        @type universe: L{MMTK.Universe.Universe}
        @param cutoff: the wavelength cutoff. A smaller value yields
                       a larger basis.
        @type cutoff: C{float}
        """
        p1, p2 = universe.boundingBox()
        p2 = p2 + Vector(cutoff, cutoff, cutoff)
        l = (p2-p1).array
        n_max = (0.5*l/cutoff+0.5).astype(N.Int)

        wave_numbers = [(nx, ny, nz)
                        for nx in range(-n_max[0], n_max[0]+1)
                        for ny in range(-n_max[1], n_max[1]+1)
                        for nz in range(-n_max[2], n_max[2]+1)
                        if (nx/l[0])**2 + (ny/l[1])**2 + (nz/l[2])**2 \
                                    < 0.25/cutoff**2]

        atoms = universe.atomList()
        natoms = len(atoms)
        basis = N.zeros((3*len(wave_numbers)+3, natoms, 3), N.Float)
        cm = universe.centerOfMass()
        i = 0
        for rotation in [Vector(1.,0.,0.), Vector(0.,1.,0.), Vector(0.,0.,1.)]:
            v = ParticleProperties.ParticleVector(universe, basis[i])
            for a in atoms:
                v[a] = rotation.cross(a.position()-cm)
            i += i
        conf = universe.configuration().array-p1.array
        for n in wave_numbers:
            k = 2.*N.pi*N.array(n)/l
            w = self._w(conf[:, 0], k[0]) * self._w(conf[:, 1], k[1]) * \
                self._w(conf[:, 2], k[2])
            basis[i, :, 0] = w
            basis[i+1, :, 1] = w
            basis[i+2, :, 2] = w
            i += 3

        self.array = basis
        self.universe = universe

    __safe_for_unpickling__ = True
    __had_initargs__ = True

    def _w(self, x, k):
        if k < 0:
            return N.sin(-k*x)
        else:
            return N.cos(k*x)

    def __len__(self):
        return self.array.shape[0]

    def __getitem__(self, item):
        return ParticleProperties.ParticleVector(self.universe,
                                                 self.array[item])


# Estimate number of basis vectors for a given cutoff

def countBasisVectors(universe, cutoff):
    """
    Estimate the number of basis vectors for a given universe and cutoff
    @param universe: the universe
    @type universe: L{MMTK.Universe.Universe}
    @param cutoff: the wavelength cutoff. A smaller value yields a larger basis.
    @type cutoff: C{float}
    @returns: the number of basis vectors in a FourierBasis
    @rtype: C{int}
    """
    p1, p2 = universe.boundingBox()
    p2 = p2 + Vector(cutoff, cutoff, cutoff)
    l = (p2-p1).array
    n_max = (0.5*l/cutoff+0.5).astype(N.Int)
    n_wave_numbers = 0
    for nx in range(-n_max[0], n_max[0]+1):
        for ny in range(-n_max[1], n_max[1]+1):
            for nz in range(-n_max[2], n_max[2]+1):
                if (nx/l[0])**2 + (ny/l[1])**2 + (nz/l[2])**2 < 0.25/cutoff**2:
                    n_wave_numbers += 1
    return 3*n_wave_numbers+3


# Estimate cutoff for a given number of basis vectors

def estimateCutoff(universe, nmodes):
    """
    Estimate the cutoff that yields a given number of basis vectors
    for a given universe.
    @param universe: the universe
    @type universe: L{MMTK.Universe.Universe}
    @param nmodes: the number of basis vectors in a FourierBasis
    @type nmodes: C{int}
    @returns: the wavelength cutoff and the precise number of basis vectors
    @rtype: C{tuple} C{(float, int)}
    """
    natoms = universe.numberOfCartesianCoordinates()
    if nmodes > natoms:
        nmodes = 3*natoms
        cutoff = None
    else:
        p1, p2 = universe.boundingBox()
        cutoff_max = (p2-p1).length()
        cutoff = 0.5*cutoff_max
        nmodes_opt = nmodes
        nmodes = countBasisVectors(universe, cutoff)
        while nmodes > nmodes_opt:
            cutoff += 0.1
            if cutoff > cutoff_max:
                cutoff = cutoff_max
                break
            nmodes = countBasisVectors(universe, cutoff)
        while nmodes < nmodes_opt:
            cutoff -= 0.1
            if cutoff < 0.1:
                cutoff = 0.1
                break
            nmodes = countBasisVectors(universe, cutoff)
    return cutoff, nmodes
