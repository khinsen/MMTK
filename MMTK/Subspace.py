# This module implements subspaces for motion analysis etc.
#
# Written by Konrad Hinsen
# last revision: 2009-2-12
#

"""
Linear subspaces for infinitesimal motions

This module implements subspaces for infinitesimal (or finite
small-amplitude) atomic motions. They can be used in normal mode
calculations or for analyzing complex motions. For an explanation, see:

 -  K. Hinsen and G.R. Kneller
    Projection methods for the analysis of complex motions in macromolecules
    Mol. Sim. 23:275-292 (2000)
"""

__docformat__ = 'epytext'

from MMTK import Utility, ParticleProperties
from Scientific.Geometry import Vector
from Scientific import N

#
# Import LAPACK routines
#
try:
    array_package = N.package
except AttributeError:
    array_package = "Numeric"

dgesdd = None
try:
    if array_package == "Numeric":
        from lapack_lite import dgesdd, LapackError
    else:
        from numpy.linalg.lapack_lite import dgesdd, LapackError
except ImportError: pass
if dgesdd is None:
    try:
        # PyLAPACK
        from lapack_dge import dgesdd
    except ImportError: pass
if dgesdd:
    n = 1
    array = N.zeros((n, n), N.Float)
    sv = N.zeros((n,), N.Float)
    u = N.zeros((n, n), N.Float)
    vt = N.zeros((n, n), N.Float)
    work = N.zeros((1,), N.Float)
    int_types = [N.Int, N.Int8, N.Int16, N.Int32]
    try:
        int_types.append(N.Int64)
        int_types.append(N.Int128)
    except AttributeError:
        pass
    for int_type in int_types:
        iwork = N.zeros((1,), int_type)
        try:
            dgesdd('S', n, n, array, n, sv, u, n, vt, n, work, -1, iwork, 0)
            break
        except LapackError:
            pass
    del n, array, sv, u, vt, work, iwork, int_types

del array_package

#
# A set of particle vectors that define a subspace
#
class ParticleVectorSet(object):

    def __init__(self, universe, data):
        self.universe = universe
        if type(data) == N.arraytype:
            self.n = data.shape[0]
            self.array = data
        else:
            self.n = data
            self.array = N.zeros((self.n, universe.numberOfAtoms(), 3),
                                 N.Float)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        if i >= self.n:
            raise IndexError
        return ParticleProperties.ParticleVector(self.universe, self.array[i])


class Subspace(object):

    """
    Subspace of infinitesimal atomic motions
    """

    def __init__(self, universe, vectors):
        """
        @param universe: the universe for which the subspace is created
        @type universe: L{MMTK.Universe.Universe}
        @param vectors: a list of L{MMTK.ParticleProperties.ParticleVector}
                        objects that define the subspace. They need not be
                        orthogonal or linearly independent.
        @type vectors: C{list}
        """
        self.universe = universe
        self.vectors = vectors
        self._basis = None

    def __add__(self, other):
        return Subspace(self.vectors+other.vectors)

    def __len__(self):
        return len(self.vectors)

    def __getitem__(self, item):
        return self.vectors[item]

    def getBasis(self):
        """
        Construct a basis for the subspace by orthonormalization of
        the input vectors using Singular Value Decomposition. The
        basis consists of a sequence of
        L{MMTK.ParticleProperties.ParticleVector}
        objects that are orthonormal in configuration space.
        @returns: the basis
        """
        if self._basis is None:
            basis = N.array([v.array for v in self.vectors], N.Float)
            shape = basis.shape
            nvectors = shape[0]
            natoms = shape[1]
            basis.shape = (nvectors, 3*natoms)
            sv = N.zeros((min(nvectors, 3*natoms),), N.Float)
            min_n_m = min(3*natoms, nvectors)
            vt = N.zeros((nvectors, min_n_m), N.Float)
            work = N.zeros((1,), N.Float)
            iwork = N.zeros((8*min_n_m,), int_type)
            if 3*natoms >= nvectors:
                result = dgesdd('O', 3*natoms, nvectors, basis, 3*natoms,
                                sv, basis, 3*natoms, vt, min_n_m,
                                work, -1, iwork, 0)
                work = N.zeros((int(work[0]),), N.Float)
                result = dgesdd('O', 3*natoms, nvectors, basis, 3*natoms,
                                sv, basis, 3*natoms, vt, min_n_m,
                                work, work.shape[0], iwork, 0)
                u = basis
            else:
                u = N.zeros((min_n_m, 3*natoms), N.Float)
                result = dgesdd('S', 3*natoms, nvectors, basis, 3*natoms,
                                sv, u, 3*natoms, vt, min_n_m,
                                work, -1, iwork, 0)
                work = N.zeros((int(work[0]),), N.Float)
                result = dgesdd('S', 3*natoms, nvectors, basis, 3*natoms,
                                sv, u, 3*natoms, vt, min_n_m,
                                work, work.shape[0], iwork, 0)
            if result['info'] != 0:
                raise ValueError('Lapack SVD: ' + `result['info']`)
            svmax = N.maximum.reduce(sv)
            nvectors = N.add.reduce(N.greater(sv, 1.e-10*svmax))
            u = u[:nvectors]
            u.shape = (nvectors, natoms, 3)
            self._basis = ParticleVectorSet(self.universe, u)
        return self._basis

    def projectionOf(self, vector):
        """
        @param vector: a particle vector
        @type vector: L{MMTK.ParticleProperties.ParticleVector}
        @returns: the projection of the vector onto the subspace.
        """
        vector = vector.array
        basis = self.getBasis().array
        p = N.zeros(vector.shape, N.Float)
        for bv in basis:
            N.add(N.add.reduce(N.ravel(bv*vector))*bv, p, p)
        return ParticleProperties.ParticleVector(self.universe, p)

    def projectionComplementOf(self, vector):
        """
        @param vector: a particle vector
        @type vector: L{MMTK.ParticleProperties.ParticleVector}
        @returns: the projection of the vector onto the orthogonal complement
                  of the subspace.
        """
        return vector - self.projectionOf(vector)


class RigidMotionSubspace(Subspace):

    """
    Subspace of rigid-body motions

    A rigid-body motion subspace is the subspace which contains
    the rigid-body motions of any number of chemical objects.
    """

    def __init__(self, universe, objects):
        """
        @param universe: the universe for which the subspace is created
        @type universe: L{MMTK.Universe.Universe}
        @param objects: a sequence of objects whose rigid-body motion is
                        included in the subspace
        """
        if not Utility.isSequenceObject(objects):
            objects = [objects]
        vectors = []
        for o in objects:
            atoms = o.atomList()
            for d in [Vector(1.,0.,0.), Vector(0.,1.,0.), Vector(0.,0.,1.)]:
                v = ParticleProperties.ParticleVector(universe)
                for a in atoms:
                    v[a] = d
                vectors.append(v/N.sqrt(len(atoms)))
            if len(atoms) > 1:
                center = o.centerOfMass()
                iv = len(vectors)-3
                for d in [Vector(1.,0.,0.),Vector(0.,1.,0.),Vector(0.,0.,1.)]:
                    v = ParticleProperties.ParticleVector(universe)
                    for a in atoms:
                        v[a] = d.cross(a.position()-center)
                    for vt in vectors[iv:]:
                        v = v - v.dotProduct(vt)*vt
                    vectors.append(v/N.sqrt(v.dotProduct(v)))
        Subspace.__init__(self, universe, vectors)
        # The vector set is already orthonormal by construction
        # (assuming that the rigid bodies have no atoms in common),
        # so we can eliminate the lengthy SVD procedure
        count = ParticleProperties.ParticleScalar(universe)
        for o in objects:
            count = count + o.booleanMask()
        if N.maximum.reduce(count.array) == 1:
            self._basis = ParticleVectorSet(universe, len(vectors))
            for i in range(len(vectors)):
                self._basis.array[i] = vectors[i].array


class PairDistanceSubspace(Subspace):

    """
    Subspace of pair-distance motions

    A pair-distance motion subspace is the subspace which contains
    the relative motions of any number of atom pairs along
    their distance vector, e.g. bond elongation between two
    bonded atoms.
    """

    def __init__(self, universe, atom_pairs):
        """
        @param universe: the universe for which the subspace is created
        @type universe: L{MMTK.Universe.Universe}
        @param atom_pairs: a sequence of atom pairs whose distance-vector
                           motion is included in the subspace
        """
        vectors = []
        for atom1, atom2 in atom_pairs:
            v = ParticleProperties.ParticleVector(universe)
            distance = atom1.position()-atom2.position()
            v[atom1] = distance
            v[atom2] = -distance
            vectors.append(v)
        Subspace.__init__(self, universe, vectors)
