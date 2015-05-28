# Common functions for PI MD integrators
#

#cython: boundscheck=False, wraparound=False, cdivision=True

from libc.stdint cimport int32_t

cimport numpy as N

from MMTK.core cimport vector3
cimport MMTK_trajectory_generator

from MMTK import Units, ParticleProperties
import MMTK_trajectory_generator

cdef double hbar = Units.hbar

cdef class PIIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

    def __init__(self, universe, options, name):
        MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator \
                   .__init__(self, universe, options, name)

    cdef fixBeadPositions(self, N.ndarray[double, ndim=2] x,
                               Py_ssize_t bead_index, int32_t nb):
        cdef Py_ssize_t i, j
        cdef vector3 *xv = <vector3 *> x.data
        cdef vector3 temp
        if self.universe_spec.is_periodic and nb > 1:
            for j in range(1, nb):
                self.universe_spec.distance_function(temp,
                    xv[bead_index+j-1], xv[bead_index+j],
                    self.universe_spec.geometry_data)
                for i in range(3):
                    xv[bead_index+j][i] = xv[bead_index+j-1][i] + temp[i]

    cdef double springEnergyCartesian(self, N.ndarray[double, ndim=2] x,
                                      N.ndarray[double, ndim=1] m,
                                      N.ndarray[N.int32_t, ndim=2] bd,
                                      double beta):
        cdef Py_ssize_t i, j, k
        cdef int32_t nb
        cdef double sumsq
        cdef double e = 0.
        for i in range(x.shape[0]):
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                sumsq = 0.
                for j in range(nb):
                    for k in range(3):
                        sumsq += (x[i+(j+1)%nb, k]-x[i+j, k])**2
                e += 0.5*nb*nb*m[i]*sumsq/(beta*beta*hbar*hbar)
        return e

    cdef double centroidVirial(self,
                               N.ndarray[double, ndim=2] x,
                               N.ndarray[double, ndim=2] g,
                               N.ndarray[N.int32_t, ndim=2] bd):
        cdef double centroid[3]
        cdef double cvirial = 0.
        cdef Py_ssize_t i, j, k
        cdef int32_t nb
        for i in range(x.shape[0]):
            # bd[i, 0] == 0 means "first bead of an atom"
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                for j in range(3):
                    centroid[j] = 0.
                    for k in range(nb):
                        centroid[j] += x[i+k, j]
                for j in range(3):
                    for k in range(nb):
                        cvirial -= (x[i+k, j]-centroid[j]/nb)*g[i+k, j]
        return cvirial

    cdef void freeze(self, N.ndarray[double, ndim=2] d, N.ndarray[double, ndim=3] ss):
        cdef int ndim = ss.shape[0]
        cdef int npoints = ss.shape[1]
        cdef double dp
        cdef Py_ssize_t i, j, k
        for i in range(ndim):
            dp = 0.
            for j in range(npoints):
                for k in range(3):
                    dp += d[j, k]*ss[i, j, k]
            for j in range(npoints):
                for k in range(3):
                    d[j, k] -= dp*ss[i, j, k]

    cdef int centroidDegreesOfFreedom(self, subspace,
                                      N.ndarray[N.int32_t, ndim=2] bd):
        cdef N.ndarray[double, ndim=2] va
        cdef Py_ssize_t i, j, k
        from MMTK.Subspace import Subspace
        vectors = []
        for k in range(3):
            for i in range(bd.shape[0]):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bd[i, 0] == 0:
                    v = ParticleProperties.ParticleVector(self.universe)
                    vectors.append(v)
                    va = v.array
                    for j in range(bd[i, 1]):
                        va[i+j, k] = 1.
        vectors = [subspace.projectionComplementOf(v) for v in vectors]
        return len(Subspace(self.universe, vectors).getBasis())

