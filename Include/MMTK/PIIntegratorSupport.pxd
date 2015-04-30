# Common functions for PI MD integrators
#

from libc.stdint cimport int32_t

cimport numpy as N

cimport MMTK_trajectory_generator


cdef class PIIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

    cdef fixBeadPositions(self, N.ndarray[double, ndim=2] x,
                               Py_ssize_t bead_index, int32_t nb)

    cdef double springEnergyCartesian(self, N.ndarray[double, ndim=2] x,
                                      N.ndarray[double, ndim=1] m,
                                      N.ndarray[N.int32_t, ndim=2] bd,
                                      double beta)

    cdef double centroidVirial(self,
                               N.ndarray[double, ndim=2] x,
                               N.ndarray[double, ndim=2] g,
                               N.ndarray[N.int32_t, ndim=2] bd)

    cdef void freeze(self, N.ndarray[double, ndim=2] d, N.ndarray[double, ndim=3] ss)

    cdef int centroidDegreesOfFreedom(self, subspace,
                                      N.ndarray[N.int32_t, ndim=2] bd)

