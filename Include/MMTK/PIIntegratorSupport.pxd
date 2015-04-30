# Common functions for PI MD integrators
#

cimport numpy as N

cimport MMTK_trajectory_generator


cdef class PIIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

    cdef fixBeadPositions(self, N.ndarray[double, ndim=2] x,
                               Py_ssize_t bead_index, int nb)

    cdef double springEnergyCartesian(self, N.ndarray[double, ndim=2] x,
                                      N.ndarray[double, ndim=1] m,
                                      N.ndarray[short, ndim=2] bd,
                                      double beta)

    cdef double centroidVirial(self,
                               N.ndarray[double, ndim=2] x,
                               N.ndarray[double, ndim=2] g,
                               N.ndarray[short, ndim=2] bd)

    cdef void freeze(self, N.ndarray[double, ndim=2] d, N.ndarray[double, ndim=3] ss)

    cdef int centroidDegreesOfFreedom(self, subspace,
                                      N.ndarray[short, ndim=2] bd)

