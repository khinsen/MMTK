# Common functions for PI MD integrators
#

cimport numpy as N

cimport MMTK_trajectory_generator

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

cdef extern from "stdlib.h":
    cdef double fabs(double)
    cdef double sqrt(double)
    cdef double sin(double)
    cdef double cos(double)
    cdef double exp(double)
    cdef double M_PI


cdef class PIIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

    cdef fixBeadPositions(self, N.ndarray[double, ndim=2] x,
                               int bead_index, int nb)

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

