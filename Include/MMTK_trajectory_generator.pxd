# Trajectory generators in Cython
#
# Written by Konrad Hinsen
#

cimport numpy as N

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

cdef extern from "stdlib.h":

    ctypedef long size_t
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

cdef class TrajectoryGenerator(object):

    cdef readonly universe
    cdef readonly options, call_options
    cdef readonly features
    cdef readonly actions
    cdef public state_accessor
    cdef PyTrajectoryVariable *tvars
    cdef PyTrajectoryOutputSpec *tspec
    cdef PyUniverseSpecObject *universe_spec
    cdef int natoms, df
    cdef N.ndarray conf_array
    cdef int lock_state

    cdef void declareTrajectoryVariable_double(self, double * var, char *name,
                                               char *text, char*unit,
                                               int data_class) except *
    cdef void declareTrajectoryVariable_int(self, int * var, char *name,
                                            char *text, char*unit,
                                            int data_class) except *
    cdef void declareTrajectoryVariable_array(self, N.ndarray array, char*name,
                                              char *text, char*unit,
                                              int data_class) except *
    cdef void declareTrajectoryVariable_box(self, double * var, int l) except *
    cdef void _addTrajectoryVariable(self, PyTrajectoryVariable v) except *

    cdef void initializeTrajectoryActions(self, char *name) except *
    cdef void finalizeTrajectoryActions(self, int last_step, int error=?) \
              except *
    cdef int trajectoryActions(self, int step) except -1

    cdef void foldCoordinatesIntoBox(self) nogil

    cdef void acquireReadLock(self) nogil
    cdef void releaseReadLock(self) nogil
    cdef void acquireWriteLock(self) nogil
    cdef void releaseWriteLock(self) nogil

    cdef start(self)

cdef class EnergyBasedTrajectoryGenerator(TrajectoryGenerator):

    cdef PyFFEvaluatorObject *evaluator
    cdef evaluator_object

    cdef void calculateEnergies(self, N.ndarray[double, ndim=2] conf_array,
                                energy_data *energy, int small_change=?) \
              except *
