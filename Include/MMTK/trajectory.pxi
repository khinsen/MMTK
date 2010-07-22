# Include some definitions for Python internals
include 'python.pxi'
# Include the definitions for Numeric arrays
include 'numeric.pxi'

cdef extern from "MMTK/core.h":

    ctypedef double vector3[3]

cdef extern from "MMTK/universe.h": 

    void import_MMTK_universe()

import_MMTK_universe()

cdef extern from "MMTK/trajectory.h":

    void import_MMTK_trajectory()

    cdef enum PyTrajectory_VariableTypes:
        PyTrajectory_Scalar
        PyTrajectory_ParticleScalar
        PyTrajectory_ParticleVector
        PyTrajectory_IntScalar
        PyTrajectory_BoxSize

    cdef enum PyTrajectory_DataClass:
        PyTrajectory_Configuration
        PyTrajectory_Velocities
        PyTrajectory_Gradients
        PyTrajectory_Energy
        PyTrajectory_Thermodynamic
        PyTrajectory_Time
        PyTrajectory_Internal
        PyTrajectory_Auxiliary

    cdef union data:
        int *ip
        double *dp
        PyArrayObject *array

    ctypedef struct PyTrajectoryVariable:
        char *name
        char *text
        char *unit
        data value
        int length
        int type
        #int class
        int modified

import_MMTK_trajectory()

