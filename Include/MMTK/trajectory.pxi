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
        int data_class "class"
        int modified

    cdef char *length_unit_name
    cdef char *volume_unit_name
    cdef char *time_unit_name
    cdef char *frequency_unit_name
    cdef char *frequency_square_unit_name
    cdef char *velocity_unit_name
    cdef char *mass_unit_name
    cdef char *energy_unit_name
    cdef char *energy_gradient_unit_name
    cdef char *temperature_unit_name
    cdef char *pressure_unit_name

    ctypedef struct PyTrajectoryOutputSpec

    cdef PyTrajectoryOutputSpec *PyTrajectory_OutputSpecification(object universe,
                                                                  object spec_list,
                                                                  char *description,
                                                                  PyTrajectoryVariable *data)
    cdef void PyTrajectory_OutputFinish(PyTrajectoryOutputSpec *spec, int step, int error_flag,
                                        int time_stamp_flag, PyTrajectoryVariable *data)
    cdef int PyTrajectory_Output(PyTrajectoryOutputSpec *spec, int step,
                                 PyTrajectoryVariable *data, PyThreadState **thread)

import_MMTK_trajectory()

