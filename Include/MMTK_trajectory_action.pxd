# Trajectory actions in Cython
#
# Written by Konrad Hinsen
#

include 'MMTK/trajectory.pxi'

cdef class TrajectoryAction:

    cdef long first
    cdef long last
    cdef long skip

    cdef int initialize(self, object universe, long generator_id,
                        PyTrajectoryVariable *dynamic_data)

    cdef int finalize(self, object universe, long generator_id,
                      PyTrajectoryVariable *dynamic_data)

    cdef int apply(self, object universe, long generator_id,
                   int step, PyTrajectoryVariable *dynamic_data)
