# Trajectory actions in Cython
#
# Written by Konrad Hinsen
#

import MMTK_trajectory

cdef extern from "cobject.h":
    object PyCObject_FromVoidPtr(void *pointer, void (*destruct)(void*))
    void *PyCObject_AsVoidPtr(object)

# The function run_action implements the C interface for trajectory action
# functions, mapping it to the nicer OO interface of the Cython class.
cdef run_action(PyTrajectoryVariable *dynamic_data,
                TrajectoryAction action,
                int step,
                void **scratch,
                object universe):
    # Every trajectory generator has its own scratch pointer,
    # therefore this can be used as an ID for the generator.
    if step == -1:
        return action.initialize(universe, <long>scratch, dynamic_data)
    elif step == -2:
        return action.finalize(universe, <long>scratch, dynamic_data)
    else:
        return action.apply(universe, <long>scratch, step, dynamic_data)

# The base class for trajectoy actions
cdef class TrajectoryAction:

    def __init__(self, long first, object last, long skip):
        self.first = first
        if last is None:
            self.last = MMTK_trajectory.maxint
        else:
            self.last = last
        self.skip = skip

    def getSpecificationList(self, object trajectory_generator, long steps):
        cdef long first
        cdef long last
        first = self.first
        if first < 0:
            first += steps
        last = self.last
        if last < 0:
            last += steps+1
        return ('function', first, last, self.skip,
                PyCObject_FromVoidPtr(<void *>run_action, NULL), self)

    def cleanup(self):
        pass

    # These three methods need to be overridden by subclasses to
    # do any useful work.
    cdef int initialize(self, object universe, long generator_id,
                        PyTrajectoryVariable *dynamic_data):
        return 1

    cdef int finalize(self, object universe, long generator_id,
                      PyTrajectoryVariable *dynamic_data):
        return 1

    cdef int apply(self, object universe, long generator_id,
                   int step, PyTrajectoryVariable *dynamic_data):
        return 1
