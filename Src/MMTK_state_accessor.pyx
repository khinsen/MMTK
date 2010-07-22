# A special trajectory action that provides access to the state of
# any trajectory generator
#
# Written by Konrad Hinsen
# last revision: 2010-1-22
#

include 'MMTK/trajectory.pxi'

cimport MMTK_trajectory_action

import MMTK
import MMTK_trajectory

cdef class StateAccessor(MMTK_trajectory_action.TrajectoryAction):

    cdef long generator_id
    cdef PyTrajectoryVariable *dynamic_data
    cdef object universe
    cdef bint accessing_state
    
    def __init__(self):
        MMTK_trajectory_action.TrajectoryAction.__init__(self,
                                   0, None, MMTK_trajectory.maxint)

    def __cinit__(self):
        # Cython zeros all of these automatically, but...
        # ... explicit is better than implicit, says The Zen of Python.
        self.generator_id = 0
        self.dynamic_data = NULL
        self.universe = None
        self.accessing_state = False

    def copyState(self):
        """
        Makes a copy of the current state of a trajectory generator
        in a thread-safe way.
        @returns: a copy of the current state, or C{None} if no trajectory
                  generator is currently active
        @rtype:  C{dict}
        """
        cdef PyTrajectoryVariable *dynamic_data = self.dynamic_data
        self.universe.acquireReadStateLock()
        state = {}
        # Timing can be tricky. The trajectory generator runs in another
        # thread and may well terminate while copyState is executed.
        # copyState sets self.accessing_state while it reads data, which
        # prevents finalize from returning. However, it is still possible
        # that finalize sets self.generator_id to 0 and returns just before
        # self.accesing_state is set to 1. Therefore the test for
        # termination must be where it is (and not before).
        self.accessing_state = True
        if self.generator_id == 0:
            self.accessing_state = False
            self.universe.releaseReadStateLock()
            return None
        while dynamic_data.name != NULL:
            name = dynamic_data.name
            dtype = dynamic_data.type
            if name == 'configuration':
                state[name] = MMTK.copy(self.universe.configuration())
            elif name == 'box_size':
                continue
            else:
                if dtype == PyTrajectory_Scalar:
                    state[name] = dynamic_data.value.dp[0]
                elif dtype == PyTrajectory_IntScalar:
                    state[name] = dynamic_data.value.ip[0]
                else:
                    array = <object>dynamic_data.value.array
                    array = array.copy()
                    if dtype == PyTrajectory_ParticleScalar:
                        state[name] = MMTK.ParticleScalar(self.universe, array)
                    elif dtype == PyTrajectory_ParticleVector:
                        state[name] = MMTK.ParticleVector(self.universe, array)
            dynamic_data += 1
        self.accessing_state = False
        self.universe.releaseReadStateLock()
        return state
    
    cdef int initialize(self, object universe, long generator_id,
                        PyTrajectoryVariable *dynamic_data):
        self.generator_id = generator_id
        self.universe = universe
        self.dynamic_data = dynamic_data
        return 1

    cdef int finalize(self, object universe, long generator_id,
                      PyTrajectoryVariable *dynamic_data):
        assert self.generator_id == generator_id
        assert self.universe == universe
        assert self.dynamic_data == dynamic_data
        while self.accessing_state:
            continue
        self.generator_id = 0
        return 1
