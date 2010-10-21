# Trajectory generators in Cython
#
# Written by Konrad Hinsen
#

import_MMTK_universe()
import_MMTK_trajectory()
import_MMTK_forcefield()

from MMTK import Features
import numpy as N
cimport numpy as N
import cython

cdef extern from "stdlib.h":

    ctypedef long size_t
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

#
# Base class for trajectory generators
#
cdef class TrajectoryGenerator(object):

    """
    Trajectory generator base class

    This base class implements the common aspects of everything that
    generates trajectories: integrators, minimizers, etc.
    """

    def __init__(self, universe, options, name):
        self.universe = universe
        self.options = options
        self.name = name
        self.call_options = {}
        self.features = []
        self.tvars = NULL
        self.tspec = NULL
        
    def setCallOptions(self, options):
        self.call_options = options

    def getActions(self):
        try:
            steps = self.getOption('steps')
        except ValueError:
            steps = None
        return [a.getSpecificationList(self, steps) for a in self.actions]

    def cleanupActions(self):
        for a in self.actions:
            a.cleanup()

    def getOption(self, option):
        try:
            value = self.call_options[option]
        except KeyError:
            try:
                value = self.options[option]
            except KeyError:
                try:
                    value = self.default_options[option]
                except KeyError:
                    raise ValueError('undefined option: ' + option)
        return value

    def optionString(self, options):
        return sum((o + '=' + `self.getOption(o)` + ', ' for o in options),
                   '')[:-2]

    def __call__(self, **options):
        self.setCallOptions(options)
        try:
            self.actions = self.getOption('actions')
        except ValueError:
            self.actions = []
        try:
            if self.getOption('background'):
                import MMTK_state_accessor
                self.state_accessor = MMTK_state_accessor.StateAccessor()
                self.actions.append(self.state_accessor)
        except ValueError:
            pass
        Features.checkFeatures(self, self.universe)
        if self.tvars != NULL:
            free(self.tvars)
            self.tvars = NULL
        configuration = self.universe.configuration()
        self.conf_array = configuration.array
        self.declareTrajectoryVariable_array(self.conf_array,
                                             "configuration",
                                             "Configuration:\n",
                                             length_unit_name,
                                             PyTrajectory_Configuration)
        self.universe_spec = <PyUniverseSpecObject *>self.universe._spec
        if self.universe_spec.geometry_data_length > 0:
            self.declareTrajectoryVariable_box(
                self.universe_spec.geometry_data,
                self.universe_spec.geometry_data_length)
        masses = self.universe.masses()
        self.declareTrajectoryVariable_array(masses.array,
                                             "masses",
                                             "Masses:\n",
                                             mass_unit_name,
                                             PyTrajectory_Internal)
        self.natoms = self.universe.numberOfAtoms()
        self.df = self.universe.degreesOfFreedom()
        self.declareTrajectoryVariable_int(&self.df,
                                           "degrees_of_freedom",
                                           "Degrees of freedom: %d\n",
                                           "", PyTrajectory_Internal)
        if self.getOption('background'):
            from MMTK import ThreadManager
            return ThreadManager.TrajectoryGeneratorThread(
                self.universe, self.start_py, (), self.state_accessor)
        else:
            self.start()
            return None

    def start_py(self):
        self.start()

    cdef start(self):
        pass

    cdef void _addTrajectoryVariable(self, PyTrajectoryVariable v) except *:
        cdef PyTrajectoryVariable *tv
        cdef int i, n
        if self.tvars == NULL:
            self.tvars = <PyTrajectoryVariable *> \
                         malloc(2*sizeof(PyTrajectoryVariable))
            if self.tvars == NULL:
                raise MemoryError
            self.tvars[0] = v
            self.tvars[1].name = NULL
        else:
            n = 0
            tv = self.tvars
            while tv.name != NULL:
                n += 1
                tv += 1
            tv = <PyTrajectoryVariable *> \
                 malloc((n+2)*sizeof(PyTrajectoryVariable))
            if tv == NULL:
                raise MemoryError
            for i in range(n):
                tv[i] = self.tvars[i]
            tv[n] = v
            tv[n+1].name = NULL
            free(self.tvars)
            self.tvars = tv

    cdef void declareTrajectoryVariable_double(self, double * var, char *name,
                                               char *text, char*unit,
                                               int data_class) except *:
        cdef PyTrajectoryVariable v
        v.name = name
        v.text = text
        v.unit = unit
        v.type = PyTrajectory_Scalar
        v.data_class = data_class
        v.value.dp = var
        self._addTrajectoryVariable(v)

    cdef void declareTrajectoryVariable_int(self, int * var, char *name,
                                            char *text, char*unit,
                                            int data_class) except *:
        cdef PyTrajectoryVariable v
        v.name = name
        v.text = text
        v.unit = unit
        v.type = PyTrajectory_IntScalar
        v.data_class = data_class
        v.value.ip = var
        self._addTrajectoryVariable(v)

    cdef void declareTrajectoryVariable_array(self, N.ndarray array, char*name,
                                              char *text, char*unit,
                                              int data_class) except *:
        cdef PyTrajectoryVariable v
        v.name = name
        v.text = text
        v.unit = unit
        v.data_class = data_class
        v.value.array = <PyArrayObject *>array
        if array.ndim == 1:
            v.type = PyTrajectory_ParticleScalar
        elif array.ndim == 2:
            v.type = PyTrajectory_ParticleVector
        self._addTrajectoryVariable(v)

    cdef void declareTrajectoryVariable_box(self, double * var, int l) except *:
        cdef PyTrajectoryVariable v
        v.name = "box_size"
        v.text = "Box size:"
        v.unit = length_unit_name
        v.type = PyTrajectory_BoxSize
        v.data_class = PyTrajectory_Configuration
        v.value.dp = var
        v.length = l
        self._addTrajectoryVariable(v)

    cdef void initializeTrajectoryActions(self) except *:
        self.tspec = PyTrajectory_OutputSpecification(self.universe,
                                                      self.getActions(),
                                                      self.name, self.tvars);
        if self.tspec == NULL:
            raise MemoryError

    cdef void finalizeTrajectoryActions(self, int last_step,
                                        int error=False) except *:
        if error:
            PyTrajectory_OutputFinish(self.tspec, last_step, 1, 1, self.tvars)
        else:
            PyTrajectory_OutputFinish(self.tspec, last_step, 0, 1, self.tvars)

    cdef int trajectoryActions(self, int step) except -1:
        return PyTrajectory_Output(self.tspec, step, self.tvars, NULL)

    cdef void foldCoordinatesIntoBox(self) nogil:
        self.universe_spec.correction_function(<vector3 *>self.conf_array.data,
                                               self.natoms,
                                               self.universe_spec.geometry_data)

    cdef void acquireReadLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, 1)
        self.lock_state = 1

    cdef void releaseReadLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, 2)
        self.lock_state = 0

    cdef void acquireWriteLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, -1)
        self.lock_state = -1

    cdef void releaseWriteLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, -2)
        self.lock_state = 0

#
# Base class for trajectory generators that call the C-level
# energy evaluators. It implements a mechanism that makes such
# generators thread-safe.
#
cdef class EnergyBasedTrajectoryGenerator(TrajectoryGenerator):

    def __init__(self, universe, options, name):
        TrajectoryGenerator.__init__(self, universe, options, name)
        evaluator_object = None
        self.evaluator = NULL

    cdef void initializeTrajectoryActions(self) except *:
        TrajectoryGenerator.initializeTrajectoryActions(self)
        # Construct a C evaluator object for the force field, using
        # the specified number of threads or the default value
        nt = self.getOption('threads')
        self.evaluator_object = \
                self.universe.energyEvaluator(threads=nt).CEvaluator()
        self.evaluator = <PyFFEvaluatorObject*>self.evaluator_object

    cdef void finalizeTrajectoryActions(self, int last_step,
                                        int error=False) except *:
        TrajectoryGenerator.finalizeTrajectoryActions(self, last_step, error)
        self.evaluator = NULL
        self.evaluator_object = None
    
    cdef int trajectoryActions(self, int step) except -1:
        cdef int ret_code
        self.evaluator.tstate_save = PyEval_SaveThread()
        if self.lock_state != 0:
            PyUniverseSpec_StateLock(self.universe_spec, 2*self.lock_state)
        ret_code = PyTrajectory_Output(self.tspec, step, self.tvars,
                                       &self.evaluator.tstate_save)
        if self.lock_state != 0:
            PyUniverseSpec_StateLock(self.universe_spec, self.lock_state)
        PyEval_RestoreThread(self.evaluator.tstate_save)
        return ret_code

    cdef void calculateEnergies(self, N.ndarray[double, ndim=2] conf_array,
                                energy_data *energy, int small_change=0) \
                            except *:
        self.evaluator.tstate_save = PyEval_SaveThread()
        if self.lock_state == -1:
            PyUniverseSpec_StateLock(self.universe_spec, -2)
        if self.lock_state != 1:
            PyUniverseSpec_StateLock(self.universe_spec, 1)
        self.evaluator.eval_func(self.evaluator, energy, <PyArrayObject *>conf_array,
                                 small_change)
        if self.lock_state != 1:
            PyUniverseSpec_StateLock(self.universe_spec, 2)
        if self.lock_state == -1:
            PyUniverseSpec_StateLock(self.universe_spec, -1)
        PyEval_RestoreThread(self.evaluator.tstate_save)
