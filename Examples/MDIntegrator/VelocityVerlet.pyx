# This module implements a Velocity Verlet integrator.
#
# Written by Konrad Hinsen
#

from MMTK import Dynamics, Environment, Features, Trajectory, \
                 Units, ParticleProperties
import numpy as N
cimport numpy as N
import cython

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

#
# Integrator base class
#
cdef class TrajectoryGenerator(object):

    cdef public universe
    cdef public options, call_options
    cdef public features
    cdef readonly actions
    cdef PyTrajectoryVariable *tvars
    cdef PyTrajectoryOutputSpec *tspec
    cdef PyUniverseSpecObject *universe_spec
    cdef int natoms, df
    cdef N.ndarray conf_array
    cdef int lock_state

    """
    Trajectory generator base class

    This base class implements the common aspects of everything that
    generates trajectories: integrators, minimizers, etc.
    """

    def __init__(self, universe, options):
        self.universe = universe
        self.options = options
        self.call_options = {}
        self.features = []
        self.tvars = NULL
        self.tspec = NULL
        
    def setCallOptions(self, options):
        self.call_options = options

    def getActions(self):
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
        s = ''
        for o in options:
            s = s + o + '=' + `self.getOption(o)` + ', '
        return s[:-2]

    def run(self, function, args):
        if self.getOption('background'):
            import ThreadManager
            return ThreadManager.TrajectoryGeneratorThread(self.universe,
                                      function, args, self.state_accessor)
        else:
            apply(function, args)

    def __call__(self, **options):
        self.setCallOptions(options)
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

    cdef void initializeTrajectoryActions(self, char *name) except *:
        self.tspec = PyTrajectory_OutputSpecification(self.universe,
                                                      self.getActions(),
                                                      name, self.tvars);
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

    cdef PyFFEvaluatorObject *evaluator
    cdef evaluator_object
    
    def __init__(self, universe, options):
        TrajectoryGenerator.__init__(universe.options)
        evaluator_object = None
        self.evaluator = NULL

    cdef void initializeTrajectoryActions(self, char *name) except *:
        TrajectoryGenerator.initializeTrajectoryActions(self, name)
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
        self.evaluator.eval_func(self.evaluator, energy, conf_array,
                                 small_change)
        if self.lock_state != 1:
            PyUniverseSpec_StateLock(self.universe_spec, 2)
        if self.lock_state == -1:
            PyUniverseSpec_StateLock(self.universe_spec, -1)
        PyEval_RestoreThread(self.evaluator.tstate_save)

#
# Velocity Verlet integrator
#
cdef class VelocityVerletIntegrator(EnergyBasedTrajectoryGenerator):

    """
    Velocity-Verlet molecular dynamics integrator

    The integrator can handle fixed atoms, distance constraints,
    a thermostat, and a barostat, as well as any combination.
    It is fully thread-safe.

    The integration is started by calling the integrator object.
    All the keyword options (see documnentation of __init__) can be
    specified either when creating the integrator or when calling it.

    The following data categories and variables are available for
    output:

     - category "time": time

     - category "configuration": configuration and box size (for
       periodic universes)

     - category "velocities": atomic velocities

     - category "gradients": energy gradients for each atom

     - category "energy": potential and kinetic energy, plus
       extended-system energy terms if a thermostat and/or barostat
       are used

     - category "thermodynamic": temperature, volume (if a barostat
       is used) and pressure

     - category "auxiliary": extended-system coordinates if a thermostat
       and/or barostat are used
    """

    def __init__(self, universe, **options):
        """
        @param universe: the universe on which the integrator acts
        @type universe: L{MMTK.Universe}
        @keyword steps: the number of integration steps (default is 100)
        @type steps: C{int}
        @keyword delta_t: the time step (default is 1 fs)
        @type delta_t: C{float}
        @keyword actions: a list of actions to be executed periodically
                          (default is none)
        @type actions: C{list}
        @keyword threads: the number of threads to use in energy evaluation
                          (default set by MMTK_ENERGY_THREADS)
        @type threads: C{int}
        @keyword background: if True, the integration is executed as a
                             separate thread (default: False)
        @type background: C{bool}
        @keyword mpi_communicator: an MPI communicator object, or C{None},
                                   meaning no parallelization (default: C{None})
        @type mpi_communicator: C{Scientific.MPI.MPICommunicator}
        """
        TrajectoryGenerator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'mpi_communicator': None, 'actions': []}

    available_data = ['configuration', 'velocities', 'gradients',
                      'energy', 'thermodynamic', 'time', 'auxiliary']

    restart_data = ['configuration', 'velocities', 'energy',
                    'thermodynamic', 'auxiliary']

    # Cython compiler directives set for efficiency:
    # - No bound checks on index operations
    # - No support for negative indices
    # - Division uses C semantics
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef start(self):
        cdef N.ndarray[double, ndim=2] x, v, g, dv
        cdef N.ndarray[double, ndim=1] m
        cdef energy_data energy
        cdef double time, delta_t, ke
        cdef int natoms, nsteps, step
        cdef Py_ssize_t i, j, k
    
        # Check if velocities have been initialized
        if self.universe.velocities() is None:
            raise ValueError("no velocities")

        # Gather state variables and parameters
        configuration = self.universe.configuration()
        velocities = self.universe.velocities()
        gradients = ParticleProperties.ParticleVector(self.universe)
        masses = self.universe.masses()
        delta_t = self.getOption('delta_t')
        nsteps = self.getOption('steps')
        natoms = self.universe.numberOfAtoms()
    
        # For efficiency, the Cython code works at the array
        # level rather than at the ParticleProperty level.
        x = configuration.array
        v = velocities.array
        g = gradients.array
        m = masses.array
        dv = N.zeros((natoms, 3), N.float)

        # Ask for energy gradients to be calculated and stored in
        # the array g. Force constants are not requested.
        energy.gradients = <void *>g
        energy.gradient_fn = NULL
        energy.force_constants = NULL
        energy.fc_fn = NULL

        # Declare the variables accessible to trajectory actions.
        self.declareTrajectoryVariable_double(
            &time, "time", "Time: %lf\n", time_unit_name, PyTrajectory_Time)
        self.declareTrajectoryVariable_array(
            v, "velocities", "Velocities:\n", velocity_unit_name,
            PyTrajectory_Velocities)
        self.declareTrajectoryVariable_array(
            g, "gradients", "Energy gradients:\n", energy_gradient_unit_name,
            PyTrajectory_Gradients)
        self.declareTrajectoryVariable_double(
            &energy.energy,"potential_energy", "Potential energy: %lf\n",
            energy_unit_name, PyTrajectory_Energy)
        self.declareTrajectoryVariable_double(
            &ke, "kinetic_energy", "Kinetic energy: %lf\n",
            energy_unit_name, PyTrajectory_Energy)
        self.initializeTrajectoryActions("Velocity Verlet")

        # Acquire the write lock of the universe. This is necessary to
        # make sure that the integrator's modifications to positions
        # and velocities are synchronized with other threads that
        # attempt to use or modify these same values.
        #
        # Note that the write lock will be released temporarily
        # for trajectory actions. It will also be converted to
        # a read lock temporarily for energy evaluation. This
        # is taken care of automatically by the respective methods
        # of class EnergyBasedTrajectoryGenerator.
        self.acquireWriteLock()

        # Preparation: Calculate initial half-step accelerations
        # and run the trajectory actions on the initial state.
        ke = 0.
        self.foldCoordinatesIntoBox()
        self.calculateEnergies(x, &energy, 0)
        for i in range(natoms):
            for j in range(3):
                dv[i, j] = -0.5*delta_t*g[i, j]/m[i]
                ke += 0.5*m[i]*v[i, j]*v[i, j]
        self.trajectoryActions(step)

        # Main integration loop
        time = 0.
        for step in range(nsteps):
            # First half-step
            for i in range(natoms):
                for j in range(3):
                    v[i, j] += dv[i, j]
                    x[i, j] += delta_t*v[i, j]
            # Mid-step energy calculation
            self.foldCoordinatesIntoBox()
            self.calculateEnergies(x, &energy, 1)
            # Second half-step
            ke = 0.
            for i in range(natoms):
                for j in range(3):
                    dv[i, j] = -0.5*delta_t*g[i, j]/m[i]
                    v[i, j] += dv[i, j]
                    ke += 0.5*m[i]*v[i, j]*v[i, j]
            time += delta_t
            self.trajectoryActions(step)

        # Release the write lock.
        self.releaseWriteLock()

        # Finalize all trajectory actions (close files etc.)
        self.finalizeTrajectoryActions(nsteps)
