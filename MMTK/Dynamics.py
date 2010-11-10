# This module implements MD integrators
#
# Written by Konrad Hinsen
#

"""
Molecular Dynamics integrators
"""

__docformat__ = 'epytext'

from MMTK import Environment, Features, Trajectory, Units
import MMTK_dynamics
from Scientific import N

#
# Integrator base class
#
class Integrator(Trajectory.TrajectoryGenerator):

    def __init__(self, universe, options):
        Trajectory.TrajectoryGenerator.__init__(self, universe, options)

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'mpi_communicator': None, 'actions': []}

    available_data = ['configuration', 'velocities', 'gradients',
                      'energy', 'thermodynamic', 'time', 'auxiliary']

    restart_data = ['configuration', 'velocities', 'energy',
                    'thermodynamic', 'auxiliary']

    def __call__(self, options):
        raise AttributeError

#
# Velocity-Verlet integrator
#
class VelocityVerletIntegrator(Integrator):

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
        Integrator.__init__(self, universe, options)
        self.features = [Features.FixedParticleFeature,
                         Features.NoseThermostatFeature,
                         Features.AndersenBarostatFeature,
                         Features.PathIntegralsWithSpringTermsFeature,
                         Features.DistanceConstraintsFeature]

    def __call__(self, **options):
        """
        Run the integrator. The keyword options are the same as described
        under L{__init__}.
        """
        self.setCallOptions(options)
        used_features = Features.checkFeatures(self, self.universe)
        configuration = self.universe.configuration()
        velocities = self.universe.velocities()
        if velocities is None:
            raise ValueError("no velocities")
        masses = self.universe.masses()
        fixed = self.universe.getAtomBooleanArray('fixed')
        nt = self.getOption('threads')
        comm = self.getOption('mpi_communicator')
        evaluator = self.universe.energyEvaluator(threads=nt,
                                                  mpi_communicator=comm)
        evaluator = evaluator.CEvaluator()
        constraints, const_distances_sq, c_blocks = \
                     _constraintArrays(self.universe)
                                                  
        type = 'NVE'
        if Features.NoseThermostatFeature in used_features:
            type = 'NVT'
            thermostat = self.universe.environmentObjectList(
                                          Environment.NoseThermostat)[0]
            t_parameters = thermostat.parameters
            t_coordinates = thermostat.coordinates
        else:
            t_parameters = N.zeros((0,), N.Float)
            t_coordinates = N.zeros((2,), N.Float)
        if Features.AndersenBarostatFeature in used_features:
            if self.universe.cellVolume() is None:
                raise ValueError("Barostat requires finite volume universe")
            if type == 'NVE':
                type = 'NPH'
            else:
                type = 'NPT'
            barostat = self.universe.environmentObjectList(
                                          Environment.AndersenBarostat)[0]
            b_parameters = barostat.parameters
            b_coordinates = barostat.coordinates
        else:
            b_parameters = N.zeros((0,), N.Float)
            b_coordinates = N.zeros((1,), N.Float)

        args = (self.universe,
                configuration.array, velocities.array,
                masses.array, fixed.array, evaluator,
                constraints, const_distances_sq, c_blocks,
                t_parameters, t_coordinates,
                b_parameters, b_coordinates,
                self.getOption('delta_t'), self.getOption('first_step'),
                self.getOption('steps'), self.getActions(),
                type + ' dynamics trajectory with ' +
                self.optionString(['delta_t', 'steps']))
        return self.run(MMTK_dynamics.integrateVV, args)

#
# Velocity scaling, removal of global translation/rotation
#
class VelocityScaler(Trajectory.TrajectoryAction):

    """
    Periodic velocity scaling action

    A VelocityScaler object is used in the action list of
    a VelocityVerletIntegrator. It rescales all atomic velocities
    by a common factor to make the temperature of the system equal
    to a predefined value.
    """

    def __init__(self, temperature, window=0., first=0, last=None, skip=1):
        """
        @param temperature: the temperature value to which the velocities
                            should be scaled
        @type temperature: C{float}
        @param window: the deviation from the ideal temperature that
                       is tolerated in either direction before rescaling
                       takes place
        @type window: C{float}
        @param first: the number of the first step at which the action
                      is executed
        @type first: C{int}
        @param last: the number of the first step at which the action
                     is no longer executed. A value of C{None} indicates
                     that the action should be executed indefinitely.
        @type last: C{int}
        @param skip: the number of steps to skip between two applications
                     of the action
        @type skip: C{int}
        """
        Trajectory.TrajectoryAction.__init__(self, first, last, skip)
        self.parameters = N.array([temperature, window], N.Float)
        self.Cfunction = MMTK_dynamics.scaleVelocities

class Heater(Trajectory.TrajectoryAction):

    """
    Periodic heating action

    A Heater object us used in the action list of a VelocityVerletIntegrator.
    It scales the velocities to a temperature that increases over time.
    """

    def __init__(self, temp1, temp2, gradient, first=0, last=None, skip=1):
        """
        @param temp1: the temperature value to which the velocities
                      should be scaled initially
        @type temp1: C{float}
        @param temp2: the final temperature value to which the velocities
                      should be scaled
        @type temp2: C{float}
        @param gradient: the temperature gradient (in K/ps)
        @type gradient: C{float}
        @param first: the number of the first step at which the action
                      is executed
        @type first: C{int}
        @param last: the number of the first step at which the action
                     is no longer executed. A value of C{None} indicates
                     that the action should be executed indefinitely.
        @type last: C{int}
        @param skip: the number of steps to skip between two applications
                     of the action
        @type skip: C{int}
        """
        Trajectory.TrajectoryAction.__init__(self, first, last, skip)
        self.parameters = N.array([temp1, temp2, gradient], N.Float)
        self.Cfunction = MMTK_dynamics.heat

class BarostatReset(Trajectory.TrajectoryAction):

    """
    Barostat reset action

    A BarostatReset object is used in the action list of a
    VelocityVerletIntegrator. It resets the barostat coordinate
    to zero.
    """

    def __init__(self, first=0, last=None, skip=1):
        """
        @param first: the number of the first step at which the action
                      is executed
        @type first: C{int}
        @param last: the number of the first step at which the action
                     is no longer executed. A value of C{None} indicates
                     that the action should be executed indefinitely.
        @type last: C{int}
        @param skip: the number of steps to skip between two applications
                     of the action
        @type skip: C{int}
        """
        Trajectory.TrajectoryAction.__init__(self, first, last, skip)
        self.parameters = N.zeros((0,), N.Float)
        self.Cfunction = MMTK_dynamics.resetBarostat

class TranslationRemover(Trajectory.TrajectoryAction):

    """
    Action that eliminates global translation

    A TranslationRemover object is used in the action list of a
    VelocityVerletIntegrator. It subtracts the total velocity
    of the system from each atomic velocity.
    """

    def __init__(self, first=0, last=None, skip=1):
        """
        @param first: the number of the first step at which the action
                      is executed
        @type first: C{int}
        @param last: the number of the first step at which the action
                     is no longer executed. A value of C{None} indicates
                     that the action should be executed indefinitely.
        @type last: C{int}
        @param skip: the number of steps to skip between two applications
                     of the action
        @type skip: C{int}
        """
        Trajectory.TrajectoryAction.__init__(self, first, last, skip)
        self.parameters = None
        self.Cfunction = MMTK_dynamics.removeTranslation

class RotationRemover(Trajectory.TrajectoryAction):

    """
    Action that eliminates global rotation

    A RotationRemover object is used in the action list of a
    VelocityVerletIntegrator. It adjusts the atomic velocities
    such that the total angular momentum is zero.
    """

    def __init__(self, first=0, last=None, skip=1):
        """
        @param first: the number of the first step at which the action
                      is executed
        @type first: C{int}
        @param last: the number of the first step at which the action
                     is no longer executed. A value of C{None} indicates
                     that the action should be executed indefinitely.
        @type last: C{int}
        @param skip: the number of steps to skip between two applications
                     of the action
        @type skip: C{int}
        """
        Trajectory.TrajectoryAction.__init__(self, first, last, skip)
        self.parameters = None
        self.Cfunction = MMTK_dynamics.removeRotation

#
# Construct constraint arrays
#
def _constraintArrays(universe):
    nc = universe.numberOfDistanceConstraints()
    constraints = N.zeros((nc, 2), N.Int)
    const_distances_sq = N.zeros((nc, ), N.Float)
    if nc > 0:
        i = 0
        c_blocks = [0]
        for o in universe.objectList():
            for c in o.distanceConstraintList():
                constraints[i, 0] = c[0].index
                constraints[i, 1] = c[1].index
                const_distances_sq[i] = c[2]*c[2]
                i = i + 1
            c_blocks.append(i)
        c_blocks = N.array(c_blocks)
    else:
        c_blocks = N.zeros((1,), N.Int)
    return constraints, const_distances_sq, c_blocks

#
# Enforce distance constraints for current configuration
#
def enforceConstraints(universe, configuration=None):
    constraints, const_distances_sq, c_blocks = _constraintArrays(universe)
    if len(constraints) == 0:
        return
    if configuration is None:
        configuration = universe.configuration()
    MMTK_dynamics.enforceConstraints(universe._spec, configuration.array,
                                     universe.masses().array,
                                     constraints, const_distances_sq, c_blocks)

#
# Project velocity vector onto the constraint surface
#
def projectVelocities(universe, velocities):
    constraints, const_distances_sq, c_blocks = _constraintArrays(universe)
    if len(constraints) == 0:
        return
    MMTK_dynamics.projectVelocities(universe._spec,
                                    universe.configuration().array,
                                    velocities.array, universe.masses().array,
                                    constraints, const_distances_sq, c_blocks)
