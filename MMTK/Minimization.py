# This module implements energy minimizers.
#
# Written by Konrad Hinsen
#

"""
Energy minimizers
"""

__docformat__ = 'epytext'

from MMTK import Features, Trajectory, Units
from MMTK_minimization import conjugateGradient, steepestDescent

#
# Minimizer base class
#
class Minimizer(Trajectory.TrajectoryGenerator):

    def __init__(self, universe, options):
	Trajectory.TrajectoryGenerator.__init__(self, universe, options)

    default_options = {'steps': 100, 'step_size': 0.02*Units.Ang,
		       'convergence': 0.01*Units.kJ/(Units.mol*Units.nm),
                       'background': False, 'threads': None,
                       'mpi_communicator': None, 'actions': []}

    available_data = ['energy', 'configuration', 'gradients']

    restart_data = ['configuration', 'energy']

    def __call__(self, options):
	raise AttributeError

#
# Steepest descent minimizer
#
class SteepestDescentMinimizer(Minimizer):

    """
    Steepest-descent minimizer

    The minimizer can handle fixed atoms, but no distance constraints.
    It is fully thread-safe.

    The minimization is started by calling the minimizer object.
    All the keyword options (see documnentation of __init__) can be
    specified either when creating the minimizer or when calling it.

    The following data categories and variables are available for
    output:

     - category "configuration": configuration and box size (for
       periodic universes)

     - category "gradients": energy gradients for each atom

     - category "energy": potential energy and
                          norm of the potential energy gradient
    """

    def __init__(self, universe, **options):
        """
        @param universe: the universe on which the integrator acts
        @type universe: L{MMTK.Universe}
        @keyword steps: the number of minimization steps (default is 100)
        @type steps: C{int}
        @keyword step_size: the initial size of a minimization step
                            (default is 0.002 nm)
        @type step_size: C{float}
        @keyword convergence: the root-mean-square gradient length at which
                              minimization stops (default is 0.01 kJ/mol/nm)
        @type convergence: C{float}
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
	Minimizer.__init__(self, universe, options)
	self.features = [Features.FixedParticleFeature,
			 Features.NoseThermostatFeature,
                         Features.PathIntegralsFeature,
                         Features.AndersenBarostatFeature]

    def __call__(self, **options):
        """
        Run the minimizer. The keyword options are the same as described
        under L{__init__}.
        """
	self.setCallOptions(options)
	Features.checkFeatures(self, self.universe)
	configuration = self.universe.configuration()
	fixed = self.universe.getAtomBooleanArray('fixed')
        nt = self.getOption('threads')
        comm = self.getOption('mpi_communicator')
	evaluator = self.universe.energyEvaluator(threads=nt,
                                                  mpi_communicator=comm)
        evaluator = evaluator.CEvaluator()
	args = (self.universe,
                configuration.array, fixed.array, evaluator,
                self.getOption('steps'), self.getOption('step_size'),
                self.getOption('convergence'), self.getActions(),
                'Steepest descent minimization with ' +
                self.optionString(['convergence', 'step_size', 'steps']))
        return self.run(steepestDescent, args)

#
# Conjugate gradient minimizer
#
class ConjugateGradientMinimizer(Minimizer):

    """
    Conjugate gradient minimizer

    The minimizer can handle fixed atoms, but no distance constraints.
    It is fully thread-safe.

    The minimization is started by calling the minimizer object.
    All the keyword options can be specified either when
    creating the minimizer or when calling it.

    The following data categories and variables are available for
    output:

     - category "configuration": configuration and box size (for
       periodic universes)

     - category "gradients": energy gradients for each atom

     - category "energy": potential energy and
                          norm of the potential energy gradient
    """

    def __init__(self, universe, **options):
        """
        @param universe: the universe on which the integrator acts
        @type universe: L{MMTK.Universe}
        @keyword steps: the number of minimization steps (default is 100)
        @type steps: C{int}
        @keyword step_size: the initial size of a minimization step
                            (default is 0.002 nm)
        @type step_size: C{float}
        @keyword convergence: the root-mean-square gradient length at which
                              minimization stops (default is 0.01 kJ/mol/nm)
        @type convergence: C{float}
        @keyword actions: a list of actions to be executed periodically
                          (default is none)
        @type actions: C{list}
        @keyword threads: the number of threads to use in energy evaluation
                          (default set by MMTK_ENERGY_THREADS)
        @type threads: C{int}
        @keyword background: if True, the integration is executed as a
                             separate thread (default: False)
        @type background: C{bool}
        """
	Minimizer.__init__(self, universe, options)
	self.features = [Features.FixedParticleFeature,
                         Features.PathIntegralsFeature,
			 Features.NoseThermostatFeature]

    def __call__(self, **options):
        """
        Run the minimizer. The keyword options are the same as described
        under L{__init__}.
        """
	self.setCallOptions(options)
	Features.checkFeatures(self, self.universe)
	configuration = self.universe.configuration()
	fixed = self.universe.getAtomBooleanArray('fixed')
        nt = self.getOption('threads')
	evaluator = self.universe.energyEvaluator(threads=nt).CEvaluator()
	args =(self.universe,
               configuration.array, fixed.array, evaluator,
               self.getOption('steps'), self.getOption('step_size'),
               self.getOption('convergence'), self.getActions(),
               'Conjugate gradient minimization with ' +
               self.optionString(['convergence', 'step_size', 'steps']))
        return self.run(conjugateGradient, args)
