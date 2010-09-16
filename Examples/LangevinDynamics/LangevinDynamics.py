# This module implements a Langevin integrator.
#
# Written by Konrad Hinsen
#

from MMTK import Dynamics, Environment, Features, Trajectory, \
                 Units, ParticleProperties
import MMTK_langevin
import MMTK_forcefield

#
# Langevin integrator
#
class LangevinIntegrator(Dynamics.Integrator):

    def __init__(self, universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)
        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
        # Get the universe variables needed by the integrator
        configuration = self.universe.configuration()
        masses = self.universe.masses()
        velocities = self.universe.velocities()
        if velocities is None:
            raise ValueError("no velocities")

        # Get the friction coefficients. First check if a keyword argument
        # 'friction' was given to the integrator. Its value can be a
        # ParticleScalar or a plain number (used for all atoms). If no
        # such argument is given, collect the values of the attribute
        # 'friction' from all atoms (default is zero).
        try:
            friction = self.getOption('friction')
        except KeyError:
            friction = self.universe.getParticleScalar('friction')
        if not ParticleProperties.isParticleProperty(friction):
            var = ParticleProperties.ParticleScalar(self.universe)
            var.array[:] = friction
            friction = var

        # Construct a C evaluator object for the force field, using
        # the specified number of threads or the default value
        nt = self.getOption('threads')
        evaluator = self.universe.energyEvaluator(threads=nt).CEvaluator()

        # Run the C integrator
        MMTK_langevin.integrateLD(self.universe,
                                  configuration.array, velocities.array,
                                  masses.array, friction.array, evaluator,
                                  self.getOption('temperature'),
                                  self.getOption('delta_t'),
                                  self.getOption('steps'),
                                  self.getActions())
