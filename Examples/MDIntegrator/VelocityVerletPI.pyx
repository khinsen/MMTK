# This module implements a Velocity Verlet integrator for path integral systems.
#
# Written by Konrad Hinsen
#

#cython: boundscheck=False, wraparound=False, cdivision=True

import numpy as N
cimport numpy as N

cimport MMTK_trajectory_generator
from MMTK import Units, ParticleProperties, Features
import MMTK_trajectory
import MMTK_forcefield

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"


#
# Velocity Verlet integrator
#
cdef class VelocityVerletIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

    """
    Velocity-Verlet molecular dynamics integrator
    with path integral support

    The integrator is fully thread-safe.

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
        """
        MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator.__init__(
            self, universe, options, "Velocity Verlet integrator")
        # Supported features: PathIntegrals
        self.features = [Features.PathIntegralsWithSpringTermsFeature]

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'actions': []}

    available_data = ['configuration', 'velocities', 'gradients',
                      'energy', 'time']

    restart_data = ['configuration', 'velocities', 'energy']

    cdef start(self):
        cdef N.ndarray[double, ndim=2] x, v, g, dv
        cdef N.ndarray[double, ndim=1] m
        cdef energy_data energy
        cdef double time, delta_t, ke
        cdef int nbeads, nsteps, step
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
        nbeads = self.universe.numberOfPoints()
    
        # For efficiency, the Cython code works at the array
        # level rather than at the ParticleProperty level.
        x = configuration.array
        v = velocities.array
        g = gradients.array
        m = masses.array
        dv = N.zeros((nbeads, 3), N.float)

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
        self.initializeTrajectoryActions()

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
        for i in range(nbeads):
            for j in range(3):
                dv[i, j] = -0.5*delta_t*g[i, j]/m[i]
                ke += 0.5*m[i]*v[i, j]*v[i, j]
        self.trajectoryActions(step)

        # Main integration loop
        time = 0.
        for step in range(nsteps):
            # First half-step
            for i in range(nbeads):
                for j in range(3):
                    v[i, j] += dv[i, j]
                    x[i, j] += delta_t*v[i, j]
            # Mid-step energy calculation
            self.foldCoordinatesIntoBox()
            self.calculateEnergies(x, &energy, 1)
            # Second half-step
            ke = 0.
            for i in range(nbeads):
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
