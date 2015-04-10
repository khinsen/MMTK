# This module implements a Velocity Verlet integrator for path integral systems.
#
# Written by Konrad Hinsen
#

"""
Velocity Verlet integrator for path integral systems
"""

__docformat__ = 'restructuredtext'

import numpy as N
cimport numpy as N
import numpy.linalg as LA
import cython

from MMTK import Units, ParticleProperties, Features, Environment
import MMTK_trajectory
import MMTK_forcefield
import MMTK.PIIntegratorSupport
cimport MMTK.PIIntegratorSupport

import MMTK.mtrand
cimport MMTK.mtrand

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

cdef extern from "stdlib.h":
    cdef double fabs(double)
    cdef double sqrt(double)
    cdef double sin(double)
    cdef double cos(double)
    cdef double exp(double)
    cdef double M_PI

cdef double hbar = Units.hbar
cdef double k_B = Units.k_B

#
# Velocity Verlet integrator
#
cdef class PICartesianIntegrator(MMTK.PIIntegratorSupport.PIIntegrator):

    """
    Velocity-Verlet molecular dynamics integrator with path integral support

    The integrator is fully thread-safe.

    The integration is started by calling the integrator object.
    All the keyword options (see documentation of __init__) can be
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

     - category "thermodynamic": temperature

     - category "auxiliary": primitive and virial quantum energy estimators

    """

    def __init__(self, universe, **options):
        """
        :param universe: the universe on which the integrator acts
        :type universe: MMTK.Universe
        :keyword steps: the number of integration steps (default is 100)
        :type steps: int
        :keyword delta_t: the time step (default is 1 fs)
        :type delta_t: float
        :keyword actions: a list of actions to be executed periodically
                          (default is none)
        :type actions: list
        :keyword threads: the number of threads to use in energy evaluation
                          (default set by MMTK_ENERGY_THREADS)
        :type threads: int
        :keyword background: if True, the integration is executed as a
                             separate thread (default: False)
        :type background: bool
        """
        MMTK.PIIntegratorSupport.PIIntegrator.__init__(
            self, universe, options, "Velocity Verlet integrator")
        # Supported features: PathIntegrals
        self.features = [Features.PathIntegralsWithSpringTermsFeature]

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'frozen_subspace': None, 'actions': []}

    available_data = ['time', 'configuration', 'velocities', 'gradients',
                      'energy', 'thermodynamic', 'auxiliary']

    restart_data = ['configuration', 'velocities', 'energy']

    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v,
                              N.ndarray[double, ndim=1] m, N.ndarray[short, ndim=2] bd,
                              double dt, double beta):
        pass

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
        cdef N.ndarray[short, ndim=2] bd
        cdef N.ndarray[double, ndim=3] ss
        cdef energy_data energy
        cdef double time, delta_t, ke, beta, temperature
        cdef double qe_prim, qe_vir, qe_cvir
        cdef int natoms, nbeads, nsteps, step, df, cdf
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
        nbeads = self.universe.numberOfPoints()
        bd = self.evaluator_object.global_data.get('bead_data')
        pi_environment = self.universe.environmentObjectList(Environment.PathIntegrals)[0]
        beta = pi_environment.beta

        # Check if there is a frozen_subspace
        subspace = self.getOption('frozen_subspace')
        if subspace is None:
            ss = N.zeros((0, nbeads, 3), N.float)
            df = 3*nbeads
            cdf = 3*natoms
        else:
            ss = subspace.getBasis().array
            df = 3*nbeads-ss.shape[0]
            cdf = self.centroidDegreesOfFreedom(subspace, bd)

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
        self.declareTrajectoryVariable_double(
            &temperature, "temperature", "Temperature: %lf\n",
            temperature_unit_name, PyTrajectory_Thermodynamic)
        self.declareTrajectoryVariable_double(
            &qe_prim, "quantum_energy_primitive",
            "Primitive quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.declareTrajectoryVariable_double(
            &qe_vir, "quantum_energy_virial",
            "Virial quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
        self.declareTrajectoryVariable_double(
            &qe_cvir, "quantum_energy_centroid_virial",
            "Centroid virial quantum energy estimator: %lf\n",
            energy_unit_name, PyTrajectory_Auxiliary)
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
        self.foldCoordinatesIntoBox()
        self.calculateEnergies(x, &energy, 0)
        self.freeze(v, ss)

        se = self.springEnergyCartesian(x, m, bd, beta)
        qe_prim = energy.energy - se + 0.5*df/beta
        qe_vir = energy.energy - 0.5*energy.virial
        qe_cvir = energy.energy \
                  - 0.5*self.centroidVirial(x, g, bd) \
                  + 0.5*cdf/beta

        ke = 0.
        for i in range(nbeads):
            for j in range(3):
                dv[i, j] = -0.5*delta_t*g[i, j]/m[i]
                ke += 0.5*m[i]*v[i, j]*v[i, j]
        temperature = 2.*ke/(df*k_B)

        # Main integration loop
        time = 0.
        self.trajectoryActions(0)
        for step in range(nsteps):
            # First application of thermostat
            self.applyThermostat(v, m, bd, delta_t, beta)
            # First integration half-step
            for i in range(nbeads):
                for j in range(3):
                    v[i, j] += dv[i, j]
            # Remove frozen subspace
            self.freeze(v, ss)
            # Update of positions
            for i in range(nbeads):
                for j in range(3):
                    x[i, j] += delta_t*v[i, j]
            # Mid-step energy calculation
            self.foldCoordinatesIntoBox()
            self.calculateEnergies(x, &energy, 1)
            # Quantum energy estimators
            se = self.springEnergyCartesian(x, m, bd, beta)
            qe_prim = energy.energy - se + 0.5*df/beta
            qe_vir = energy.energy - 0.5*energy.virial
            qe_cvir = energy.energy \
                      - 0.5*self.centroidVirial(x, g, bd) \
                      + 0.5*cdf/beta
            # Second integration half-step
            for i in range(nbeads):
                for j in range(3):
                    dv[i, j] = -0.5*delta_t*g[i, j]/m[i]
                    v[i, j] += dv[i, j]
            # Second application of thermostat
            self.applyThermostat(v, m, bd, delta_t, beta)
            # Remove frozen subspace
            self.freeze(v, ss)
            # Calculate kinetic energy
            ke = 0.
            for i in range(nbeads):
                for j in range(3):
                    ke += 0.5*m[i]*v[i, j]*v[i, j]
            temperature = 2.*ke/(df*k_B)
            time += delta_t
            self.trajectoryActions(step+1)

        # Release the write lock.
        self.releaseWriteLock()

        # Finalize all trajectory actions (close files etc.)
        self.finalizeTrajectoryActions(nsteps)


#
# Velocity Verlet integrator with a colored-noise Langevin thermostat
#
cdef class PILangevinCartesianIntegrator(PICartesianIntegrator):

    """
    Velocity-Verlet molecular dynamics integrator with path integral support
    and a Langevin thermostat.

    This integrator works like PICartesianIntegrator, but has
    an additional open "friction_matrix", which is a square array.

    """

    cdef N.ndarray C1, C2
    cdef s, xi, temp

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v,
                              N.ndarray[double, ndim=1] m, N.ndarray[short, ndim=2] bd,
                              double dt, double beta):
        cdef N.ndarray[double, ndim=2] C1 = self.C1
        cdef N.ndarray[double, ndim=2] C2 = self.C2
        cdef N.ndarray[double, ndim=3] s = self.s
        cdef N.ndarray[double, ndim=1] xi = self.xi
        cdef N.ndarray[double, ndim=1] temp = self.temp
        cdef double mb
        cdef int nbeads = v.shape[0]
        cdef int ns = self.xi.shape[0]
        cdef int nb
        cdef int i, j, k, l
        for i in range(nbeads):
            nb = bd[i, 1]
            mb = sqrt(1./(beta*m[i]))
            for j in range(3):
                for k in range(ns):
                    temp[k] = C1[k, 0]*v[i, j] + mb*C2[k, 0]*MMTK.mtrand.standard_normal()
                    for l in range(1, ns):
                        temp[k] += C1[k, l]*s[l-1, i, j] + mb*C2[k, l]*MMTK.mtrand.standard_normal()
                v[i, j] = temp[0]
                for k in range(1, ns):
                    s[k-1, i, j] = temp[k]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef start(self):
        friction = self.getOption('friction_matrix')
        assert isinstance(friction, N.ndarray)
        assert len(friction.shape) == 2
        assert friction.shape[0] == friction.shape[1]
        nbeads = self.universe.numberOfPoints()
        delta_t = self.getOption('delta_t')
        e, v = LA.eig(-0.5*delta_t*N.transpose(friction))
        self.C1 = N.dot(v*N.exp(e), LA.inv(v))
        C2_sq = N.identity(len(self.C1))-N.dot(N.transpose(self.C1), self.C1)
        self.C2 = N.transpose(LA.cholesky(C2_sq))
        self.s = N.zeros((len(self.C1)-1, nbeads, 3), N.float)
        self.xi = N.zeros((len(self.C1),), N.float)
        self.temp = N.zeros((len(self.C1),), N.float)
        PICartesianIntegrator.start(self)

