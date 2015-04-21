# This module implements path integral MD integrator using normal mode coordinates
#
# Written by Konrad Hinsen
#

#cython: boundscheck=False, wraparound=False, cdivision=True

"""
Path integral MD integrator using normal-mode coordinates
"""

__docformat__ = 'restructuredtext'

import numpy as N
cimport numpy as N

from MMTK import Units, ParticleProperties, Features, Environment
import MMTK_trajectory
import MMTK_forcefield
import MMTK_universe
import MMTK.PIIntegratorSupport
cimport MMTK.PIIntegratorSupport
import numbers

import MMTK.mtrand
cimport MMTK.mtrand

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

cdef extern from "fftw3.h":
    ctypedef struct fftw_complex
    ctypedef void *fftw_plan
    cdef int FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
    cdef void fftw_execute(fftw_plan p)
    cdef fftw_plan fftw_plan_dft_1d(int n, fftw_complex *data_in, fftw_complex *data_out,
                                    int sign, int flags)
    cdef void fftw_destroy_plan(fftw_plan p)

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
# Velocity Verlet integrator in normal-mode coordinates
#
cdef class PINormalModeIntegrator(MMTK.PIIntegratorSupport.PIIntegrator):

    """
    Molecular dynamics integrator for path integral systems using
    normal-mode coordinates.

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

    cdef N.ndarray workspace1, workspace2
    cdef double *workspace_ptr_1
    cdef double *workspace_ptr_2

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
            self, universe, options, "Path integral normal-mode integrator")
        # Supported features: PathIntegrals
        self.features = [Features.PathIntegralsFeature]

    default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                       'background': False, 'threads': None,
                       'frozen_subspace': None, 'actions': []}

    available_data = ['time', 'configuration', 'velocities', 'gradients',
                      'energy', 'thermodynamic', 'auxiliary']

    restart_data = ['configuration', 'velocities', 'energy']

    # The implementation of the equations of motion follows the article
    #   Ceriotti et al., J. Chem. Phys. 133, 124104 (2010)
    # with the following differences:
    # 1) All the normal mode coordinates are larger by a factor sqrt(nbeads),
    #    and the non-real ones (k != 0, k != n/2) are additionally smaller by
    #    sqrt(2).
    # 2) The spring energy is smaller by a factor of nbeads to take
    #    into account the factor nbeads in Eq. (3) of the paper cited above.
    #    The potential energy of the system is also smaller by a factor of
    #    nbeads compared to the notation in this paper.
    # 3) Velocities are used instead of momenta in the integrator.
    # 4) Eq. (18) is also used for odd n, ignoring the k = n/2 case.

    cdef cartesianToNormalMode(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] nmc,
                               int bead_index, int nb):
        cdef double *w1 = self.workspace_ptr_1
        cdef double *w2 = self.workspace_ptr_2
        cdef fftw_plan p
        cdef int i, j
        if nb == 1:
            for i in range(3):
                nmc[i, bead_index] = x[bead_index, i]
        else:
            p = fftw_plan_dft_1d(nb, <fftw_complex *>w1, <fftw_complex *>w2,
                                 FFTW_FORWARD, FFTW_ESTIMATE)
            for i in range(3):
                for j in range(nb):
                    w1[2*j] = x[bead_index+j, i]
                    w1[2*j+1] = 0.
                fftw_execute(p)
                nmc[i, bead_index+0] = w2[0]
                for j in range(1, (nb+1)/2):
                    nmc[i, bead_index+j] = w2[2*j]
                    nmc[i, bead_index+nb-j] = w2[2*j+1]
                if nb % 2 == 0:
                    nmc[i, bead_index+nb/2] = w2[nb]
            fftw_destroy_plan(p)

    cdef normalModeToCartesian(self, N.ndarray[double, ndim=2] x, N.ndarray[double, ndim=2] nmc,
                               int bead_index, int nb):
        cdef double *w1 = self.workspace_ptr_1
        cdef double *w2 = self.workspace_ptr_2
        cdef fftw_plan p
        cdef int i, j
        if nb == 1:
            for i in range(3):
                x[bead_index, i] = nmc[i, bead_index]
        else:
            p = fftw_plan_dft_1d(nb, <fftw_complex *>w1, <fftw_complex *>w2,
                                 FFTW_BACKWARD, FFTW_ESTIMATE)
            for i in range(3):
                w1[0] = nmc[i, bead_index+0]
                w1[1] = 0.
                for j in range(1, (nb+1)/2):
                    w1[2*j] = nmc[i, bead_index+j]
                    w1[2*j+1] = nmc[i, bead_index+nb-j]
                    w1[2*nb-2*j] = w1[2*j]
                    w1[2*nb-2*j+1] = -w1[2*j+1]
                if nb % 2 == 0:
                    w1[nb] = nmc[i, bead_index+nb/2]
                    w1[nb+1] = 0.
                fftw_execute(p)
                for j in range(nb):
                    x[bead_index+j, i] = w2[2*j]/nb
            fftw_destroy_plan(p)

    cdef void propagateOscillators(self, N.ndarray[double, ndim=2] nmc,
                                   N.ndarray[double, ndim=2] nmv,
                                   int bead_index, int nb, double beta, double dt):
        cdef double omega_n = nb/(beta*hbar)
        cdef double omega_k, omega_k_dt, s, c
        cdef double temp
        cdef int i, k
        for i in range(3):
            nmc[i, bead_index] += dt*nmv[i, bead_index]
            for k in range(1, nb):
                omega_k = 2.*omega_n*sin(k*M_PI/nb)
                omega_k_dt = omega_k*dt
                s = sin(omega_k_dt)
                c = cos(omega_k_dt)
                temp = c*nmv[i, bead_index+k]-omega_k*s*nmc[i, bead_index+k]
                nmc[i, bead_index+k] = s*nmv[i, bead_index+k]/omega_k + c*nmc[i, bead_index+k]
                nmv[i, bead_index+k] = temp

    cdef double springEnergyNormalModes(self, N.ndarray[double, ndim=2] nmc,
                                        N.ndarray[double, ndim=1] m,
                                        N.ndarray[short, ndim=2] bd,
                                        double beta):
        cdef int i, j, k, nb
        cdef double sumsq
        cdef double omega_n, omega_k
        cdef double e = 0.
        for i in range(nmc.shape[1]):
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                omega_n = nb/(beta*hbar)
                # Start at j=1 because the contribution from the centroid is zero
                for j in range(1, nb):
                    omega_k = 2.*omega_n*sin(j*M_PI/nb)
                    sumsq = 0.
                    for k in range(3):
                        sumsq += nmc[k, i+j]*nmc[k, i+j]
                    # j=nb/2 corresponds to the real-valued coordinate at
                    # the maximal frequency.
                    if nb % 2 == 0 and j == nb/2:
                        sumsq *= 0.5
                    e += m[i]*sumsq*omega_k*omega_k/nb
        return e

    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v, N.ndarray[double, ndim=2] nmv,
                              N.ndarray[double, ndim=1] m, N.ndarray[short, ndim=2] bd,
                              double dt, double beta):
        pass

    cdef start(self):
        cdef N.ndarray[double, ndim=2] x, v, g, dv, nmc, nmv
        cdef N.ndarray[double, ndim=1] m
        cdef N.ndarray[short, ndim=2] bd
        cdef N.ndarray[double, ndim=3] ss
        cdef energy_data energy
        cdef double time, delta_t, ke, ke_nm, se, beta, temperature
        cdef double qe_prim, qe_vir, qe_cvir
        cdef int natoms, nbeads, nsteps, step, df, cdf, nb
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
        nmc = N.zeros((3, nbeads), N.float)
        nmv = N.zeros((3, nbeads), N.float)

        # Allocate workspace for Fourier transforms
        nb_max = bd[:, 1].max()
        self.workspace1 = N.zeros((2*nb_max,), N.float)
        self.workspace_ptr_1 = <double *>self.workspace1.data
        self.workspace2 = N.zeros((2*nb_max,), N.float)
        self.workspace_ptr_2 = <double *>self.workspace2.data

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
            &se, "spring_energy", "Spring energy: %lf\n",
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

        for i in range(nbeads):
            if bd[i, 0] == 0:
                self.fixBeadPositions(x, i, bd[i, 1])
                self.cartesianToNormalMode(x, nmc, i, bd[i, 1])

        se = self.springEnergyNormalModes(nmc, m, bd, beta)
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

        # Check FFT
        if False:
            x_test = N.zeros((nbeads, 3), N.float)
            v_test = N.zeros((nbeads, 3), N.float)
            for i in range(nbeads):
                if bd[i, 0] == 0:
                    self.cartesianToNormalMode(x, nmc, i, bd[i, 1])
                    self.normalModeToCartesian(x_test, nmc, i, bd[i, 1])
                    self.cartesianToNormalMode(v, nmv, i, bd[i, 1])
                    self.normalModeToCartesian(v_test, nmv, i, bd[i, 1])
            for i in range(nbeads):
                for j in range(3):
                    assert fabs(x[i, j]-x_test[i, j]) < 1.e-7
                    assert fabs(v[i, j]-v_test[i, j]) < 1.e-7

        # Main integration loop
        time = 0.
        self.trajectoryActions(0)
        for step in range(nsteps):
            # First application of thermostat
            self.applyThermostat(v, nmv, m, bd, delta_t, beta)
            # First integration half-step
            for i in range(nbeads):
                for j in range(3):
                    v[i, j] += dv[i, j]
            # Remove frozen subspace
            self.freeze(v, ss)
            # Conversion to normal mode coordinates
            for i in range(nbeads):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bd[i, 0] == 0:
                    self.fixBeadPositions(x, i, bd[i, 1])
                    self.cartesianToNormalMode(x, nmc, i, bd[i, 1])
                    self.cartesianToNormalMode(v, nmv, i, bd[i, 1])
            # Harmonic oscillator time propagation
            for i in range(nbeads):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bd[i, 0] == 0:
                    self.propagateOscillators(nmc, nmv, i, bd[i, 1], beta, delta_t)
            # Conversion back to Cartesian coordinates
            for i in range(nbeads):
                # bd[i, 0] == 0 means "first bead of an atom"
                if bd[i, 0] == 0:
                    self.normalModeToCartesian(x, nmc, i, bd[i, 1])
                    self.normalModeToCartesian(v, nmv, i, bd[i, 1])
            # Mid-step energy calculation
            self.calculateEnergies(x, &energy, 1)
            # Quantum energy estimators
            se = self.springEnergyNormalModes(nmc, m, bd, beta)
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
            self.applyThermostat(v, nmv, m, bd, delta_t, beta)
            for i in range(nbeads):
                if bd[i, 0] == 0:
                    self.cartesianToNormalMode(v, nmv, i, bd[i, 1])
            # Remove frozen subspace
            self.freeze(v, ss)
            # Calculate kinetic energy
            ke = 0.
            for i in range(nbeads):
                for j in range(3):
                    ke += 0.5*m[i]*v[i, j]*v[i, j]
            temperature = 2.*ke/(df*k_B)
            if False:
                ke_nm = 0.
                for i in range(nbeads):
                    if bd[i, 0] == 0:
                        nb = bd[i, 1]
                        for j in range(3):
                            for k in range(nb):
                                if k == 0 or (nb % 2 == 0 and k == nb/2):
                                    ke_nm += 0.5*m[i]*nmv[j, i+k]*nmv[j, i+k]/nb
                                else:
                                    ke_nm += m[i]*nmv[j, i+k]*nmv[j, i+k]/nb
                assert fabs(ke-ke_nm) < 1.e-7
            # End of time step
            time += delta_t
            self.foldCoordinatesIntoBox()
            self.trajectoryActions(step+1)

        # Release the write lock.
        self.releaseWriteLock()

        # Finalize all trajectory actions (close files etc.)
        self.finalizeTrajectoryActions(nsteps)

        # Deallocate the Fourier transform workspace
        self.workspace_ptr_1 = NULL
        self.workspace_ptr_2 = NULL
        self.workspace1 = None
        self.workspace2 = None


#
# Velocity Verlet integrator in normal-mode coordinates
# with a Langevin thermostat
#
cdef class PILangevinNormalModeIntegrator(PINormalModeIntegrator):

    """
    Molecular dynamics integrator for path integral systems using
    normal-mode coordinates and a Langevin thermostat.

    This integrator works like PINormalModeIntegrator, but has
    an additional option "centroid_friction", which is a ParticleScalar
    (one friction constant per atom) or a plain number.

    """

    cdef N.ndarray gamma

    cdef void applyThermostat(self, N.ndarray[double, ndim=2] v, N.ndarray[double, ndim=2] nmv,
                              N.ndarray[double, ndim=1] m, N.ndarray[short, ndim=2] bd,
                              double dt, double beta):
        cdef N.ndarray[double, ndim=1] g = self.gamma
        cdef int nbeads = v.shape[0]
        cdef double f, c1, c2
        cdef double omega_n, mb
        cdef int i, j, k, nb
        for i in range(nbeads):
            # bd[i, 0] == 0 means "first bead of an atom"
            if bd[i, 0] == 0:
                nb = bd[i, 1]
                # Conversion to normal mode coordinates
                self.cartesianToNormalMode(v, nmv, i, nb)
                # Modify velocities
                omega_n = nb/(beta*hbar)
                mb = sqrt(nb/(beta*m[i]))
                for k in range(nb):
                    if k == 0:
                        f = g[i]
                    else:
                        f = 4.*omega_n*sin(k*M_PI/nb)
                    c1 = exp(-0.5*dt*f)
                    c2 = sqrt(1-c1*c1)
                    for j in range(3):
                        if k == 0 or (nb % 2 == 0 and k == nb/2):
                            nmv[j, i+k] = c1*nmv[j, i+k] + c2*mb*MMTK.mtrand.standard_normal()
                        else:
                            nmv[j, i+k] = c1*nmv[j, i+k] + sqrt(0.5)*c2*mb*MMTK.mtrand.standard_normal()
                # Conversion back to Cartesian coordinates
                self.normalModeToCartesian(v, nmv, i, nb)

    cdef start(self):
        friction = self.getOption('centroid_friction')
        if isinstance(friction, ParticleProperties.ParticleScalar):
            self.gamma = friction.array
        else:
            assert isinstance(friction, numbers.Number)
            nbeads = self.universe.numberOfPoints()
            self.gamma = N.zeros((nbeads,), N.float)+friction
        PINormalModeIntegrator.start(self)

