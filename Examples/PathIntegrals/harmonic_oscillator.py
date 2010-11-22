# A one-dimensional harmonic oscillator made up of
# two atoms, an harmonic restraint, and a frozen subspace
#

from MMTK import *
from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint
from MMTK_PINormalModeIntegrator import PILangevinNormalModeIntegrator
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.Subspace import RigidMotionSubspace
from Scientific import N
from Scientific.Statistics import mean, standardDeviation

# Parameters
temperature = 30.*Units.K  # temperature
nb = 32                    # number of beads
k = 10.                    # force constant
dt = 0.005                 # time step for integration

# Define the system
universe = InfiniteUniverse()
universe.addObject(Environment.PathIntegrals(temperature))
universe.h1 = Atom('H', position = Vector(-0.3, 0., 0.))
universe.h1.setNumberOfBeads(nb)
universe.h2 = Atom('H', position = Vector(0.3, 0., 0.))
universe.h2.setNumberOfBeads(nb)
universe.setForceField(HarmonicDistanceRestraint(universe.h1,
                                                 universe.h2,
                                                 0., k))
universe.initializeVelocitiesToTemperature(temperature)

# Calculate the exact average quantum energy
m1 = universe.h1.mass()
m2 = universe.h2.mass()
m_eff = m1*m2/(m1+m2)
omega = N.sqrt(k/m_eff)
beta = 1./(Units.k_B*temperature)
e = Units.hbar*omega*(0.5+1./(N.exp(beta*Units.hbar*omega)-1))
print "Quantum energy:", e

#  Define the subspace of rigid-body motions of the universe
subspace = RigidMotionSubspace(universe, [universe])

# Create the integrator
integrator = PILangevinNormalModeIntegrator(universe, delta_t=dt,
                                            centroid_friction = 1./Units.ps,
                                            frozen_subspace = subspace)

# Equilibration
for i in range(10):
    integrator(steps=100)
    universe.initializeVelocitiesToTemperature(temperature)
integrator(steps=100)

# Production run
trajectory = Trajectory(universe, "harmonic_oscillator.nc", "w")
integrator(steps=10000,
           actions = [TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "auxiliary",
                                                    "configuration"),
                                       0, None, 5)])

# Print averages of the quantu energy estimator
data = trajectory.quantum_energy_primitive
print "Primitive estimator:", mean(data), "+/-",  standardDeviation(data)
data = trajectory.quantum_energy_virial
print "Virial estimator:", mean(data), "+/-",  standardDeviation(data)
data = trajectory.quantum_energy_centroid_virial
print "Centroid virial estimator:", mean(data), "+/-",  standardDeviation(data)

# Close the trajectory
trajectory.close()
