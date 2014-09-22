# A sugar molecule
#

from MMTK import *
from MMTK.ForceFields import Amber99ForceField
from MMTK.PICartesianIntegrator import PILangevinCartesianIntegrator
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.NormalModes import VibrationalModes

from Scientific import N
from Scientific.Statistics import mean, standardDeviation

# Add the Database directory for this example to the database search path
import MMTK.Database
import os
MMTK.Database.path.append(os.path.abspath('./Database'))

# A function that estimates the integration time step
# for a universe with path integrals.
def timeStep(universe):
    # We want to choose the time step based on the frequency of the
    # fastest vibrations of the classical analog system, i.e. the
    # beads with harmonic interactions. This requires an explicit
    # representation of the spring terms in the potential.
    modes = VibrationalModes(universe)

    # Make timestep 1/20 of the period of the fastest vibration
    return 0.05/modes[-1].frequency

# Simulation parameters
temperature = 200.*Units.K
number_of_beads = 32

# Make a universe using the Amber99 force field with a mod file for sugars
universe = InfiniteUniverse(Amber99ForceField(mod_files=['glycam_06g.dat']))

# Add a sugar molecule
sugar = Molecule('beta-D-arabinose-OMe')
universe.sugar = sugar

# The PathIntegrals environment object defines the temperature
# used in the spring terms between the beads of each atom.
# The thermostat in the PILangevinCartesianIntegrator uses
# the same temperature.
universe.addObject(Environment.PathIntegrals(temperature, True))

# Use path integrals with 32 beads for all atoms.
for atom in universe.atomIterator():
    atom.setNumberOfBeads(number_of_beads)

# Compute the integration time step.
dt = timeStep(universe)

# Initialize the velocities from a Gaussian distribution.
universe.initializeVelocitiesToTemperature(temperature)

# Create the integrator
integrator = PILangevinCartesianIntegrator(universe, delta_t=dt,
                                           friction_matrix = N.array([[1./Units.ps]]))

# Produce a trajectory.
trajectory = Trajectory(universe, "sugar.nc", "w")
integrator(steps=100,
           actions = [TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "auxiliary",
                                                    "configuration"),
                                       0, None, 5)])

# Print averages of the quantum energy estimator
data = trajectory.quantum_energy_primitive
print "Primitive estimator:", mean(data), "+/-",  standardDeviation(data)
data = trajectory.quantum_energy_virial
print "Virial estimator:", mean(data), "+/-",  standardDeviation(data)
data = trajectory.quantum_energy_centroid_virial
print "Centroid virial estimator:", mean(data), "+/-",  standardDeviation(data)
data = trajectory.temperature
print "Temperature:", mean(data), "+/-",  standardDeviation(data)

# Close the trajectory
trajectory.close()


