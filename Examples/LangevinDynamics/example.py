# An example for LangevinDynamics
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.Dynamics import TranslationRemover, RotationRemover, VelocityScaler
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
from LangevinDynamics import LangevinIntegrator

# Define system
universe = InfiniteUniverse(Amber94ForceField())
universe.protein = Protein('bala1')
temperature = 300.*Units.K

# Initialize velocities
universe.initializeVelocitiesToTemperature(temperature)

# Set friction coefficients
friction = universe.masses()*0.01/Units.ps

# Create integrator
integrator = LangevinIntegrator(universe, delta_t=1.*Units.fs,
                                friction=friction, temperature=temperature)

# Equilibration
integrator(steps=10000,
	   actions=[# Scale velocities every 50 steps.
                    VelocityScaler(temperature, 0.1*temperature,
                                   0, None, 50),
                    # Remove global translation every 50 steps.
		    TranslationRemover(0, None, 50),
                    # Remove global rotation every 50 steps.
		    RotationRemover(0, None, 50),
                    # Log output to screen every 500 steps.
                    StandardLogOutput(500)])

# "Production" run
trajectory = Trajectory(universe, "langevin.nc", "w", "Langevin test")
integrator(steps=10000,
           actions = [# Remove global translation every 50 steps.
                      TranslationRemover(0, None, 50),
                      # Remove global rotation every 50 steps.
                      RotationRemover(0, None, 50),
                      # Write every fifth step to the trajectory file.
                      TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 10),
                      # Log output to screen every 100 steps.
                      StandardLogOutput(100)])
trajectory.close()
