# A simple Molecular Dynamics simulation.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from MMTK.Dynamics import VelocityVerletIntegrator, Heater, \
                          TranslationRemover, RotationRemover
from MMTK.Visualization import view
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput, \
                            trajectoryInfo

# Define system
universe = InfiniteUniverse(Amber99ForceField(mod_files=['frcmod.ff99SB']))
universe.protein = Protein('bala1')

# Initialize velocities
universe.initializeVelocitiesToTemperature(50.*Units.K)
print 'Temperature: ', universe.temperature()
print 'Momentum: ', universe.momentum()
print 'Angular momentum: ', universe.angularMomentum()

# Create integrator
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs)

# Heating and equilibration
integrator(steps=1000,
                    # Heat from 50 K to 300 K applying a temperature
                    # change of 0.5 K/fs; scale velocities at every step.
	   actions=[Heater(50.*Units.K, 300.*Units.K, 0.5*Units.K/Units.fs,
                           0, None, 1),
                    # Remove global translation every 50 steps.
		    TranslationRemover(0, None, 50),
                    # Remove global rotation every 50 steps.
		    RotationRemover(0, None, 50),
                    # Log output to screen every 100 steps.
                    StandardLogOutput(100)])

# "Production" run
trajectory = Trajectory(universe, "bala1.nc", "w", "A simple test case")
integrator(steps=100,
                      # Remove global translation every 50 steps.
           actions = [TranslationRemover(0, None, 50),
                      # Remove global rotation every 50 steps.
                      RotationRemover(0, None, 50),
                      # Write every second step to the trajectory file.
                      TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 2),
                      # Write restart data every fifth step.
                      RestartTrajectoryOutput("restart.nc", 5),
                      # Log output to screen every 10 steps.
                      StandardLogOutput(10)])
trajectory.close()

# Print information about the trajectory file
print "Information about the trajectory file 'bala1.nc':"
print trajectoryInfo('bala1.nc')

# Reopen trajectory file
trajectory = Trajectory(universe, "bala1.nc", "r")

# Read step 10 and display configuration
step10 = trajectory[10]
view(universe, step10['configuration'])

# Print the kinetic energy along the trajectory
print "Kinetic energy along trajectory:"
print trajectory.kinetic_energy

trajectory.close()
