# Continuation of the simulation from protein.py using the restart
# trajectory.
#
# Note: Restart trajectories don't differ much from ordinary
# trajectories; they have a fixed number of steps that are reused
# cyclically, but otherwise they are plain trajectory files.
# As a consequence, you can restart a molecular dynamics run
# from any trajectory that contains positions and velocities,
# with only a minor modification mentioned below.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.Dynamics import VelocityVerletIntegrator, Heater, \
                          TranslationRemover, RotationRemover
from MMTK.Visualization import view
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput, \
                            trajectoryInfo

# Define system
universe = InfiniteUniverse(Amber94ForceField())
universe.protein = Protein('bala1')

# Initialize system state from the restart trajectory
universe.setFromTrajectory(Trajectory(universe, "restart.nc"))

# Note: for restarting from a standard trajectory, use:
#
#   universe.setFromTrajectory(Trajectory(universe, "trajectory.nc"), -1)
#
# in order to load the last step (default would be the first one).

# Create integrator
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs)

# Do some more MD steps.
trajectory = Trajectory(universe, "bala1.nc", "a")
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
