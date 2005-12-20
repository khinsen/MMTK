# A simple Molecular Dynamics simulation using MPI.
#
# Note: the system simulated here is much too small to profit from
# parallelization, and multiple processors will actually slow down the
# simulation. This program is meant to be an illustration, not a
# real-life application.
#
# This example does exactly the same as the one in
# MolecularDynamics/protein.py. Compare the two to see what must be
# changed to use MPI parallelization.
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
from Scientific.MPI import world

# Define the system and write it to a file, on one processor only.
if world.rank == 0:
    universe = InfiniteUniverse(Amber94ForceField())
    universe.protein = Protein('bala1')
    universe.initializeVelocitiesToTemperature(50.*Units.K)
    save(universe, 'md_mpi.setup')

# Synchronize and load the universe definition on all processors.
world.barrier()
universe = load('md_mpi.setup')

# Note: Creating the universe in parallel on all processors is possible
# if the construction is purely deterministic, but with the assignment
# of random velocities, it cannot be guaranteed that all processors will
# end up with the same system state.


# Create integrator, specifying the MPI communicator to be used
# in energy evaluations
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs,
                                      mpi_communicator = world)

# Define output actions for one processor only
if world.rank == 0:
    # Log output to console every 100 steps for processor 0.
    output_actions = [StandardLogOutput(100)]
else:
    # Nothing for the other processors
    output_actions = []

# Heating and equilibration
integrator(steps=1000,
                    # Heat from 50 K to 300 K applying a temperature
                    # change of 0.5 K/fs; scale velocities at every step.
	   actions=[Heater(50.*Units.K, 300.*Units.K, 0.5*Units.K/Units.fs,
                           0, None, 1),
                    # Remove global translation every 50 steps.
		    TranslationRemover(0, None, 50),
                    # Remove global rotation every 50 steps.
		    RotationRemover(0, None, 50)] + output_actions)

# Synchronize
world.barrier()

# Define output actions for one processor only
if world.rank == 0:
    trajectory = Trajectory(universe, "bala1.nc", "w", "A simple test case")
    # Trajectory, restart, and log output for processor 0.
    output_actions = [RestartTrajectoryOutput("restart.nc", 5),
                      TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 2),
                      StandardLogOutput(10)]
else:
    # Nothing for the other processors.
    output_actions = []

# "Production" run
integrator(steps=100,
                      # Remove global translation every 50 steps.
           actions = [TranslationRemover(0, None, 50),
                      # Remove global rotation every 50 steps.
                      RotationRemover(0, None, 50)] + output_actions)

# Close trajectory
if world.rank == 0:
    trajectory.close()
