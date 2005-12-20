# This program reads a trajecory and calculates the average
# configuration, which is then visualized. This is not a particularly
# useful operation, but it illustrates the use of trajectories and
# configuration variables.
#
# Note: In order to create the trajectory file "bala1.nc" which is
#       read by this example, you must run the example
#       MolecularDynamics/protein.py!
#

from MMTK import *
from MMTK.Trajectory import Trajectory
from MMTK.Visualization import view

# Open the input trajectory.
trajectory = Trajectory(None, 'bala1.nc')
universe = trajectory.universe

# Create a zeroed particle vector object.
sum = ParticleVector(universe)

# Loop over all steps and add the configurations to the sum.
for configuration in trajectory.configuration:
    sum = sum + configuration

# Divide by the number of steps.
sum = sum/len(trajectory)

# Visualize the average.
view(universe, sum)
