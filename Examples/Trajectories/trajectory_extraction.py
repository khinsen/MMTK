# This program reads a trajecory for a solvated protein and
# extracts just the protein data, which it stores in a new
# trajectory. This trajectory is smaller and easier to analyze
# or visualize.
#
# Note: You cannot run this example without having a suitable
#       trajectory file, whose name you should put in place of
#       "full_trajectory.nc" below.
#

from MMTK import *
from MMTK.Trajectory import Trajectory, TrajectoryOutput, SnapshotGenerator

# Open the input trajectory.
full = Trajectory(None, 'full_trajectory.nc')

# Collect the items you want to keep. Here we keep everything but
# water molecules.
universe = full.universe
keep = Collection()
for object in universe:
    try:
        is_water = object.type.name == 'water'
    except AttributeError:
        is_water = 0
    if not is_water:
        keep.addObject(object)

# Open the new trajectory for just the interesting objects.
subset = Trajectory(keep, 'subset_trajectory.nc', 'w')

# Make a snapshot generator for saving.
snapshot = SnapshotGenerator(universe,
                             actions = [TrajectoryOutput(subset, None,
                                                         0, None, 1)])

# Loop over all steps and save them to the new trajectory.
for configuration in full.configuration:
    universe.setConfiguration(configuration)
    snapshot()

# Close both trajectories.
full.close()
subset.close()
