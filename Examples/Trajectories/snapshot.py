# "Manual" creation of a trajectory using the snapshot generator.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput

# Construct system
universe = InfiniteUniverse()
universe.protein = Protein('bala1')

# Transformation to be applied between the steps
transformation = Rotation(Vector(0.,0.,1.), 1.*Units.deg)

# Create trajectory
trajectory = Trajectory(universe, "rotation.nc", "w", "A rotating molecule")

# Create the snapshot generator
snapshot = SnapshotGenerator(universe,
                             actions = [TrajectoryOutput(trajectory,
                                                         ["all"], 0, None, 1)])

# Perform rotations and write the configurations
snapshot()
for i in range(20):
    universe.protein.applyTransformation(transformation)
    snapshot()

# Close trajectory
trajectory.close()
