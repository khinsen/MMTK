# Modify sidechain dihedral angles
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
import Numeric; N = Numeric

# Construct system: lysozyme in vaccuum
universe = InfiniteUniverse()
universe.protein = Protein('~/hao/proteins/PDB/193l.pdb')

# Select residues to rotate
# (this particular choice here is completely arbitrary)
residues = [universe.protein[0][i] for i in [11, 35, 68, 110]]

# Create trajectory
trajectory = Trajectory(universe, "rotamers.nc", "w", "Sidechain rotations")

# Create the snapshot generator
snapshot = SnapshotGenerator(universe,
                             actions = [TrajectoryOutput(trajectory,
                                                         ["all"], 0, None, 1)])

# Perform sidechain rotations and write the configurations
snapshot()
for residue in residues:
    print residue
    chi = residue.chiAngle()
    for angle in N.arange(-N.pi, N.pi, N.pi/10.):
        chi.setValue(angle)
        print angle
        snapshot()

# Close trajectory
trajectory.close()
