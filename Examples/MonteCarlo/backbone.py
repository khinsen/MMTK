# Generate an ensemble of backbone configurations for a protein using
# a simplified protein model which consists only of the C_alpha atoms.
# The ensemble is written to a trajectory file.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.NormalModes import NormalModes
from MMTK.Random import gaussian
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput

# Construct the system. The scaling factor of 0.1 for the CalphaForceField
# was determined empirically for a C-phycocyanin dimer at 300 K; its value
# is expected to depend at least on the protein and the temperature,
# but the order of magnitude should be right for mid-sized proteins.
universe = InfiniteUniverse(CalphaForceField(2.5, 0.1))
universe.protein = Protein('insulin.pdb', model='calpha')

# Set all masses to the same value, since we don't want
# mass-weighted sampling.
for atom in universe.atomList():
    atom.setMass(1.)

# Calculate normal modes. Note that the normal modes are automatically
# scaled to their vibrational amplitudes at a given temperature.
modes = NormalModes(universe, 300.*Units.K)

# Create trajectory
trajectory = Trajectory(universe, "insulin_backbone.nc", "w",
                        "Monte-Carlo sampling for insulin backbone")

# Create the snapshot generator
snapshot = SnapshotGenerator(universe,
                             actions = [TrajectoryOutput(trajectory,
                                                         ["all"], 0, None, 1)])

# Generate an ensemble of 100 configurations. The scaling factor for
# each mode is half the vibrational amplitude, which explains the
# factor 0.5.
minimum = copy(universe.configuration())
for i in range(100):
    conf = minimum
    for j in range(6, len(modes)):
        conf = conf + gaussian(0., 0.5)*modes[j]
    universe.setConfiguration(conf)
    snapshot()

# Close trajectory
trajectory.close()
