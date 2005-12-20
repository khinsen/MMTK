# This program reads a trajecory containing one or more proteins
# and extracts a trajectory for a C-alpha model, which it stores in
# a new trajectory. This trajectory is smaller and easier to analyze
# or visualize. Moreover, the global (rigid-body) motion of the
# protein(s) is eliminated from the trajectory.
#
# Note: You cannot run this example without having a suitable
#       trajectory file, whose name you should put in place of
#       "full_trajectory.nc" below.
#

from MMTK import *
from MMTK.Trajectory import Trajectory, TrajectoryOutput, SnapshotGenerator
from MMTK.Proteins import Protein, PeptideChain
import Numeric

# Open the input trajectory.
trajectory = Trajectory(None, 'full_trajectory.nc')
universe = trajectory.universe

# Choose all objects of class "Protein"
proteins = universe.objectList(Protein)

# Construct an initial configuration in which the proteins are contiguous.
conf = universe.contiguousObjectConfiguration(proteins)

# Construct a C_alpha representation and a mapping from the "real"
# C_alpha atoms to their counterparts in the reduced model.
universe_calpha = InfiniteUniverse()
index = 0
map = []
for protein in proteins:
    chains_calpha = []
    for chain in protein:
        chains_calpha.append(PeptideChain(chain.sequence(),
                                          model='calpha'))
    protein_calpha = Protein(chains_calpha)
    universe_calpha.addObject(protein_calpha)
    for i in range(len(protein)):
        chain = protein[i]
        chain_calpha = protein_calpha[i]
        for j in range(len(chain)):
            r = conf[chain[j].peptide.C_alpha]
            chain_calpha[j].peptide.C_alpha.setPosition(r)
            chain_calpha[j].peptide.C_alpha.index = index
            map.append(chain[j].peptide.C_alpha.index)
            index = index + 1
universe_calpha.configuration()

# Open the new trajectory for just the interesting objects.
trajectory_calpha = Trajectory(universe_calpha, 'calpha_trajectory.nc', 'w')

# Make a snapshot generator for saving.
snapshot = SnapshotGenerator(universe_calpha,
                             actions = [TrajectoryOutput(trajectory_calpha,
                                                         None, 0, None, 1)])

# Loop over all steps, eliminate rigid-body motion with reference to
# the initial configuration, and save the configurations to the new
# trajectory.
first = 1
for step in trajectory:
    conf = universe.contiguousObjectConfiguration(proteins,
                                                  step['configuration'])
    conf_calpha = Configuration(universe_calpha, Numeric.take(conf.array, map),
                                None)
    universe_calpha.setConfiguration(conf_calpha)
    if first:
        initial_conf_calpha = copy(conf_calpha)
        first = 0
    else:
        tr, rms = universe_calpha.findTransformation(initial_conf_calpha)
        universe_calpha.applyTransformation(tr)
    del step['step']
    del step['configuration']
    try:
        del step['box_size']
    except KeyError: pass
    try:
        del step['velocities']
    except KeyError: pass
    snapshot(data = step)

# Close both trajectories.
trajectory.close()
trajectory_calpha.close()
