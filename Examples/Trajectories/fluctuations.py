# Calculate atomic fluctuations from a trajectory.
#
# The trajectory file contains a protein. The fluctuations are calculated
# for all atoms, and then the average over the C-alpha atoms is determined.
# Note that calculating the fluctuations for only the C-alpha atoms is
# more complicated and no faster.
#
# This example illustrates:
# 1) Reading trajectory files
# 2) Selecting parts of a system
# 3) Calculating trajectory averages
#

from MMTK import *
from MMTK.Trajectory import Trajectory
from MMTK.Proteins import Protein

# Open the trajectory, use every tenth step
trajectory = Trajectory(None, 'lysozyme.nc')[::10]
universe = trajectory.universe

# Calculate the average conformation
average = ParticleVector(universe)
for step in trajectory:
    average += step['configuration']
average /= len(trajectory)

# Calculate the fluctiations for all atoms
fluctuations = ParticleScalar(universe)
for step in trajectory:
    d = step['configuration']-average
    fluctuations += d*d
fluctuations /= len(trajectory)

# Collect the C-alpha atoms
calphas = Collection()
for protein in universe.objectList(Protein):
    for chain in protein:
        for residue in chain:
            calphas.addObject(residue.peptide.C_alpha)

# Average over the C-alpha atoms
calpha_mask = calphas.booleanMask()
print (fluctuations*calpha_mask).sumOverParticles()/calphas.numberOfAtoms()
