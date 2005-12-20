# A normal mode calculation in the subspace of residue rigid-body motion.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.NormalModes import VibrationalModes
from MMTK.Subspace import RigidMotionSubspace
from MMTK.Minimization import ConjugateGradientMinimizer
from MMTK.Trajectory import StandardLogOutput

import Numeric

# Construct system
universe = InfiniteUniverse(Amber94ForceField())
universe.protein = Protein('bala1')

# Minimize
minimizer = ConjugateGradientMinimizer(universe,
                                       actions=[StandardLogOutput(50)])
minimizer(convergence = 1.e-3, steps = 10000)

# Set up the subspace: rigid-body translation and rotation for each residue
subspace = RigidMotionSubspace(universe, universe.protein.residues())

# Calculate normal modes
modes = VibrationalModes(universe, subspace=subspace)

# Calculate full  modes for comparison
full_modes = VibrationalModes(universe, subspace=None)

# Compare the modes. For each reduced mode, find the full mode with
# the most similar displacement vector. Print both frequencies and
# the overlap of the displacement vectors.
masses = universe.masses()
for i in range(6, len(modes)):
    m1 = modes[i]
    overlap = []
    for m2 in full_modes:
	o = m1.massWeightedDotProduct(m2) / \
	    Numeric.sqrt(m1.massWeightedDotProduct(m1)) / \
	    Numeric.sqrt(m2.massWeightedDotProduct(m2))
	overlap.append(abs(o))
    best = Numeric.argsort(overlap)[-1]
    print m1.frequency, full_modes[best].frequency, overlap[best]
