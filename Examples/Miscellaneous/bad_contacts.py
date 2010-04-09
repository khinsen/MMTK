# Find bad contacts in a protein structure
#
# This script identifies that atoms in a protein structure
# that have the highest forces acting on them after a few
# steps of minimization.
#
# Note that there are no bad contacts in the example protein.
# The gradient norm values printed by this script can thus be
# used as a reference for what is normal.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from MMTK.Minimization import SteepestDescentMinimizer
from MMTK.Trajectory import StandardLogOutput

# Construct system
universe = InfiniteUniverse(Amber99ForceField())
universe.protein = Protein('bala1')

# Minimize
minimizer = SteepestDescentMinimizer(universe,
                                     actions=[StandardLogOutput(20)])
minimizer(convergence = 1.e-3, steps = 100)

# Calculate gradients and their lengths
energy, gradients = universe.energyAndGradients()
gradient_lengths = gradients.length()

# Make a list of all atoms sorted by decreasing gradient length
atoms = sorted(universe.atomList(),
               key = lambda a: gradient_lengths[a],
               reverse = True)

# Print the five atoms with the highest forces acting on them
for a in atoms[:5]:
    print a, gradient_lengths[a]

