# This example shows how to calculate approximate low-frequency
# modes for big proteins. For a description of the techniques,
# see
#
#    K. Hinsen
#    Analysis of domain motions by approximate normal mode calculations
#    Proteins 33, 417 (1998)
#
# and
#
#    K. Hinsen, A.J. Petrescu, S. Dellerue, M.C. Bellissent-Funel, G.R. Kneller
#    Harmonicity in slow protein dynamics
#    Chem. Phys. 261, 25 (2000)
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.FourierBasis import FourierBasis, estimateCutoff
from MMTK.NormalModes import EnergeticModes
from MMTK.Visualization import view

# Construct system
universe = InfiniteUniverse(CalphaForceField(2.5))
universe.protein = Protein('insulin.pdb', model='calpha')

# Find a reasonable basis set size and cutoff
nbasis = max(10, universe.numberOfAtoms()/5)
cutoff, nbasis = estimateCutoff(universe, nbasis)
print "Calculating %d low-frequency modes." % nbasis

if cutoff is None:
    # Do full normal mode calculation
    modes = EnergeticModes(universe, 300.*Units.K)
else:
    # Do subspace mode calculation with Fourier basis
    subspace = FourierBasis(universe, cutoff)
    modes = EnergeticModes(universe, 300.*Units.K, subspace)

# Show animation of the first non-trivial mode
view(modes[6], 15.)
