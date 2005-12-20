# This example shows how to calculate Brownian modes
# for big proteins. For a description of the techniques, see
#
#    K. Hinsen, A.J. Petrescu, S. Dellerue, M.C. Bellissent-Funel, G.R. Kneller
#    Harmonicity in slow protein dynamics
#    Chem. Phys. 261, 25 (2000)
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.ProteinFriction import calphaFrictionConstants
from MMTK.NormalModes import BrownianModes

# Construct system
universe = InfiniteUniverse(CalphaForceField(2.5))
universe.protein = Protein('insulin.pdb', model='calpha')

# Calculate Brownian modes
modes = BrownianModes(universe, calphaFrictionConstants(universe.protein),
                      300.*Units.K)

# Print relaxation times for the slowest ten modes
#
# Note: realistic time scales can be expected only if the potential
# has been scaled by a protein-specific factor such that the fluctuation
# amplitudes are approximated well. That factor must in practice be
# obtained from experiment. See the publication cited above for details.
#
for i in range(6, 16):
    print i, 1./modes[i].inv_relaxation_time

