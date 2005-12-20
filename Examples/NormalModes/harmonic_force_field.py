# Normal mode calculation using a simplified harmonic force field.
# For a description of the force field, see
#
#     K. Hinsen & G.R. Kneller
#     A simplified force field for describing vibrational protein dynamics
#     over the whole frequency range
#     J. Chem. Phys. 111, 10766 (1999)
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import HarmonicForceField
from MMTK.NormalModes import VibrationalModes
from MMTK.Visualization import view

# Construct system
universe = InfiniteUniverse(HarmonicForceField())
universe.protein = Protein('bala1')

# Calculate normal modes
modes = VibrationalModes(universe)

# Print frequencies
for mode in modes:
    print mode

# Show animation of the first non-trivial mode 
view(modes[6])
