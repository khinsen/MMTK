# Finding the partial charges of the atoms in a universe
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber12SBForceField

# Define system
universe = InfiniteUniverse(Amber12SBForceField())
universe.protein = Protein('bala1')

# Print positions and charges of all atoms
charges = universe.charges()
for a in universe.atomIterator():
    print a.position(), charges[a]
