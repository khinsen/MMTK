# Example: water molecules on a lattice
#
# This example creates a system of water molecules whose centers
# of mass are positioned simple cubic lattice consisting of
# 3x4x5 cells, i.e. there are 60 water molecules in total.
# The molecules are put into a suitably sized periodic universe.
#

from MMTK import *
from Scientific.Geometry.Objects3D import SCLattice
from MMTK.Visualization import view

# Define parameters of the lattice
edge_length = 0.5*Units.nm
lattice_size = (3, 4, 5)

# Construct the universe
universe = OrthorhombicPeriodicUniverse((edge_length*lattice_size[0],
                                         edge_length*lattice_size[1],
                                         edge_length*lattice_size[2]))

# Add the water molecules
for point in SCLattice(edge_length, lattice_size):
    universe.addObject(Molecule('water', position = point))

# Visualization
view(universe)
