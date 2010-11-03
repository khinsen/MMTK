# Example: Vector field visualization of the slowest normal mode of insulin.
#
# This example shows the use of vector fields in the analysis and
# visualization of collective motions. For a more elaborate application,
# see
#       A. Thomas, K. Hinsen, M.J. Field, D. Perahia
#       Tertiary and quaternary confrmational changes in aspartate
#       transcarbamylase: a normal mode study.
#       Proteins, 34(1), 96-112 (1999)
#
from MMTK import *
from MMTK.Proteins import Protein
from MMTK.NormalModes import EnergeticModes
from MMTK.ForceFields import CalphaForceField
from MMTK.Field import AtomicVectorField
from Scientific.Visualization import VMD

# Create a simulation universe containing the protein
universe = InfiniteUniverse(CalphaForceField())
universe.addObject(Protein("insulin.pdb", model='calpha'))

# Calculate normal modes
modes = EnergeticModes(universe)

# Create the vector field corresponding to the first non-zero mode
# using a grid spacing of 0.5 nm.
field = AtomicVectorField(universe, 0.5, modes[6])

# Create graphics objects for the vector fields.
# The arrows are yellow and their lengths are multiplied by 80.
# Only displacement vectors of lengths between 0. and 0.01 will be shown.
# (This restriction is not necessary for this example, but used for
#  illustration.)
graphics = field.graphicsObjects(color='yellow',
                                 scale=80., range=(0., 0.01),
                                 graphics_module=VMD)

# Create graphics object for the protein.
graphics = graphics + [VMD.Molecules(universe)]

# Launch VMD.
VMD.Scene(graphics, scale=1./Units.Ang).view()
