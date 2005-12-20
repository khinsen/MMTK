# Example: Vector field analysis of the slowest normal mode of insulin.
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
from MMTK.ForceFields import DeformationForceField
from MMTK.NormalModes import NormalModes
from MMTK.Field import AtomicVectorField
from Scientific.Visualization import VRML2 as VRML

# Construct system
universe = InfiniteUniverse(DeformationForceField())
universe.protein = Protein('insulin.pdb', model='calpha')

# Calculate normal modes
modes = NormalModes(universe)

# Create the vector field corresponding to the first non-zero mode
field = AtomicVectorField(universe, 0.5, modes[6])

# Calculate the curl of the vector field. The curl indicated the
# local rotational motion; the direction of the curl vector is
# the rotation axis, and its length is proportional to the amount
# of the rotation.
curl = field.curl()

# Create graphics objects for the protein and the vector fields
graphics = universe.graphicsObjects(model='backbone', graphics_module=VRML) + \
           field.graphicsObjects(color='yellow',
                                 scale=80., range=(0., 0.01),
                                 graphics_module=VRML) + \
           curl.graphicsObjects(color='red',
                                 scale=80., range=(0., 0.01),
                                 graphics_module=VRML)

# Show graphics in a VRML browser
VRML.Scene(graphics).view()
