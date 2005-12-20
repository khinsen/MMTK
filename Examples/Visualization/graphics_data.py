# This example shows how the graphics generator functions in MMTK
# can be redirected to give access to the numerical values.
# It extracts the data about the arrows in an AtomicVectorField.
# The trick is to write a highly specialized graphics "module",
# which here is a fake module - in fact any object with the right
# attributes will do.

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.NormalModes import NormalModes
from MMTK.ForceFields import CalphaForceField
from MMTK.Field import AtomicVectorField

# Use generic color handling routines
import Scientific.Visualization.Color

# The fake graphics module - a class.
class DummyGraphics:

    def __init__(self):
        self.arrows = []

    def Color(self, rgb):
        return Scientific.Visualization.Color.Color(rgb)

    def ColorByName(self, name):
        return Scientific.Visualization.Color.ColorByName(name)

    def ColorScale(self, range):
        return Scientific.Visualization.Color.ColorScale(range)
    
    def SymmetricColorScale(self, range):
        return Scientific.Visualization.Color.SymmetricColorScale(range)

    def Sphere(self, center, radius, **attributes):
        pass

    def Cube(self, center, edge_length, **attributes):
        pass

    def Cylinder(self, point1, point2, radius, **attributes):
        pass

    def Cone(self, point1, point2, radius, **attributes):
        pass

    def Line(self, point1, point2, **attributes):
        pass

    # Only arrow data is stored for later extraction
    def Arrow(self, point1, point2, radius, **attributes):
        self.arrows.append((point1, point2))

    def Material(self, **attributes):
        return None

    def DiffuseMaterial(self, color):
        return None

    def EmissiveMaterial(self, color):
        return None

# Create the system
universe = InfiniteUniverse(CalphaForceField())
universe.protein = Protein('insulin', model='calpha')
# Calculate the normal modes
modes = NormalModes(universe)

# Generate the vector field for the first non-zero mode
vector_field = AtomicVectorField(universe, 0.5, modes[6])

# Generate the graphics data
graphics = DummyGraphics()
vector_field.graphicsObjects(graphics_module = graphics)

# Print the arrow coordinates
for point1, point2 in graphics.arrows:
    print point1, point2
