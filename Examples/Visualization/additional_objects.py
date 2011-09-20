# Generate a protein backbone representation plus an indication
# of the principal axes of inertia by arrows.
#
from MMTK import *
from MMTK.Proteins import Protein
from Scientific import N, LA

# Import the graphics module. Substitute any other graphics
# module name to make the example use that module.
from Scientific.Visualization import VRML2; module = VRML2

# Create the protein and find its center of mass and tensor of inertia.
protein = Protein('insulin')
center, inertia = protein.centerAndMomentOfInertia()

# Diagonalize the inertia tensor and scale the axes to a suitable length.
mass = protein.mass()
diagonal, directions = LA.eigenvectors(inertia.array)
diagonal = N.sqrt(diagonal/mass)

# Generate the backbone graphics objects.
graphics = protein.graphicsObjects(graphics_module = module,
                                   model = 'backbone', color = 'red')

# Add an atomic wireframe representation of all valine sidechains
valines = protein.residuesOfType('val')
sidechains = valines.map(lambda r: r.sidechain)
graphics = graphics + sidechains.graphicsObjects(graphics_module = module,
                                                 model='wireframe',
                                                 color = 'blue')

# Add three arrows corresponding to the principal axes of inertia.
for length, axis in map(None, diagonal, directions):
    graphics.append(module.Arrow(center, center+length*Vector(axis), 0.02,
                                 material=module.EmissiveMaterial('green')))

# Display everything using a VRML browser.
scene = module.Scene(graphics)
scene.view()
