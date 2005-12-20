# This example demonstrates some of the analysis features of MMTK.
# Two protein structures from the PDB are compared in some detail.
#
# Disclaimer: I know nothing about this protein, and the analysis
# I present here may be of no biological interest.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain

# First we read the two PDB files.

configuration1 = PDBConfiguration('4q21.pdb.gz')
configuration2 = PDBConfiguration('6q21.pdb.gz')

# The first file contains a monomer, the second one a tetramer in
# which each chain is almost identical to the monomer from the first
# file. We have to cut off the last (incomplete) residue from the
# monomer and the last three residues of each chain of the tetramers
# to get matching sequences. We'll just deal with one of the chains of
# the tetramer here.

monomer = configuration1.peptide_chains[0]
monomer.removeResidues(-1, None)
tetramer_1 = configuration2.peptide_chains[0]
tetramer_1.removeResidues(-3, None)

# Next we build a PeptideChain for the monomer. We have to specify
# that we have no C terminus (missing) and no hydrogens.

chain = PeptideChain(monomer, model='no_hydrogens', c_terminus=0)

# How big is this beast? For a quick guess, print the size of the smallest
# rectangular box containing it:

p1, p2 = chain.boundingBox()
print "Size of a bounding box: %4.2f x %4.2f x %4.2f nm" % tuple(p2-p1)

# Formalities: we define a universe and put the chain inside.
# Then we make a copy of the configuration of the universe for later use.
# This is necessary because configurations are defined only for
# universes (for various good reasons...).

universe = InfiniteUniverse()
universe.addObject(chain)
monomer_configuration = copy(universe.configuration())

# Next we change the conformation of the protein to that of the first
# tetramer chain.

tetramer_1.applyTo(chain)

# Now we want to get rid of the difference in position and orientation
# of the two configurations, so we calculate the transformation that
# minimizes the RMS difference.

transformation, rms = chain.findTransformation(monomer_configuration)

# That's the transformation *from* the current *to* the monomer
# configuration. Let's analyze it:

print "Translation/rotation fit:"
print "  RMS difference:", rms

translation = transformation.translation()
rotation = transformation.rotation()
axis, angle = rotation.axisAndAngle()

print "  Rotation axis: ", axis
print "  Rotation angle: ", angle
print "  Translation: ", translation.displacement()

# Note that order is important: the rotation (around the origin) is
# applied first, and then the translation. If you want to do it the
# other way round, the translation displacement must be different.

# Note also the units: distances in nanometers, angles in radians!
# If you insist on Angstrom and degrees, you'll have to type a bit more:

print "  Rotation angle in degrees: ", angle/Units.deg
print "  Translation in Angstrom: ", translation.displacement()/Units.Ang

# Just for fun we'll do it again, but minimizing the RMS for the
# C_alpha atoms only:

c_alphas = chain.residues().map(lambda residue: residue.peptide.C_alpha)
calpha_transformation, rms = c_alphas.findTransformation(monomer_configuration)

print "C_alpha fit:"
print "  RMS difference:", rms

translation = calpha_transformation.translation()
rotation = calpha_transformation.rotation()
axis, angle = rotation.axisAndAngle()

print "  Rotation axis: ", axis
print "  Rotation angle: ", angle
print "  Translation: ", translation.displacement()

# As expected there is little difference. Now we apply the transformation
# to the current configuration to align it with the reference configuration:

chain.applyTransformation(transformation)

# Let's calculate the local deformation due to the change in conformation
# for each atom. Small values indicate rigid regions and large values
# indicate flexible regions. The visualization requires a VRML viewer.
# Reference: K. Hinsen, A. Thomas, M.J. Field:  Proteins 34, 369 (1999)

from MMTK.Deformation import FiniteDeformationFunction
from Scientific.Visualization import VRML
difference = monomer_configuration-universe.configuration()
deformation_function = FiniteDeformationFunction(universe, 0.3, 0.6)
deformation = deformation_function(difference)
graphics = universe.graphicsObjects(color_values = deformation,
                                    graphics_module = VRML)
VRML.Scene(graphics).view()

# The deformation shows a small rigid region and a much larger
# more flexible region.

# Action: We show an animation of the two configurations. This requires
# an animation-capable external viewer, e.g. VMD.

from MMTK.Visualization import viewSequence
viewSequence(chain, [universe.configuration(), monomer_configuration])

# Interesting: there's a floppy loop with residues Gln61, Glu62, and Glu63
# at the tip. It seems to do a funny rotation relative to a closeby helix
# (His94-Lys101) when switching between the two configurations. Let's
# find out what this rotation is!

# First, we define the parts of the chain we are interested in. We don't
# want the sidechains to mess up the picture, so we use just the backbone:

tip = chain[60:64].backbone()
helix = chain[93:101].backbone()

# Then we align the two configurations to make the helix match as well
# as possible:

transformation, rms = helix.findTransformation(monomer_configuration)
chain.applyTransformation(transformation)

# Note that we *fit* only the helix, but *apply* the transformation
# to the whole chain.

# And now another fit for the tip of the loop. This time we don't want
# to apply the transformation, but print out its characteristics:

transformation, rms = tip.findTransformation(monomer_configuration)

print "Movement of the tip of the loop:"
print "  RMS difference:", rms

translation = transformation.translation()
rotation = transformation.rotation()
axis, angle = rotation.axisAndAngle()

print "  Rotation axis: ", axis
print "  Rotation angle: ", angle
print "  Translation: ", translation.displacement()
