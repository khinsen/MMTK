# This example shows how to construct a periodic universe containing
# the crystallographic unit cell from the information in a PDB file.
#
# Note that this will not necessarily work with any PDB file. Many files
# use non-crystallographic symmetry information in a non-standard way.
# This is usually explained in REMARK records, but those cannot be
# evauluated automatically.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import Protein

# A utility function that creates an image of an object by making
# a copy and applying a transformation to the copy.
def makeImage(object, transformation):
    image = deepcopy(object) 
    for atom in image.atomList():
        atom.setPosition(transformation(atom.position()))
    return image

# Read PDB configuration and create MMTK objects for all peptide chains.
# A C-alpha model is used to reduce the system size. You can remove
# 'model="calpha"' to get an all-atom model, but for insulin this will
# create more than 380000 atoms for the 27-unit-cell crystal!
conf = PDBConfiguration('insulin.pdb')
chains = Collection(conf.createPeptideChains(model="calpha"))

# Apply non-crystallographic symmetries to construct the asymmetric unit
asu = Collection(chains)
for so in conf.ncs_transformations:
    if not so.given:
        image = makeImage(chains, so)
        asu.addObject(image)

# Apply crystallographic symmetries to construct the unit cell
# Note that the list of crystallographic symmetries includes the
# identity transformation, so the unmodified asu is not added
# to the unit cell.
cell = Collection()
for so in conf.cs_transformations:
    image = makeImage(asu, so)
    cell.addObject(image)
    # Translate the image such that its center of mass is
    # located inside the unit cell.
    cm = image.centerOfMass()
    cm_fract = conf.to_fractional(cm)
    cm_fract = Vector(cm_fract[0] % 1., cm_fract[1] % 1., cm_fract[2] % 1.) \
               - Vector(0.5, 0.5, 0.5)
    cm = conf.from_fractional(cm_fract)
    image.translateTo(cm)

# Make a universe corresponding to the unit cell
e1 = conf.from_fractional(Vector(1., 0., 0.))
e2 = conf.from_fractional(Vector(0., 1., 0.))
e3 = conf.from_fractional(Vector(0., 0., 1.))

# Use OrthorhombicPeriodicUniverse when possible as it is more efficient
if abs(e1.normal()*Vector(1., 0., 0.)-1.) < 1.e-15 \
       and abs(e2.normal()*Vector(0., 1., 0.)-1.) < 1.e-15 \
       and abs(e3.normal()*Vector(0., 0., 1.)-1.) < 1.e-15:
    universe = OrthorhombicPeriodicUniverse((e1.length(),
                                             e2.length(),
                                             e3.length()))
else:
    universe = ParallelepipedicPeriodicUniverse((e1, e2, e3))

# Add each chain as one protein. If the unit cell contains multimers,
# the chains must be combined into protein objects by hand,
# as no corresponding information can be extracted from the PDB file.
for chain in cell:
    universe.addObject(Protein(chain))

# If VMD has been defined as the PDB viewer, this will not only show
# the molecules, but also the shape of the simulation box (which
# equals the unit cell).
universe.view()
