# This example shows how to construct a periodic universe containing
# the crystallographic unit cell from the information in a PDB file.
#
# Note that this will not necessarily work with any PDB file. Many files
# use non-crystallographic symmetry information in a non-standard way.
# This is usually explained in REMARK records, but those cannot be
# evaluated automatically.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import Protein

# Read PDB configuration and create MMTK objects for all peptide chains.
# A C-alpha model is used to reduce the system size. You can remove
# 'model="calpha"' to get an all-atom model.
conf = PDBConfiguration('insulin.pdb')
chains = Collection(conf.createPeptideChains(model="calpha"))

# Copy and transform the objects representing the asymmetric unit in order
# to obtain the contents of the unit cell.
chains = conf.asuToUnitCell(chains)

# Construct a periodic universe representing the unit cell.
universe = conf.createUnitCellUniverse()

# Add each chain as one protein. If the unit cell contains multimers,
# the chains must be combined into protein objects by hand,
# as no corresponding information can be extracted from the PDB file.
for chain in chains:
    universe.addObject(Protein(chain))

# If VMD has been defined as the PDB viewer, this will not only show
# the molecules, but also the shape of the universe (which
# equals the unit cell).
universe.view()
