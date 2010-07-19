# Set up an all-atom model and a C-alpha model for the same protein,
# with a facility for copying configurations back and forth.
#

from MMTK import *
from MMTK.Proteins import Protein

# Create an all-atom universe and a c-alpha universe
universe_aa = InfiniteUniverse()
universe_ca = InfiniteUniverse()

# Add the protein to both universes
protein_aa = Protein('insulin.pdb')
universe_aa.addObject(protein_aa)
protein_ca = Protein('insulin.pdb', model='calpha')
universe_ca.addObject(protein_ca)

# Set up a dictionary mapping aa atoms to ca atoms
atom_map = {}
for ichain in range(len(protein_aa)):
    chain_aa = protein_aa[ichain]
    chain_ca = protein_ca[ichain]
    for iresidue in range(len(chain_aa)):
        ca_aa = chain_aa[iresidue].peptide.C_alpha
        ca_ca = chain_ca[iresidue].peptide.C_alpha
        atom_map[ca_aa] = ca_ca

# Copy the aa configuration to the ca configuration
for ca_aa, ca_ca in atom_map.items():
    ca_ca.setPosition(ca_aa.position())
