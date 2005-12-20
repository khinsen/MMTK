# This file contains some examples for non-standard generation of
# protein objects from PDB files.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain, Protein

#
# First problem: construct an all-atom model from a structure without
# hydrogens. This is the standard problem when using an all-atom force
# field with crystallographic structures.
#
# Note: the simple solution in this case is just
#       insulin = Protein('insulin.pdb')
# but the explicit form shown below is necessary when any kind of
# modification is required.
#
# Load the PDB file.
configuration = PDBConfiguration('insulin.pdb')

# Construct the peptide chain objects. This also constructs positions
# for any missing hydrogens, using geometrical criteria.
chains = configuration.createPeptideChains()

# Make the protein object.
insulin = Protein(chains)

# Write out the structure with hydrogens to a new file - we will use
# it as an input example later on.
insulin.writeToFile('insulin_with_h.pdb')


#
# Second problem: read a file with hydrogens and create a structure
# without them. This is useful for analysis; if you don't need the
# hydrogens, processing is faster without them. Or you might want
# to compare structures with and without hydrogens.
#
# Load the PDB file.
configuration = PDBConfiguration('./insulin_with_h.pdb')

# Delete hydrogens.
configuration.deleteHydrogens()

# Construct the peptide chain objects without hydrogens.
chains = configuration.createPeptideChains(model = 'no_hydrogens')

# Make the protein object
insulin = Protein(chains)

#
# Third problem: cut off three residues from the start of the second
# chain before constructing the protein. This is useful for comparing
# incomplete structures.
#
# Load the PDB file.
configuration = PDBConfiguration('insulin.pdb')

# Cut off the first three residues of the third chain.
configuration.peptide_chains[2].removeResidues(0, 3)

# Mark the second chain as modified
for chain in configuration.peptide_chains:
    chain.modified = 0
configuration.peptide_chains[2].modified = 1

# Construct the peptide chain objects without hydrogens. For
# modified chains, don't use the N-terminal form for the fist residue.
chains = []
for chain in configuration.peptide_chains:
    if chain.modified:
	chains.append(PeptideChain(chain, model='no_hydrogens', n_terminus=0))
    else:
	chains.append(PeptideChain(chain, model='no_hydrogens'))

# Make the protein object.
insulin = Protein(chains)


#
# Stop here, since the last "problem" is just an illustration,
# you can't run it directly because there is no "special residue"
# definition.
#
if 0:

    #
    # Fourth problem: construct a protein with a non-standard residue.
    # The new residue's name is 'XXX' and its definition is stored in the
    # MMTK group database under the name 'special_residue'.
    #
    from MMTK.Proteins import defineAminoAcidResidue
    defineAminoAcidResidue('special_residue', 'xxx')
    protein = Protein('special_protein.pdb')
