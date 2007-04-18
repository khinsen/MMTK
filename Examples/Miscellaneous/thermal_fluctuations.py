# Retrieve atomic fluctuation information from a PDBConfiguration object.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from Scientific import N

# Read a PDB file containing ANISOU records.
conf = PDBConfiguration('1G66.pdb')

# By passing the applyTo methods of MMTK.PDB a dictionary argument
# atom_map, one obtains a dictionary from MMTK atom objects to the
# corresponding PDB atom objects. This dictionary can be used to
# retrieve additional atom data from the PDB file.
atom_map = {}
for c in conf.peptide_chains:
    # Create a PeptideChain object
    chain = c.createPeptideChain()
    # Retrieve the atom_map dictionary
    # Note: this also redefines the configuration, but in this
    # application that makes no difference since it is the same
    # that was defined in the previous line.
    c.applyTo(chain, atom_map=atom_map)

    # Print the B factor and the trace of the anisotropic displacement
    # tensor for each atom for which both are available. They should
    # be equal, but due to the limited precision in PDB files
    # there can be small differences.
    for atom in chain.atomList():
        try:
            pdb_atom = atom_map[atom]
        except KeyError:
            print atom, " not in PDB file"
            continue
        b_factor = pdb_atom.properties['temperature_factor'] * Units.Ang**2
        try:
            u_trace = pdb_atom.properties['u'].trace()/3. * Units.Ang**2
        except KeyError:
            u_trace = "[no ANISOU record]"
        print atom, b_factor/(8*N.pi**2), u_trace
