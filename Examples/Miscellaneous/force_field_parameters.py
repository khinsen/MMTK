# Obtain the force field parameters for a given universe.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField

# Define system
universe = InfiniteUniverse(Amber99ForceField())
universe.protein = Protein('bala1')

# Get energy term parameters
parameters = universe.energyEvaluatorParameters()

documentation = {
    "nonbonded":
    """
    excluded_pairs and one_four_pairs: lists of pairs of atom indices
    atom_subset: list of atom indices, or None for 'all atoms'
    """,
    "lennard_jones":
    """
    type: array of atom type numbers
    epsilon_sigma: 3D array. The first two indices are the atom type
                   numbers of the two atoms, the third index is 0
                   for epsilon and 1 for sigma
    type_names: the atom type names from the force field (for documentation)
    """,
    "electrostatic":
    """
    charge: array of partial charges
    """,
    "harmonic_distance_term":
    """
    a list of tuples (atom_index1, atom_index2, length, force_constant)
    """,
    "harmonic_angle_term":
    """
    a list of tuples (atom_index1, atom_index2, atom_index3, angle, force_constant)
    """,
    "cosine_dihedral_term":
    """
    a list of tuples (atom_index1, atom_index2, atom_index3, atom_index4,
                      multiplicity, angle, force_constant)
    """,
    }

# Print the parameters in a readable format
for p_type in parameters:
    p = parameters[p_type]
    print p_type
    print len(p_type)*'-'
    print documentation[p_type]
    if isinstance(p, list):
        print "["
        for data in p:
            print " ", data, ","
        print " ]"
    elif isinstance(p, dict):
        print "{"
        for key, value in p.items():
            s = repr(value)
            if len(key) + len(s) < 75:
                print " %s: %s," % (key, value)
            else:
                print " %s:" % key
                print "   %s," % value
        print "}"
    else:
        print p
    print
