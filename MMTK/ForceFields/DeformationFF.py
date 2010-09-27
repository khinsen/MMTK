# Deformation force field
#
# Written by Konrad Hinsen
#

"""
Deformation force field (elastic network model)
For proteins, CalphaForceField is usually a better choice.
"""

__docformat__ = 'epytext'

from MMTK.ForceFields.ForceField import ForceField, ForceFieldData
from MMTK import Utility
from MMTK_forcefield import NonbondedList, NonbondedListTerm
from MMTK_deformation import DeformationTerm
from Scientific import N

#
# The deformation force field class
#
class DeformationForceField(ForceField):

    """
    Deformation force field for protein normal mode calculations

    The pair interaction force constant has the form
    k(r)=factor*exp(-(r**2-0.01)/range**2). The default value
    for |range| is appropriate for a C-alpha model of a protein.
    The offset of 0.01 is a convenience for defining factor;
    with this definition, factor is the force constant for a
    distance of 0.1nm.
    """

    def __init__(self, fc_length = 0.7, cutoff = 1.2, factor = 46402.):
        """
        @param fc_length: a range parameter
        @type fc_length: C{float}
        @param cutoff: the cutoff for pair interactions, should be
                       at least 2.5 nm. Pair interactions in periodic
                       systems are calculated using the minimum-image
                       convention; the cutoff should therefore never be
                       larger than half the smallest edge length of the
                       elementary cell.
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        """
        self.arguments = (fc_length, cutoff, factor)
        ForceField.__init__(self, 'deformation')
        self.arguments = (fc_length, cutoff, factor)
        self.fc_length = fc_length
        self.cutoff = cutoff
        self.factor = factor

    def ready(self, global_data):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        nothing = N.zeros((0, 2), N.Int)
        if subset1 is not None:
            set1 = set(a.index for a in subset1.atomList())
            set2 = set(a.index for a in subset2.atomList())
            excluded_pairs = set(Utility.orderedPairs(list(set1-set2))) \
                             | set(Utility.orderedPairs(list(set2-set1)))
            excluded_pairs = N.array(list(excluded_pairs))
            atom_subset = list(set1 | set2)
            atom_subset.sort()
            atom_subset = N.array(atom_subset)
        else:
            atom_subset = N.array([], N.Int)
            excluded_pairs = nothing
        nbl = NonbondedList(excluded_pairs, nothing, atom_subset,
                            global_data.get('nbeads'), global_data.get('bead_data'),
                            universe._spec, self.cutoff)
        update = NonbondedListTerm(nbl)
        cutoff = self.cutoff
        if cutoff is None:
            cutoff = 0.
        ev = DeformationTerm(universe._spec, nbl, self.fc_length,
                             self.cutoff, self.factor, 1.)
        return [update, ev]
