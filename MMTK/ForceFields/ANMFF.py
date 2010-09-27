# Anisotropic network force field
#
# Written by Konrad Hinsen
#

"""
Anisotropic Network Model
"""

__docformat__ = 'epytext'

from MMTK.ForceFields.ForceField import ForceField, ForceFieldData
from MMTK import Utility
from MMTK_forcefield import NonbondedList, NonbondedListTerm
from MMTK_deformation import ANTerm
from Scientific import N

#
# The deformation force field class
#
class AnisotropicNetworkForceField(ForceField):

    """
    Effective harmonic force field for an ANM protein model
    """

    def __init__(self, cutoff = None, scale_factor = 1.):
        """
        @param cutoff: the cutoff for pair interactions. Pair interactions
                       in periodic systems are calculated using
                       the minimum-image convention; the cutoff should
                       therefore never be larger than half the smallest
                       edge length of the elementary cell.
        @type cutoff: C{float}
        @param scale_factor: a global scaling factor
        @type scale_factor: C{float}
        """
        ForceField.__init__(self, 'anisotropic_network')
        self.arguments = (cutoff,)
        self.cutoff = cutoff
        self.scale_factor = scale_factor

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
        ev = ANTerm(universe._spec, nbl, cutoff, self.scale_factor)
        return [update, ev]
