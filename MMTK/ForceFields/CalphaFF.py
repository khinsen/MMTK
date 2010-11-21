# C-alpha force field
#
# Written by Konrad Hinsen
#

"""
C-alpha force field for protein normal mode analysis
(Elastic Network Model)
"""

__docformat__ = 'restructuredtext'

from MMTK.ForceFields.ForceField import ForceField, ForceFieldData
from MMTK import Utility
from MMTK_forcefield import NonbondedList, NonbondedListTerm
from MMTK_deformation import CalphaTerm
from Scientific import N

#
# The deformation force field class
#
class CalphaForceField(ForceField):

    """
    Effective harmonic force field for a C-alpha protein model
    """

    def __init__(self, cutoff = None, scale_factor = 1., version=1):
        """
        :param cutoff: the cutoff for pair interactions, should be
                       at least 2.5 nm. Pair interactions in periodic
                       systems are calculated using the minimum-image
                       convention; the cutoff should therefore never be
                       larger than half the smallest edge length of the
                       elementary cell.
        :type cutoff: float
        :param scale_factor: a global scaling factor
        :type scale_factor: float
        """
        ForceField.__init__(self, 'calpha')
        self.arguments = (cutoff,)
        self.cutoff = cutoff
        self.scale_factor = scale_factor
        self.version = version

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
        nbl = NonbondedList(excluded_pairs, nothing, atom_subset, universe._spec,
                            self.cutoff)
        update = NonbondedListTerm(nbl)
        cutoff = self.cutoff
        if cutoff is None:
            cutoff = 0.
        ev = CalphaTerm(universe._spec, nbl, cutoff,
                        self.scale_factor, self.version)
        return [update, ev]
