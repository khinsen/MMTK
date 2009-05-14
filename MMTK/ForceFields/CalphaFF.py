# C-alpha force field
#
# Written by Konrad Hinsen
# last revision: 2009-5-13
#

"""
C-alpha force field for protein normal mode analysis
(Elastic Network Model)
"""

__docformat__ = 'epytext'

from MMTK.ForceFields.ForceField import ForceField, ForceFieldData
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
        @param cutoff: the cutoff for pair interactions, should be
                       at least 2.5 nm. Pair interactions in periodic
                       systems are calculated using the minimum-image
                       convention; the cutoff should therefore never be
                       larger than half the smallest edge length of the
                       elementary cell.
        @type cutoff: C{float}
        @param scale_factor: a global scaling factor
        @type scale_factor: C{float}
        """
        ForceField.__init__(self, 'calpha')
        self.arguments = (cutoff,)
        self.cutoff = cutoff
        self.scale_factor = scale_factor
        self.version = version

    def ready(self, global_data):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None:
            for s1, s2 in [(subset1, subset2), (subset2, subset1)]:
                set = {}
                for a in s1.atomList():
                    set[a.index] = None
                for a in s2.atomList():
                    try:
                        del set[a.index]
                    except KeyError: pass
            set = {}
            for a in subset1.atomList():
                set[a.index] = None
            for a in subset2.atomList():
                set[a.index] = None
            atom_subset = set.keys()
            atom_subset.sort()
            atom_subset = N.array(atom_subset)
        else:
            atom_subset = N.array([], N.Int)
        nothing = N.zeros((0, 2), N.Int)
        nbl = NonbondedList(nothing, nothing, atom_subset, universe._spec,
                            self.cutoff)
        update = NonbondedListTerm(nbl)
        cutoff = self.cutoff
        if cutoff is None:
            cutoff = 0.
        ev = CalphaTerm(universe._spec, nbl, cutoff,
                        self.scale_factor, self.version)
        return [update, ev]
