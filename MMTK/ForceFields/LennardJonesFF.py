# A Lennard-Jones force fields for simple liquids
#
# Written by Konrad Hinsen
#

"""
Lennard-Jones force field for simple liquids
"""

__docformat__ = 'epytext'

from MMTK.ForceFields.NonBondedInteractions import LJForceField
from Scientific.Geometry import Vector
import copy

#
# The force field class
#
class LennardJonesForceField(LJForceField):

    """
    Lennard-Jones force field

    The Lennard-Jones parameters are taken from the atom attributes
    LJ_radius and LJ_energy. The pair interaction energy has the form
    U(r)=4*LJ_energy*((LJ_radius/r)**12-(LJ_radius/r)**6).
    """

    def __init__(self, cutoff = None):
        """
        @param cutoff: a cutoff value or C{None}, meaning no cutoff.
                       Pair interactions in periodic systems are calculated
                       using the minimum-image convention; the cutoff should
                       therefore never be larger than half the smallest edge
                       length of the elementary cell.
        @type cutoff: C{float}
        """
        self.arguments = (cutoff, )
        LJForceField.__init__(self, 'LJ', cutoff)
        self.lj_14_factor = 1.

    def ready(self, global_data):
        return True

    def collectParameters(self, universe, global_data):
        if not hasattr(global_data, 'lj_parameters'):
            parameters = {}
            for o in universe:
                for a in o.atomList():
                    parameters[a.symbol] = (a.LJ_energy, a.LJ_radius, 0)
            global_data.lj_parameters = parameters

    def _atomType(self, object, atom, global_data):
        return atom.symbol

    def _ljParameters(self, type, global_data):
        return global_data.lj_parameters[type]

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        self.collectParameters(universe, global_data)
        return LJForceField.evaluatorParameters(self, universe, subset1,
                                                subset2, global_data)
