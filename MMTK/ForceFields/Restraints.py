# Harmonic restraint terms that can be added to a force field.
#
# Written by Konrad Hinsen
#

"""
Harmonic restraint terms that can be added to any force field

Example::

 from MMTK import *
 from MMTK.ForceFields import Amber99ForceField
 from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint
 
 universe = InfiniteUniverse()
 universe.protein = Protein('bala1')
 force_field = Amber99ForceField() + \
               HarmonicDistanceRestraint(universe.protein[0][1].peptide.N,
                                         universe.protein[0][1].peptide.O,
                                         0.5, 10.)
 universe.setForceField(force_field)
"""

__docformat__ = 'restructuredtext'

from MMTK.ForceFields.ForceField import ForceField
from MMTK_forcefield import HarmonicDistanceTerm, HarmonicAngleTerm, \
                            CosineDihedralTerm
from MMTK_restraints import HarmonicTrapTerm
from Scientific import N

class HarmonicDistanceRestraint(ForceField):

    """
    Harmonic distance restraint between two atoms
    """

    def __init__(self, atom1, atom2, distance, force_constant):
        """
        :param atom1: first atom
        :type atom1: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom2: second atom
        :type atom2: :class:`~MMTK.ChemicalObjects.Atom`
        :param distance: the distance at which the restraint is zero
        :type distance: float
        :param force_constant: the force constant of the restraint term.
                               The functional form of the restraint is
                               force_constant*((r1-r2).length()-distance)**2,
                               where r1 and r2 are the positions of the
                               two atoms.
        """
        self.index1, self.index2 = self.getAtomParameterIndices((atom1, atom2))
        self.arguments = (self.index1, self.index2, distance, force_constant) 
        self.distance = distance
        self.force_constant = force_constant
        ForceField.__init__(self, 'harmonic distance restraint')

    def supportsPathIntegrals(self):
        return True

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        if subset1 is not None:
            s1 = subset1.atomList()
            s2 = subset2.atomList()
            if not ((self.atom1 in s1 and self.atom2 in s2) or \
                    (self.atom1 in s2 and self.atom2 in s1)):
                raise ValueError("restraint outside subset")
        f, offsets = self.beadOffsetsAndFactor([self.index1, self.index2],
                                               global_data)
        return {'harmonic_distance_term':
                [(self.index1+o1, self.index2+o2,
                  self.distance, f*self.force_constant)
                 for o1, o2 in offsets]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['harmonic_distance_term']
        indices = N.array([params[0][:2]])
        parameters = N.array([params[0][2:]])
        return [HarmonicDistanceTerm(universe._spec, indices, parameters,
                                     self.name)]

    def description(self):
        return 'ForceFields.Restraints.' + self.__class__.__name__ + \
               `self.arguments`


class HarmonicAngleRestraint(ForceField):

    """
    Harmonic angle restraint between three atoms
    """

    def __init__(self, atom1, atom2, atom3, angle, force_constant):
        """
        :param atom1: first atom
        :type atom1: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom2: second (central) atom
        :type atom2: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom3: third atom
        :type atom3: :class:`~MMTK.ChemicalObjects.Atom`
        :param angle: the angle at which the restraint is zero
        :type angle: float
        :param force_constant: the force constant of the restraint term.
                               The functional form of the restraint is
                               force_constant*(phi-angle)**2, where
                               phi is the angle atom1-atom2-atom3.
        """
        self.index1, self.index2, self.index3 = \
                    self.getAtomParameterIndices((atom1, atom2, atom3))
        self.arguments = (self.index1, self.index2, self.index3,
                          angle, force_constant) 
        self.angle = angle
        self.force_constant = force_constant
        ForceField.__init__(self, 'harmonic angle restraint')

    def supportsPathIntegrals(self):
        return True

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        f, offsets = self.beadOffsetsAndFactor([self.index1, self.index2, self.index3],
                                               global_data)
        return {'harmonic_angle_term':
                [(self.index1+o1, self.index2+o2, self.index3+o3,
                   self.angle, f*self.force_constant)
                 for o1, o2, o3 in offsets]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['harmonic_angle_term']
        indices = N.array([params[0][:3]])
        parameters = N.array([params[0][3:]])
        return [HarmonicAngleTerm(universe._spec, indices, parameters,
                                  self.name)]

class HarmonicDihedralRestraint(ForceField):

    """
    Harmonic dihedral angle restraint between four atoms
    """

    def __init__(self, atom1, atom2, atom3, atom4, dihedral, force_constant):
        """
        :param atom1: first atom
        :type atom1: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom2: second (axis) atom
        :type atom2: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom3: third (axis)atom
        :type atom3: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom4: fourth atom
        :type atom4: :class:`~MMTK.ChemicalObjects.Atom`
        :param dihedral: the dihedral angle at which the restraint is zero
        :type dihedral: float
        :param force_constant: the force constant of the restraint term.
                               The functional form of the restraint is
                               force_constant*(phi-abs(dihedral))**2, where
                               phi is the dihedral angle
                               atom1-atom2-atom3-atom4.
        """
        self.index1, self.index2, self.index3, self.index4 = \
                   self.getAtomParameterIndices((atom1, atom2, atom3, atom4))
        self.dihedral = dihedral
        self.force_constant = force_constant
        self.arguments = (self.index1, self.index2, self.index3, self.index4,
                          dihedral, force_constant) 
        ForceField.__init__(self, 'harmonic dihedral restraint')

    def supportsPathIntegrals(self):
        return True

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        f, offsets = self.beadOffsetsAndFactor([self.index1, self.index2,
                                                self.index3, self.index4],
                                               global_data)
        return {'cosine_dihedral_term':
                [(self.index1+o1, self.index2+o2,
                  self.index3+o3, self.index4+o4,
                  0., self.dihedral,
                  0., f*self.force_constant)
                 for o1, o2, o3, o4 in offsets]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['cosine_dihedral_term']
        indices = N.array([params[0][:4]])
        parameters = N.array([params[0][4:]])
        return [CosineDihedralTerm(universe._spec, indices, parameters,
                                   self.name)]

class HarmonicTrapForceField(ForceField):

    """
    Harmonic potential with respect to a fixed point in space
    """

    def __init__(self, obj, center, force_constant):
        """
        :param obj: the object on whose center of mass the force field acts
        :type obj: :class:`~MMTK.Collections.GroupOfAtoms`
        :param center: the point to which the atom is attached by
                        the harmonic potential
        :type center: Scientific.Geometry.Vector
        :param force_constant: the force constant of the harmonic potential
        :type force_constant: float
        """
        self.atom_indices = self.getAtomParameterIndices(obj.atomList())
        self.arguments = (self.atom_indices, center, force_constant)
        ForceField.__init__(self, 'harmonic_trap')
        self.center = center
        self.force_constant = force_constant

    def ready(self, global_data):
        return True

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        raise NotImplementedError

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        for subset in [subset1, subset2]:
            if subset is not None:
                subset = set(a.index for a in subset.atomIterator())
                diff = set(self.atom_indices).difference(subset)
                if diff:
                    if len(diff) == len(self.atom_indices):
                        # The subset doesn't contain the restrained atoms.
                        return []
                    else:
                        # The subset contains some but not all of the
                        # restrained atoms.
                        raise ValueError("Restrained atoms partially "
                                         "in a subset")
        return [HarmonicTrapTerm(universe,
                                 N.array(self.atom_indices),
                                 universe.masses().array,
                                 self.center,
                                 self.force_constant)]
