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
from MMTK import Utility
from MMTK_forcefield import HarmonicDistanceTerm, HarmonicAngleTerm, \
                            CosineDihedralTerm
from MMTK_restraints import HarmonicCMTrapTerm, HarmonicCMDistanceTerm
from Scientific import N

class HarmonicDistanceRestraint(ForceField):

    """
    Harmonic distance restraint between two atoms
    """

    def __init__(self, obj1, obj2, distance, force_constant,
                 nb_exclusion=False):
        """
        :param obj1: the object defining center-of-mass 1
        :type obj1: :class:`~MMTK.Collections.GroupOfAtoms`
        :param obj2: the object defining center-of-mass 2
        :type obj2: :class:`~MMTK.Collections.GroupOfAtoms`
        :param distance: the distance between cm 1 and cm2 at which
                         the restraint is zero
        :type distance: float
        :param force_constant: the force constant of the restraint term.
                               The functional form of the restraint is
                               force_constant*((cm1-cm2).length()-distance)**2,
                               where cm1 and cm2 are the centrer-of-mass
                               positions of the two objects.
        :type force_constant: float
        :param nb_exclussion: if True, non-bonded interactions between
                              the restrained atoms are suppressed, as
                              for a chemical bond
        :type nb_exclussion: bool
        """
        if isinstance(obj1, int) and isinstance(obj2, int):
            # Older MMTK versions admitted only single atoms and
            # stored single indices. Support this mode for opening
            # trajectories made with those versions
            self.atom_indices_1 = [obj1]
            self.atom_indices_2 = [obj2]
        self.atom_indices_1 = self.getAtomParameterIndices(obj1.atomList())
        self.atom_indices_2 = self.getAtomParameterIndices(obj2.atomList())
        if nb_exclusion and (len(self.atom_indices_1) > 1
                             or len(self.atom_indices_2) > 1):
            raise ValueError("Non-bonded exclusion possible only "
                             "between single-atom objects")
        self.arguments = (self.atom_indices_1, self.atom_indices_2,
                          distance, force_constant, nb_exclusion)
        self.distance = distance
        self.force_constant = force_constant
        self.nb_exclusion = nb_exclusion
        ForceField.__init__(self, 'harmonic distance restraint')

    def supportsPathIntegrals(self):
        return True

    def declareDependencies(self, global_data):
        global_data.add('nb_exclusions', self.__class__)

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        if universe.is_periodic and \
                (len(self.atom_indices_1) > 1 or len(self.atom_indices_1) > 1):
            raise ValueError("Center-of-mass restraints not implemented"
                             " for periodic universes")
        ok = False
        for s1, s2 in [(subset1, subset2), (subset2, subset1)]:
            if s1 is None and s1 is None:
                ok = True
                break
            s1 = set(a.index for a in s1.atomIterator())
            diff1 = set(self.atom_indices_1).difference(s1)
            s2 = set(a.index for a in s2.atomIterator())
            diff2 = set(self.atom_indices_2).difference(s2)
            if not diff1 and not diff2:
                # Each object is in one of the subsets
                ok = True
                break
            if (diff1 and len(diff1) != len(self.atom_indices_1)) \
                    or (diff2 and len(diff2) != len(self.atom_indices_2)):
                # The subset contains some but not all of the
                # restrained atoms.
                raise ValueError("Restrained atoms partially "
                                 "in a subset")
        global_data.add('initialized', self.__class__)
        if not ok:
            # The objects are not in the subsets, so there is no
            # contribution to the total energy.
            return {'harmonic_distance_cm': []}
        if self.nb_exclusion:
            assert len(self.atom_indices_1) == 1 
            assert len(self.atom_indices_2) == 1
            i1 = self.atom_indices_1[0]
            i2 = self.atom_indices_2[0]
            global_data.add('excluded_pairs', Utility.normalizePair((i1, i2)))
        if len(self.atom_indices_1) == 1 and len(self.atom_indices_2) == 1:
            # Keep the old format for the single-atom case for best
            # compatibility with older MMTK versions.
            f, offsets = self.beadOffsetsAndFactor([self.atom_indices_1[0],
                                                    self.atom_indices_2[0]],
                                                   global_data)
            return {'harmonic_distance_term':
                    [(self.atom_indices_1[0]+o1, self.atom_indices_2[0]+o2,
                      self.distance, f*self.force_constant)
                     for o1, o2 in offsets]}
        else:
            f, offsets = self.beadOffsetsAndFactor(list(self.atom_indices_1)+
                                                   list(self.atom_indices_2),
                                                   global_data)
            n1 = len(self.atom_indices_1)
            offsets = [(o[:n1], o[n1:]) for o in offsets]
            return {'harmonic_distance_cm':
                    [(self.atom_indices_1+o1, self.atom_indices_2+o2,
                      self.distance, f*self.force_constant)
                     for o1, o2 in offsets]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)
        if params.has_key('harmonic_distance_term'):
            params = params['harmonic_distance_term']
            return [HarmonicDistanceTerm(universe._spec,
                                         N.array([p[:2]]),
                                         N.array([p[2:]]),
                                         self.name)
                    for p in params]
        else:
            params = params['harmonic_distance_cm']
            return [HarmonicCMDistanceTerm(universe,
                                           N.array(p[0]),
                                           N.array(p[1]),
                                           universe.masses().array,
                                           p[2],
                                           p[3])
                    for p in params]

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
        if universe.is_periodic and len(self.atom_indices) > 1:
            raise ValueError("Center-of-mass restraints not implemented"
                             " for periodic universes")
        for subset in [subset1, subset2]:
            if subset is not None:
                subset = set(a.index for a in subset.atomIterator())
                diff = set(self.atom_indices).difference(subset)
                if diff:
                    if len(diff) == len(self.atom_indices):
                        # The subset doesn't contain the restrained atoms.
                        return {'harmonic_trap_cm': []}
                    else:
                        # The subset contains some but not all of the
                        # restrained atoms.
                        raise ValueError("Restrained atoms partially "
                                         "in a subset")
        return {'harmonic_trap_cm':
                [(self.atom_indices, self.center, self.force_constant)]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['harmonic_trap_cm']
        return [HarmonicCMTrapTerm(universe,
                                   N.array(p[0]),
                                   universe.masses().array,
                                   p[1],
                                   p[2])
                for p in params]
