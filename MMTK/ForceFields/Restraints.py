# Harmonic restraint terms that can be added to a force field.
#
# Written by Konrad Hinsen
# last revision: 2005-8-30
#

"""This module contains harmonic restraint terms that can be added
to any force field.

Example:

from MMTK import *
from MMTK.ForceFields import Amber94ForceField
from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint

universe = InfiniteUniverse()
universe.protein = Protein('bala1')
force_field = Amber94ForceField() + \
              HarmonicDistanceRestraint(universe.protein[0][1].peptide.N,
                                        universe.protein[0][1].peptide.O,
                                        0.5, 10.)
universe.setForceField(force_field)
"""

from ForceField import ForceField
from MMTK_forcefield import HarmonicDistanceTerm, HarmonicAngleTerm, \
                            CosineDihedralTerm
import Numeric

class HarmonicDistanceRestraint(ForceField):

    """Harmonic distance restraint between two atoms

    Constructor: HarmonicDistanceRestraint(|atom1|, |atom2|,
                                           |distance|, |force_constant|)

    Arguments:

    |atom1|, |atom2| -- the two atoms whose distance is restrained

    |distance| -- the distance at which the restraint is zero

    |force_constant| -- the force constant of the restraint term

    The functional form of the restraint is
    |force_constant|*((r1-r2).length()-|distance|)**2, where
    r1 and r2 are the positions of the two atoms.
    """

    def __init__(self, atom1, atom2, distance, force_constant):
        if type(atom1) is type(0) and type(atom2) is type(0):
            self.index1 = atom1
            self.index2 = atom2
            self.atom1 = None
            self.atom2 = None
            self.universe = None
        else:
            self.atom1 = atom1
            self.atom2 = atom2
            self.universe = atom1.universe()
            if self.universe is not atom2.universe():
                raise ValueError("Atoms " + `self.atom1` + 'and' +
                                  `self.atom2` + 'are not in the same universe')
            self.universe.configuration()
            self.index1 = self.atom1.index
            self.index2 = self.atom2.index
        self.arguments = (self.index1, self.index2, distance, force_constant) 
        self.distance = distance
        self.force_constant = force_constant
        ForceField.__init__(self, 'harmonic distance restraint')

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if self.universe is not None and self.universe is not universe:
            raise ValueError("Atoms are not in universe " + `universe`)
        if subset1 is not None:
            s1 = subset1.atomList()
            s2 = subset2.atomList()
            if not ((self.atom1 in s1 and self.atom2 in s2) or \
                    (self.atom1 in s2 and self.atom2 in s1)):
                raise ValueError("restraint outside subset")
        indices = Numeric.array([[self.index1, self.index2]])
        parameters = Numeric.array([[self.distance, self.force_constant]])
        return [HarmonicDistanceTerm(universe._spec, indices, parameters,
                                     self.name)]

    def description(self):
        return 'ForceFields.Restraints.' + self.__class__.__name__ + \
               `self.arguments`


class HarmonicAngleRestraint(ForceField):

    """Harmonic angle restraint between three atoms

    Constructor: HarmonicAngleRestraint(|atom1|, |atom2|, |atom3|,
                                        |angle|, |force_constant|)

    Arguments:

    |atom1|, |atom2|, |atom3| -- the three atoms whose angle is restrained;
    |atom2| is the central atom

    |angle| -- the angle at which the restraint is zero

    |force_constant| -- the force constant of the restraint term

    The functional form of the restraint is
    |force_constant|*(phi-|angle|)**2, where
    phi is the angle |atom1|-|atom2|-|atom3|.
    """

    def __init__(self, atom1, atom2, atom3, angle, force_constant):
        if type(atom1) is type(0) and type(atom2) is type(0) and \
               type(atom3) is type(0):
            self.index1 = atom1
            self.index2 = atom2
            self.index3 = atom3
            self.atom1 = None
            self.atom2 = None
            self.atom3 = None
            self.universe = None
        else:
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.universe = atom1.universe()
            if self.universe is not atom2.universe() or \
                   self.universe is not atom3.universe():
                raise ValueError("Atoms " + `self.atom1` + ', ' +
                                  `self.atom2` + 'and' +
                                  `self.atom3` + 'are not in the same universe')
            self.universe.configuration()
            self.index1 = self.atom1.index
            self.index2 = self.atom2.index
            self.index3 = self.atom3.index
        self.arguments = (self.index1, self.index2, self.index3,
                          angle, force_constant) 
        self.angle = angle
        self.force_constant = force_constant
        ForceField.__init__(self, 'harmonic angle restraint')

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if self.universe is not None and self.universe is not universe:
            raise ValueError("Atoms are not in universe " + `universe`)
        indices = Numeric.array([[self.index1, self.index2, self.index3]])
        parameters = Numeric.array([[self.angle, self.force_constant]])
        return [HarmonicAngleTerm(universe._spec, indices, parameters,
                                  self.name)]

class HarmonicDihedralRestraint(ForceField):

    """Harmonic dihedral angle restraint between three atoms

    Constructor: HarmonicDihedralRestraint(|atom1|, |atom2|, |atom3|, |atom4|,
                                           |angle|, |force_constant|)

    Arguments:

    |atom1|, |atom2|, |atom3|, |atom4| -- the four atoms whose dihedral angle
    is restrained; |atom2| and |atom3| are on the common axis

    |angle| -- the dihedral angle at which the restraint is zero

    |force_constant| -- the force constant of the restraint term

    The functional form of the restraint is
    |force_constant|*(phi-|distance|)**2, where
    phi is the dihedral angle |atom1|-|atom2|-|atom3|-|atom4|.
    """

    def __init__(self, atom1, atom2, atom3, atom4, dihedral, force_constant):
        if type(atom1) is type(0) and type(atom2) is type(0) and \
               type(atom3) is type(0) and type(atom4) is type(0):
            self.index1 = atom1
            self.index2 = atom2
            self.index3 = atom3
            self.index4 = atom4
            self.atom1 = None
            self.atom2 = None
            self.atom3 = None
            self.atom4 = None
            self.universe = None
        else:
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.atom4 = atom4
            self.universe = atom1.universe()
            if self.universe is not atom2.universe() or \
                   self.universe is not atom3.universe() or \
                   self.universe is not atom3.universe():
                raise ValueError("Atoms " + `self.atom1` + ', ' +
                                  `self.atom2` + ', ' + `self.atom3` + 'and' +
                                  `self.atom4` + 'are not in the same universe')
            self.universe.configuration()
            self.index1 = self.atom1.index
            self.index2 = self.atom2.index
            self.index3 = self.atom3.index
            self.index4 = self.atom4.index
        self.dihedral = dihedral
        self.force_constant = force_constant
        self.arguments = (self.index1, self.index2, self.index3, self.index4,
                          dihedral, force_constant) 
        ForceField.__init__(self, 'harmonic dihedral restraint')

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if self.universe is not None and self.universe is not universe:
            raise ValueError("Atoms are not in universe " + `universe`)
        indices = Numeric.array([[self.index1, self.index2,
                                  self.index3, self.index4]])
        parameters = Numeric.array([[0., self.dihedral,
                                     0., self.force_constant]])
        return [CosineDihedralTerm(universe._spec, indices, parameters,
                                   self.name)]
