# This module contains code for charge fitting.
#
# Written by Konrad Hinsen
#

"""
Fit of point chages to an electrostatic potential surface

This module implements a numerically stable method (based on
Singular Value Decomposition) to fit point charges to values of an
electrostatic potential surface. Two types of constraints are
avaiable: a constraint on the total charge of the system or a subset
of the system, and constraints that force the charges of several atoms
to be equal. There is also a utility function that selects suitable
evaluation points for the electrostatic potential surface. For the
potential evaluation itself, some quantum chemistry program is needed.

The charge fitting method is described in:

  - K. Hinsen and B. Roux,
    An accurate potential for simulating proton transfer in acetylacetone,
    J. Comp. Chem. 18, 1997: 368

See also Examples/Miscellaneous/charge_fit.py.
"""

__docformat__ = 'restructuredtext'

from MMTK import Random, Units, Utility
from Scientific.Geometry import Vector
from Scientific import N
from Scientific import LA


class ChargeFit(object):

    """
    Fit of point charges to an electrostatic potential surface

    A ChargeFit object acts like a dictionary that stores the fitted charge
    value for each atom in the system.
    """

    def __init__(self, system, points, constraints = None):
        """
        :param system: any chemical object (usually a molecule)
        :param points: a list of point/potential pairs (a vector for the
                       evaluation point, a number for the potential),
                       or a dictionary whose keys are Configuration objects
                       and whose values are lists of point/potential pairs.
                       The latter case permits combined fits for several
                       conformations of the system.
        :param constraints: an optional list of constraint objects
                            (:class:`~MMTK.ChargeFit.TotalChargeConstraint`
                            and/or
                            :class:`~MMTK.ChargeFit.EqualityConstraint` objects). 
                            If the constraints are inconsistent, a warning is
                            printed and the result will satisfy the
                            constraints only in a least-squares sense.
        """
        self.atoms = system.atomList()
        if type(points) != type({}):
            points = {None: points}
        if constraints is not None:
            constraints = ChargeConstraintSet(self.atoms, constraints)
        npoints = sum([len(v) for v in points.values()])
        natoms = len(self.atoms)
        if npoints < natoms:
            raise ValueError("Not enough data points for fit")

        m = N.zeros((npoints, natoms), N.Float)
        phi = N.zeros((npoints,), N.Float)
        i = 0
        for conf, pointlist in points.items():
            for r, p in pointlist:
                for j in range(natoms):
                    m[i, j] = 1./(r-self.atoms[j].position(conf)).length()
                phi[i] = p
                i = i + 1
        m = m*Units.electrostatic_energy

        m_test = m
        phi_test = phi

        if constraints is not None:
            phi -= N.dot(m, constraints.bi_c)
            m = N.dot(m, constraints.p)
            c_rank = constraints.rank
        else:
            c_rank = 0

        u, s, vt = LA.singular_value_decomposition(m)
        s_test = s[:len(s)-c_rank]
        cutoff = 1.e-10*N.maximum.reduce(s_test)
        nonzero = N.repeat(s_test, N.not_equal(s_test, 0.))
        self.rank = len(nonzero)
        self.condition = N.maximum.reduce(nonzero) / \
                         N.minimum.reduce(nonzero)
        self.effective_rank = N.add.reduce(N.greater(s, cutoff))
        if self.effective_rank < self.rank:
            self.effective_condition = N.maximum.reduce(nonzero) / cutoff
        else:
            self.effective_condition = self.condition
        if self.effective_rank < natoms-c_rank:
            Utility.warning('Not all charges are uniquely determined' +
                            ' by the available data')

        for i in range(natoms):
            if s[i] > cutoff:
                s[i] = 1./s[i]
            else:
                s[i] = 0.
        q = N.dot(N.transpose(vt),
                  s*N.dot(N.transpose(u)[:natoms, :], phi))
        if constraints is not None:
            q = constraints.bi_c + N.dot(constraints.p, q)

        deviation = N.dot(m_test, q)-phi_test
        self.rms_error = N.sqrt(N.dot(deviation, deviation))
        deviation = N.fabs(deviation/phi_test)
        self.relative_rms_error = N.sqrt(N.dot(deviation, deviation))

        self.charges = {}
        for i in range(natoms):
            self.charges[self.atoms[i]] = q[i]

    def __getitem__(self, item):
        return self.charges[item]


class TotalChargeConstraint(object):

    """
    Constraint on the total system charge

    To be used with :class:`~MMTK.ChargeFit.ChargeFit`
    """

    def __init__(self, system, charge):
        """
        :param system: any chamical object whose total charge
                       is to be constrained
        :param charge: the total charge value
        :type charge: number
        """
        self.atoms = system.atomList()
        self.charge = charge

    def __len__(self):
        return 1

    def setCoefficients(self, atoms, b, c, i):
        for a in self.atoms:
            j = atoms.index(a)
            b[i, j] = 1.
        c[i] = self.charge


class EqualityConstraint(object):

    """
    Constraint forcing two charges to be equal

    To be used with :class:`~MMTK.ChargeFit.ChargeFit`

    Any atom may occur in more than one EqualityConstraint object,
    in order to keep the charges of more than two atoms equal.
    """

    def __init__(self, atom1, atom2):
        """
        :param atom1: the first atom in the equality relation
        :type atom1: :class:`~MMTK.ChemicalObjects.Atom`
        :param atom2: the second atom in the equality relation
        :type atom2: :class:`~MMTK.ChemicalObjects.Atom`
        """
        self.a1 = atom1
        self.a2 = atom2

    def __len__(self):
        return 1

    def setCoefficients(self, atoms, b, c, i):
        b[i, atoms.index(self.a1)] = 1.
        b[i, atoms.index(self.a2)] = -1.
        c[i] = 0.


class ChargeConstraintSet(object):

    def __init__(self, atoms, constraints):
        self.atoms = atoms
        natoms = len(self.atoms)
        nconst = sum([len(c) for c in constraints])
        b = N.zeros((nconst, natoms), N.Float)
        c = N.zeros((nconst,), N.Float)
        i = 0
        for cons in constraints:
            cons.setCoefficients(self.atoms, b, c, i)
            i = i + len(cons)
        u, s, vt = LA.singular_value_decomposition(b)
        self.rank = 0
        for i in range(min(natoms, nconst)):
            if s[i] > 0.:
                self.rank = self.rank + 1
        self.b = b
        self.bi = LA.generalized_inverse(b)
        self.p = N.identity(natoms)-N.dot(self.bi, self.b)
        self.c = c
        self.bi_c = N.dot(self.bi, c)
        c_test = N.dot(self.b, self.bi_c)
        if N.add.reduce((c_test-c)**2)/nconst > 1.e-12:
            Utility.warning("The charge constraints are inconsistent."
                            " They will be applied as a least-squares"
                            " condition.")


def evaluationPoints(system, n, smallest = 0.3, largest = 0.5):
    """
    Generate points in space around a molecule that are suitable
    for potential evaluation in view of a subsequent charge fit.
    The points are chosen at random and uniformly in a shell around the system.

    :param system: the chemical object for which the charges
                   will be fitted
    :param n: the number of evaluation points to be generated
    :param smallest: the smallest allowed distance of any evaluation
                     point from any non-hydrogen atom
    :param largest: the largest allowed value for the distance
                    from an evaluation point to the nearest atom
    :returns: a list of evaluation points
    :rtype: list of Scientific.Geometry.Vector
    """
    atoms = system.atomList()
    p1, p2 = system.boundingBox()
    margin = Vector(largest, largest, largest)
    p1 -= margin
    p2 += margin
    a, b, c = tuple(p2-p1)
    offset = 0.5*Vector(a, b, c)
    points = []
    while len(points) < n:
        p = p1 + Random.randomPointInBox(a, b, c) + offset
        m = 2*largest
        ok = 1
        for atom in atoms:
            d = (p-atom.position()).length()
            m = min(m, d)
            if d < smallest and atom.symbol != 'H':
                ok = 0
            if not ok: break
        if ok and m <= largest:
            points.append(p)
    return points


if __name__ == '__main__':

    from MMTK import *

    a1 = Atom('C', position=Vector(-0.05,0.,0.))
    a2 = Atom('C', position=Vector( 0.05,0.,0.))
    system = Collection(a1, a2)

    a1.charge = -0.75
    a2.charge = 0.15

    points = []
    for r in evaluationPoints(system, 50):
        p = 0.
        for atom in system.atomList():
            p = p + atom.charge/(r-atom.position()).length()
        points.append((r, p*Units.electrostatic_energy))

    constraints = [TotalChargeConstraint(system, 0.)]
    constraints = [EqualityConstraint(a1, a2)]
    constraints = None

    f = ChargeFit(system, points, constraints)
    print f[a1], a1.charge
    print f[a2], a2.charge
