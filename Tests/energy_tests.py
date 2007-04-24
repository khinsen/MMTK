# Energy tests
#
# Written by Konrad Hinsen
# last revision: 2007-4-24
#

import unittest
from MMTK import *
from MMTK_forcefield import NonbondedList
from MMTK.Random import randomPointInBox
from Scientific import N


def sorted_tuple(pair):
    if pair[0] < pair[1]:
        return (pair[0], pair[1])
    else:
        return (pair[1], pair[0])


class NonbondedListTest:

    def test_nonbondedList(self):

        self.universe.configuration()
        atoms = self.universe.atomList()
        atom_indices = N.array([a.index for a in self.universe.atomList()])
        empty = N.zeros((0, 2), N.Int)

        for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]:

            nblist = NonbondedList(empty, empty, atom_indices,
                                   self.universe._spec, cutoff)
            nblist.update(self.universe.configuration().array)
            distances = nblist.pairDistances()
            pairs1 = nblist.pairIndices()
            pairs1 = [sorted_tuple(pairs1[i]) for i in range(len(pairs1))
                      if distances[i] < cutoff]
            pairs1.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            pairs2 = []
            for i in range(len(atoms)):
                for j in range(i+1, len(atoms)):
                    d = self.universe.distance(atoms[i], atoms[j])
                    if d < cutoff:
                        pairs2.append(sorted_tuple((atoms[i].index,
                                                    atoms[j].index)))
            pairs2.sort(lambda a, b: cmp(a[0], b[0]) or cmp(a[1], b[1]))

            self.assertEqual(pairs1, pairs2)


class InfiniteUniverseNonbondedListTest(unittest.TestCase,
                                        NonbondedListTest):
    def setUp(self):
        self.universe = InfiniteUniverse()
        for i in range(100):
            p = randomPointInBox(1.)
            self.universe.addObject(Atom('C', position = p))

class OrthorhombicUniverseNonbondedListTest(unittest.TestCase,
                                            NonbondedListTest):
    def setUp(self):
        self.universe = OrthorhombicPeriodicUniverse((1.5, 0.8, 1.1))
        for i in range(100):
            p = self.universe.boxToRealCoordinates(randomPointInBox(1.))
            self.universe.addObject(Atom('C', position = p))

class ParallelepipedicUniverseNonbondedListTest(unittest.TestCase,
                                                NonbondedListTest):
    def setUp(self):
        a = Vector(1.5, 0.3, 0.)
        b = Vector(-0.1, 0.8, 0.2)
        c = Vector(0., -0.2, 1.1)
        self.universe = ParallelepipedicPeriodicUniverse((a, b, c))
        for i in range(100):
            p = self.universe.boxToRealCoordinates(randomPointInBox(1.))
            self.universe.addObject(Atom('C', position = p))


if __name__ == '__main__':
    unittest.main()
