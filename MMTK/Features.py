# This module contains feature management classes.
#
# Written by Konrad Hinsen
# last revision: 2009-2-5
#

from MMTK import Environment
import string

# Each feature class represents a feature that a certain universe
# might or might not have: fixed atoms, constraints, thermostats, etc.
# Each non-trivial algorithm that works on a universe keeps a list of
# features it can handle. This arrangement ensures a minimal
# compatibility test between universes and algorithms.

# The feature list stores all defined features.

_all = []

# The feature base class ensures that each feature is a singleton.

class Feature(object):

    def __init__(self):
        for f in _all:
            if f.__class__ == self.__class__:
                raise ValueError("feature alredy defined")
        _all.append(self)

    # Method to be redefined by subclasses
    def isInUniverse(self, universe):
        raise TypeError("must be defined in subclass")

#
# Fixed particle feature
#
class FixedParticleFeatureClass(Feature):

    def isInUniverse(self, universe):
        fixed = universe.getAtomBooleanArray('fixed')
        return fixed.sumOverParticles() > 0

    description = 'fixed particles'

FixedParticleFeature = FixedParticleFeatureClass()

#
# Distance constraints feature
#
class DistanceConstraintsFeatureClass(Feature):

    def isInUniverse(self, universe):
        return universe.numberOfDistanceConstraints() > 0

    description = 'distance constraints'

DistanceConstraintsFeature = DistanceConstraintsFeatureClass()

#
# Path Integrals feature
#
class PathIntegralsFeatureClass(Feature):

    def isInUniverse(self, universe):
        for o in universe._environment:
            if o.__class__ is Environment.PathIntegrals:
                return True
        return False

    description = 'Path Integrals'

PathIntegralsFeature = PathIntegralsFeatureClass()

#
# Nose thermostat feature
#
class NoseThermostatFeatureClass(Feature):

    def isInUniverse(self, universe):
        for o in universe._environment:
            if o.__class__ is Environment.NoseThermostat:
                return True
        return False

    description = 'Nose thermostat'

NoseThermostatFeature = NoseThermostatFeatureClass()

#
# Andersen barostat feature
#
class AndersenBarostatFeatureClass(Feature):

    def isInUniverse(self, universe):
        for o in universe._environment:
            if o.__class__ is Environment.AndersenBarostat:
                return True
        return False

    description = 'Andersen barostat'

AndersenBarostatFeature = AndersenBarostatFeatureClass()

#
# Return feature list for a universe.
#
def getFeatureList(universe):
    return [f for f in _all if f.isInUniverse(universe)]

#
# Check that a feature list contains everything necessary for a universe.
#
def checkFeatures(algorithm, universe):
    universe_features = set(getFeatureList(universe))
    unsupported = universe_features.difference(set(algorithm.features))
    if unsupported:
        f = '\n'.join([f.description for f in unsupported])
        raise ValueError(algorithm.__class__.__name__ +
                          " does not support the following features:\n" + f)
    return universe_features
