# This module defines environment objects for universes.
#
# Written by Konrad Hinsen
#

"""
Environment objects

Environment objects are objects that define a simulation system without
being composed of atoms. Examples are thermostats, barostats, external
fields, etc.
"""

__docformat__ = 'restructuredtext'

from Scientific import N

#
# The environment object base class
#
class EnvironmentObject(object):

    is_environment_object = 1

    def checkCompatibilityWith(self, other):
        pass

    def description(self):
        return "o('Environment." + self.__class__.__name__ + \
               `tuple(self.parameters)` + "')"

# Type check

def isEnvironmentObject(object):
    return hasattr(object, 'is_environment_object')

#
# Nose thermostat class
#
class NoseThermostat(EnvironmentObject):

    """
    Nose thermostat for Molecular Dynamics

    A thermostat object can be added to a universe and will then
    modify the integration algorithm to a simulation of an NVT
    ensemble.
    """

    def __init__(self, temperature, relaxation_time = 0.2):
        """
        :param temperature: the temperature set by the thermostat
        :type temperature: float
        :param relaxation_time: the relaxation time of the
                                thermostat coordinate
        :type relaxation_time: float
        """
        self.arguments = (temperature, relaxation_time)
        self.parameters = N.array([temperature, relaxation_time])
        self.coordinates = N.array([0., 0.])

    def setTemperature(self, temperature):
        self.parameters[0] = temperature

    def setRelaxationTime(self, t):
        self.parameters[1] = t

    def checkCompatibilityWith(self, other):
        if other.__class__ is NoseThermostat:
            raise ValueError("the universe already has a thermostat")

#
# Andersen barostat class
#
class AndersenBarostat(EnvironmentObject):

    """
    Andersen barostat for Molecular Dynamics

    A barostat object can be added to a universe and will then
    together with a thermostat object modify the integration algorithm
    to a simulation of an NPT ensemble.
    """

    def __init__(self, pressure, relaxation_time = 1.5):
        """
        :param pressure: the pressure set by the barostat
        :type pressure: float
        :param relaxation_time: the relaxation time of the
                                barostat coordinate
        :type relaxation_time: float
        """
        self.arguments = (pressure, relaxation_time)
        self.parameters = N.array([pressure, relaxation_time])
        self.coordinates = N.array([0.])

    def setPressure(self, pressure):
        self.parameters[0] = pressure

    def setRelaxationTime(self, t):
        self.parameters[1] = t

    def checkCompatibilityWith(self, other):
        if other.__class__ is AndersenBarostat:
            raise ValueError("the universe already has a barostat")

