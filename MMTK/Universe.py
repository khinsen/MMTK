# This module implements the various types of universes
# (infinite, periodic etc.). A universe defines the
# geometry of space, the force field, and external interactions
# (boundary conditions, external fields, etc.)
#
# Written by Konrad Hinsen
# last revision: 2010-2-11
#

"""
Universes
"""

__docformat__ = 'epytext'

from MMTK import Bonds, ChemicalObjects, Collections, Environment, \
                 Random, Utility, ParticleProperties, Visualization
from Scientific.Geometry import Transformation
from Scientific.Geometry import Vector, isVector
from Scientific import N
import copy

try:
    import threading
    if not hasattr(threading, 'Thread'):
        threading = None
except ImportError:
    threading = None

#
# The base class for all universes.
#
class Universe(Collections.GroupOfAtoms, Visualization.Viewable):

    """
    Universe

    A universe represents a complete model of a chemical system, i.e.
    the molecules, their environment (topology, boundary conditions,
    thermostats, etc.), and optionally a force field.

    The class Universe is an abstract base class that defines
    properties common to all kinds of universes. To create universe
    objects, use one of its subclasses.

    In addition to the methods listed below, universe objects support
    the following operations (C{u} is any universe object, C{o} is any
    chemical object):

     - C{len(u)} yields the number of chemical objects in the universe
     - C{u[i]} returns object number C{i}
     - C{u.name = o} adds C{o} to the universe and also makes it accessible as
       an attribute
     - C{del u.name} removes the object that was assigned to C{u.name} from
        the universe
    """

    def __init__(self, forcefield, properties):
        self._forcefield = forcefield
        self._evaluator = {}
        self.name = ''
        if properties.has_key('name'):
            self.name = properties['name']
            del properties['name']
        self._objects = Collections.Collection()
        self._environment = []
        self._configuration = None
        self._masses = None
        self._atom_properties = {}
        self._atoms = None
        self._bond_database = None
        self._bond_pairs = None
        self._version = 0
        self._np = None

    is_universe = True
    is_periodic = False
    is_orthogonal = False

    def __getstate__(self):
        state = copy.copy(self.__dict__)
        state['_evaluator'] = {}
        state['_configuration'] = None
        del state['_masses']
        del state['_bond_database']
        del state['_bond_pairs']
        del state['_np']
        del state['_spec']
        return state

    def __setstate__(self, state):
        state['_np'] = None
        state['_atoms'] = None
        state['_bond_database'] = None
        state['_bond_pairs'] = None
        self.__dict__['_environment'] = []
        if state.has_key('atom_properties'):
            self.__dict__['_atom_properties'] = state['atom_properties']
            del state['atom_properties']
        for attr, value in state.items():
            self.__dict__[attr] = value
        self._evaluator = {}
        self._masses = None
        self._createSpec()

    def __len__(self):
        return len(self._objects)

    def __getitem__(self, item):
        return self._objects[item]

    def __setattr__(self, attr, value):
        if attr[0] != '_' and self.__dict__.has_key(attr):
            try:
                self.removeObject(self.__dict__[attr])
            except ValueError:
                pass
        self.__dict__[attr] = value
        if attr[0] != '_' and (ChemicalObjects.isChemicalObject(value)
                               or Environment.isEnvironmentObject(value)):
            self.addObject(value)

    def __delattr__(self, attr):
        try:
            self.removeObject(self.__dict__[attr])
        except ValueError:
            pass
        del self.__dict__[attr]

    def __repr__(self):
        return self.__class__.__name__ + ' ' + self.name + ' containing ' + \
               `len(self._objects)` + ' objects.'
    __str__ = __repr__

    def __copy__(self):
        return copy.deepcopy(self)

    def objectList(self, klass = None):
        """
        @param klass: an optional class argument
        @type klass: C{class}
        @returns: a list of all chemical objects in the universe.
                  If klass is given, only objects that are instances
                  of klass are returned.
        @rtype: C{list}
        """
        return self._objects.objectList(klass)

    def environmentObjectList(self, klass = None):
        """
        @param klass: an optional class argument
        @type klass: C{class}
        @returns: a list of all environment objects in the universe.
                  If klass is given, only objects that are instances
                  of klass are returned.
        @rtype: C{list}
        """
        if klass is None:
            return self._environment
        else:
            return filter(lambda o, k=klass: o.__class__ is k,
                          self._environment)

    def atomList(self):
        """
        @returns: a list of all atoms in the universe
        @rtype: C{list}
        """
        if self._atoms is None:
            self._atoms = self._objects.atomList()
        return self._atoms

    def bondedUnits(self):
        return self._objects.bondedUnits()

    def universe(self):
        """
        @returns: the universe itself
        """
        return self

    def addObject(self, object, steal = False):
        """
        Adds object to the universe. If object is a Collection,
        all elements of the Collection are added to the universe.
        @param object: the object (chemical or environment) to be added
        @param steal: if C{True}, permit stealing the object from another
                      universe, otherwise the object must not yet be
                      attached to any universe.
        @type steal: C{bool}
        """
        if ChemicalObjects.isChemicalObject(object):
            if (not steal) and object.parent is not None:
                if isUniverse(object.parent):
                    raise ValueError(`object` +
                                      ' is already in another universe')
                else:
                    raise ValueError(`object` + ' is part of another object')
            object.parent = self
            self._objects.addObject(object)
            self._changed(True)
        elif Environment.isEnvironmentObject(object):
            for o in self._environment:
                o.checkCompatibilityWith(object)
            self._environment.append(object)
            self._changed(False)
        elif Collections.isCollection(object) \
             or Utility.isSequenceObject(object):
            for o in object:
                self.addObject(o, steal)
        else:
            raise TypeError(repr(object) + ' cannot be added to a universe')

    def removeObject(self, object):
        """
        Removes object from the universe. If object is a Collection,
        each of its elements is removed. The object to be removed must
        be in the universe.
        @param object: the object (chemical or environment) to be removed
        """
        if ChemicalObjects.isChemicalObject(object):
            if object.parent != self:
                raise ValueError(`object` + ' is not in this universe.')
            object.parent = None
            self._objects.removeObject(object)
            self._changed(True)
        elif Collections.isCollection(object) \
                 or Utility.isSequenceObject(object):
            for o in object:
                self.removeObject(o)
        elif Environment.isEnvironmentObject(object):
            self._environment.remove(object)
            self._changed(False)
        else:
            raise ValueError(`object` + ' is not in this universe.')

    def selectShell(self, point, r1, r2=0.):
        """
        @param point: a point in space
        @type point: C{Scientific.Geometry.Vector}
        @param r1: one of the radii of a spherical shell
        @type r1: C{float}
        @param r2: the other of the two radii of a spherical shell
        @type r2: C{float}
        @returns: a Collection of all objects in the universe whose
                  distance from point lies between r1 and r2.
        """
        return self._objects.selectShell(point, r1, r2)

    def selectBox(self, p1, p2):
        """
        @param p1: one corner of a box in space
        @type p1: C{Scientific.Geometry.Vector}
        @param p2: the other corner of a box in space
        @type p2: C{Scientific.Geometry.Vector}
        @returns: a Collection of all objects in the universe that lie
                   within the box whose diagonally opposite corners are
                   given by p1 and p2.
        """
        return self._objects.selectBox(p1, p2)

    def _changed(self, system_size_changed):
        self._evaluator = {}
        self._bond_database = None
        self._version += 1
        if system_size_changed:
            for a in self.atomList():
                a.unsetArray()
            self._configuration = None
            self._masses = None
            self._atom_properties = {}
            self._atoms = None
            self._np = None
            self._bond_pairs = None
        else:
            if self._configuration is not None:
                self._configuration.version = self._version
            if self._masses is not None:
                self._masses.version = self._version

    def acquireReadStateLock(self):
        """
        Acquire the universe read state lock. Any application that
        uses threading must acquire this lock prior to accessing the
        current state of the universe, in particular its configuration
        (particle positions). This guarantees the consistency of the
        data; while any thread holds the read state lock, no other
        thread can obtain the write state lock that permits modifying
        the state. The read state lock should be released as soon as
        possible.

        The read state lock can be acquired only if no thread holds
        the write state lock. If the read state lock cannot be
        acquired immediately, the thread will be blocked until
        it becomes available. Any number of threads can acquire
        the read state lock simultaneously.
        """
        return self._spec.stateLock(1)

    def acquireWriteStateLock(self):
        """
        Acquire the universe write state lock. Any application that
        uses threading must acquire this lock prior to modifying the
        current state of the universe, in particular its configuration
        (particle positions). This guarantees the consistency of the
        data; while any thread holds the write state lock, no other
        thread can obtain the read state lock that permits accessing
        the state. The write state lock should be released as soon as
        possible.

        The write state lock can be acquired only if no other thread
        holds either the read state lock or the write state lock. If
        the write state lock cannot be acquired immediately, the
        thread will be blocked until it becomes available.
        """
        return self._spec.stateLock(-1)

    def releaseReadStateLock(self, write=False):
        """
        Release the universe read state lock.
        """
        return self._spec.stateLock(2)

    def releaseWriteStateLock(self, write=False):
        """
        Release the universe write state lock.
        """
        return self._spec.stateLock(-2)

    def acquireConfigurationChangeLock(self, waitflag=True):
        """
        Acquire the configuration change lock. This lock should be
        acquired before starting an algorithm that changes the
        configuration continuously, e.g. minimization or molecular dynamics
        algorithms. This guarantees the proper order of execution when
        several such operations are started in succession. For example,
        when a minimization should be followed by a dynamics run,
        the use of this flag permits both operations to be started
        as background tasks which will be executed one after the other,
        permitting other threads to run in parallel.

        The configuration change lock should not be confused with
        the universe state lock. The former guarantees the proper
        sequence of long-running algorithms, whereas the latter
        guarantees the consistency of the data. A dynamics algorithm,
        for example, keeps the configuration change lock from the
        beginning to the end, but acquires the universe state lock
        only immediately before modifying configuration and velocities,
        and releases it immediately afterwards.

        @param waitflag: if true, the method waits until the lock
                         becomes available; this is the most common mode.
                         If false, the method returns immediately even
                         if another thread holds the lock.
        @type waitflag: C{bool}
        @returns: a flag indicating if the lock was successfully
                  acquired (1) or not (0).
        @rtype: C{int}
        """
        if waitflag:
            return self._spec.configurationChangeLock(1)
        else:
            return self._spec.configurationChangeLock(0)

    def releaseConfigurationChangeLock(self):
        """
        Releases the configuration change lock.
        """
        self._spec.configurationChangeLock(2)

    def setForceField(self, forcefield):
        """
        @param forcefield: the new forcefield for this universe
        @type forcefield: L{MMTK.ForceField.ForceField.ForceField}
        """
        self._forcefield = forcefield
        self._evaluator = {}
        self._bond_database = None

    def position(self, object, conf):
        if ChemicalObjects.isChemicalObject(object):
            return object.position(conf)
        elif isVector(object):
            return object
        else:
            return Vector(object)

    def numberOfAtoms(self):
        return self._objects.numberOfAtoms()

    def numberOfPoints(self):
        if self._np is None:
            self._np = Collections.GroupOfAtoms.numberOfPoints(self)
        return self._np

    numberOfCartesianCoordinates = numberOfPoints

    def configuration(self):
        """
        @returns: the configuration object describing the current
                  configuration of the universe. Note that this is not a
                  copy of the current state, but a reference: the positions
                  in the configuration object will change when coordinate
                  changes are applied to the universe in whatever way.
        @rtype: L{MMTK.ParticleProperties.Configuration}
        """
        if self._configuration is None:
            np = self.numberOfAtoms()
            coordinates = N.zeros((np, 3), N.Float)
            index_map = {}
            redef = []
            for a in self.atomList():
                if a.index is None or a.index >= np:
                    redef.append(a)
                else:
                    if index_map.get(a.index, None) is None:
                        index_map[a.index] = a
                    else:
                        redef.append(a)
            free_indices = [i for i in xrange(np)
                            if index_map.get(i, None) is None]
            assert len(free_indices) == len(redef)
            for a, i in zip(redef, free_indices):
                a.index = i
            # At this point a.index runs from 0 to np-1 in the universe.
            for a in self.atomList():
                if a.array is None:
                    try:
                        coordinates[a.index, :] = a.pos.array
                        del a.pos
                    except AttributeError:
                        coordinates[a.index, :] = Utility.undefined
                else:
                    coordinates[a.index, :] = a.array[a.index, :]
                a.array = coordinates
            # Define configuration object.
            self._configuration = 1 # a hack to prevent endless recursion
            self._configuration = \
                         ParticleProperties.Configuration(self, coordinates)
        return self._configuration

    def copyConfiguration(self):
        """
        This operation is thread-safe; it won't return inconsistent
        data even when another thread is modifying the configuration.

        @returns: a copy of the current configuration
        @rtype: L{MMTK.ParticleProperties.Configuration}
        """
        self.acquireReadStateLock()
        try:
            conf = copy.copy(self.configuration())
        finally:
            self.releaseReadStateLock()
        return conf

    def atomNames(self):
        self.configuration()
        names = self.numberOfAtoms()*[None]
        for a in self.atomList():
            names[a.index] = a.fullName()
        return names

    def setConfiguration(self, configuration, block=True):
        """
        Update the current configuration of the universe by copying
        the given input configuration.
        
        @param configuration: the new configuration
        @type configuration: L{MMTK.ParticleProperties.Configuration}
        @param block: if C{True}, the operation blocks other threads
                      from accessing the configuration before the update
                      is completed. If C{False}, it is assumed that the
                      caller takes care of locking.
        @type block: C{bool}
        """
        if not ParticleProperties.isConfiguration(configuration):
            raise TypeError('not a universe configuration')
        conf = self.configuration()
        if block:
            self.acquireWriteStateLock()
        try:
            conf.assign(configuration)
            self.setCellParameters(configuration.cell_parameters)
        finally:
            if block:
                self.releaseWriteStateLock()

    def addToConfiguration(self, displacement, block=True):
        """
        Update the current configuration of the universe by adding
        the given displacement vector.
        
        @param displacement: the displacement vector for each atom
        @type displacement: L{MMTK.ParticleProperties.ParticleVector}
        @param block: if C{True}, the operation blocks other threads
                      from accessing the configuration before the update
                      is completed. If C{False}, it is assumed that the
                      caller takes care of locking.
        @type block: C{bool}
        """
        conf = self.configuration()
        if block:
            self.acquireWriteStateLock()
        try:
            conf.assign(conf+displacement)
        finally:
            if block:
                self.releaseWriteStateLock()

    def getParticleScalar(self, name, datatype = N.Float):
        """
        @param name: the name of an atom attribute
        @type name: C{str}
        @param datatype: the datatype of the array allocated to hold the data
        @returns: the values of the attribute 'name' for each atom
                  in the universe.
        @rtype: L{MMTK.ParticleProperties.ParticleScalar}
        """
        conf = self.configuration()
        array = N.zeros((len(conf),), datatype)
        for a in self.atomList():
            array[a.index] = getattr(a, name)
        return ParticleProperties.ParticleScalar(self, array)
    getAtomScalarArray = getParticleScalar

    def getParticleBoolean(self, name):
        """
        @param name: the name of an atom attribute
        @type name: C{str}
        @returns: the values of the boolean attribute 'name' for each atom
                  in the universe, or C{False} for atoms that do not have
                  the attribute.
        @rtype: L{MMTK.ParticleProperties.ParticleScalar}
        """
        conf = self.configuration()
        array = N.zeros((len(conf),), N.Int)
        for a in self.atomList():
            try:
                array[a.index] = getattr(a, name)
            except AttributeError: pass
        return ParticleProperties.ParticleScalar(self, array)
    getAtomBooleanArray = getParticleBoolean

    def masses(self):
        """
        @returns: the masses of all atoms in the universe
        @rtype: L{MMTK.ParticleProperties.ParticleScalar}
        """
        if self._masses is None:
            self._masses = self.getParticleScalar('_mass')
        return self._masses

    def charges(self):
        """
        Return the atomic charges defined by the universe's
        force field.
        @returns: the charges of all atoms in the universe
        @rtype: L{MMTK.ParticleProperties.ParticleScalar}
        """
        ff = self._forcefield
        if ff is None:
            raise ValueError("no force field defined")
        return ff.charges(self)

    def velocities(self):
        """
        @returns: the current velocities of all atoms, or C{None} if
                  no velocities are defined. Note that this is not a
                  copy of the current state but a reference to it;
                  its data will change whenever any changes are made
                  to the current velocities.
        @rtype: L{MMTK.ParticleProperties.ParticleVector}
        """
        try:
            return self._atom_properties['velocity']
        except KeyError:
            return None

    def setVelocities(self, velocities, block=True):
        """
        Update the current velocities of the universe by copying
        the given input velocities.
        
        @param velocities: the new velocities, or C{None} to remove
                           the velocity definition from the universe
        @type velocities: L{MMTK.ParticleProperties.ParticleVector}
        @param block: if C{True}, the operation blocks other threads
                      from accessing the configuration before the update
                      is completed. If C{False}, it is assumed that the
                      caller takes care of locking.
        @type block: C{bool}
        """
        if velocities is None:
            try:
                del self._atom_properties['velocity']
            except KeyError:
                pass
        else:
            try:
                v = self._atom_properties['velocity']
            except KeyError:
                v = ParticleProperties.ParticleVector(self)
                self._atom_properties['velocity'] = v
            if block:
                self.acquireWriteStateLock()
            try:
                v.assign(velocities)
            finally:
                if block:
                    self.releaseWriteStateLock()

    def initializeVelocitiesToTemperature(self, temperature):
        """
        Generate random velocities for all atoms from a Boltzmann
        distribution.
        @param temperature: the reference temperature for the Boltzmann
                            distribution
        @type temperature: C{float}
        """
        self.configuration()
        masses = self.masses()
        if self._atom_properties.has_key('velocity'):
            del self._atom_properties['velocity']
        fixed = self.getParticleBoolean('fixed')
        np = self.numberOfPoints()
        velocities = N.zeros((np, 3), N.Float)
        for i in xrange(np):
            m = masses[i]
            if m > 0. and not fixed[i]:
                velocities[i] = Random.randomVelocity(temperature,
                                                           m).array
        self._atom_properties['velocity'] = \
                          ParticleProperties.ParticleVector(self, velocities)
        self.adjustVelocitiesToConstraints()

    def scaleVelocitiesToTemperature(self, temperature, block=True):
        """
        Scale all velocities by a common factor in order to obtain
        the specified temperature.
        @param temperature: the reference temperature
        @type temperature: C{float}
        @param block: if C{True}, the operation blocks other threads
                      from accessing the configuration before the update
                      is completed. If C{False}, it is assumed that the
                      caller takes care of locking.
        @type block: C{bool}
        """
        velocities = self.velocities()
        factor = N.sqrt(temperature/self.temperature())
        if block:
            self.acquireWriteStateLock()
        try:
            velocities.scaleBy(factor)
        finally:
            if block:
                self.releaseWriteStateLock()

    def degreesOfFreedom(self):
        return GroupOfAtoms.degreesOfFreedom(self) \
               - self.numberOfDistanceConstraints()

    def distanceConstraintList(self):
        """
        @returns: the list of distance constraints
        @rtype: C{list}
        """
        return self._objects.distanceConstraintList()

    def numberOfDistanceConstraints(self):
        """
        @returns: the number of distance constraints
        @rtype: C{int}
        """
        return self._objects.numberOfDistanceConstraints()

    def setBondConstraints(self):
        """
        Sets distance constraints for all bonds.
        """
        self.configuration()
        self._objects.setBondConstraints(self)
        self.enforceConstraints()

    def removeDistanceConstraints(self):
        """
        Removes all distance constraints.
        """
        self._objects.removeDistanceConstraints(self)

    def enforceConstraints(self, configuration=None, velocities=None):
        """
        Enforces the previously defined distance constraints
        by modifying the configuration and velocities.
        @param configuration: the configuration in which the
                              constraints are enforced
                              (C{None} for current configuration)
        @type configuration: L{MMTK.ParticleProperties.Configuration}
        @param velocities: the velocities in which the
                              constraints are enforced
                              (C{None} for current velocities)
        @type velocities: L{MMTK.ParticleProperties.ParticleVector}
        """
        from MMTK import Dynamics
        Dynamics.enforceConstraints(self, configuration)
        self.adjustVelocitiesToConstraints(velocities)

    def adjustVelocitiesToConstraints(self, velocities=None, block=True):
        """
        Modifies the velocities to be compatible with
        the distance constraints, i.e. projects out the velocity
        components along the constrained distances.
        @param velocities: the velocities in which the
                              constraints are enforced
                              (C{None} for current velocities)
        @type velocities: L{MMTK.ParticleProperties.ParticleVector}
        @param block: if C{True}, the operation blocks other threads
                      from accessing the configuration before the update
                      is completed. If C{False}, it is assumed that the
                      caller takes care of locking.
        @type block: C{bool}
        """
        from MMTK import Dynamics
        if velocities is None:
            velocities = self.velocities()
        if velocities is not None:
            if block:
                self.acquireWriteStateLock()
            try:
                Dynamics.projectVelocities(self, velocities)
            finally:
                if block:
                    self.releaseWriteStateLock()

    def bondLengthDatabase(self):
        if self._bond_database is None:
            self._bond_database = None
            if self._bond_database is None:
                ff = self._forcefield
                try:
                    self._bond_database = ff.bondLengthDatabase(self)
                except AttributeError:
                    pass
            if self._bond_database is None:
                self._bond_database = Bonds.DummyBondLengthDatabase(self)
        return self._bond_database

    def forcefield(self):
        """
        @returns: the force field
        @rtype: L{MMTK.ForceField.ForceField.ForceField}
        """
        return self._forcefield

    def energyEvaluatorParameters(self, subset1 = None, subset2 = None):
        self.configuration()
        from MMTK.ForceFields import ForceField
        ffdata = ForceField.ForceFieldData()
        return self._forcefield.evaluatorParameters(self, subset1, subset2,
                                                    ffdata)

    def energyEvaluator(self, subset1 = None, subset2 = None,
                        threads=None, mpi_communicator=None):
        if self._forcefield is None:
            raise ValueError("no force field defined")
        try:
            eval = self._evaluator[(subset1, subset2, threads)]
        except KeyError:
            from MMTK.ForceFields import ForceField
            eval = ForceField.EnergyEvaluator(self, self._forcefield,
                                              subset1, subset2,
                                              threads, mpi_communicator)
            self._evaluator[(subset1, subset2, threads)] = eval
        return eval

    def energy(self, subset1 = None, subset2 = None, small_change=False):
        """
        @param subset1: a subset of a universe, or C{None}
        @type subset1: L{MMTK.ChemicalObjects.ChemicalObject}
        @param subset2: a subset of a universe, or C{None}
        @type subset2: L{MMTK.ChemicalObjects.ChemicalObject}
        @param small_change: if C{True}, algorithms optimized for small
                             configurational changes relative to the last
                             evaluation may be used.
        @type small_change: C{bool}
        @returns: the potential energy of interaction between the atoms
                  in subset1 and the atoms in subset2. If subset2 is C{None},
                  the interactions within subset1 are calculated. It both
                  subsets are C{None}, the potential energy of the whole
                  universe is returned.
        @rtype: C{float}
        """
        eval = self.energyEvaluator(subset1, subset2)
        return eval(0, 0, small_change)

    def energyAndGradients(self, subset1 = None, subset2 = None,
                           small_change=False):
        """
        @returns: the energy and the energy gradients
        @rtype: (C{float}, L{MMTK.ParticleProperties.ParticleVector})
        """
        eval = self.energyEvaluator(subset1, subset2)
        return eval(1, 0, small_change)

    def energyAndForceConstants(self, subset1 = None, subset2 = None,
                                small_change=False):
        """
        @returns: the energy and the force constants
        @rtype: (C{float}, L{MMTK.ParticleProperties.SymmetricPairTensor})
        """
        eval = self.energyEvaluator(subset1, subset2)
        e, g, fc = eval(0, 1, small_change)
        return e, fc

    def energyGradientsAndForceConstants(self, subset1 = None, subset2 = None,
                                         small_change=False):
        """
        @returns: the energy, its gradients, and the force constants
        @rtype: (C{float}, L{MMTK.ParticleProperties.ParticleVector},
                 L{MMTK.ParticleProperties.SymmetricPairTensor})
        """
        eval = self.energyEvaluator(subset1, subset2)
        return eval(1, 1, small_change)

    def energyTerms(self, subset1 = None, subset2 = None, small_change=False):
        """
        @returns: a dictionary containing the energy values for each
                  energy term separately. The energy terms are defined by the
                  force field.
        @rtype: C{dict}
        """
        eval = self.energyEvaluator(subset1, subset2)
        eval(0, 0, small_change)
        return eval.lastEnergyTerms()

    def configurationDifference(self, conf1, conf2):
        """
        @param conf1: a configuration
        @type conf1: L{MMTK.ParticleProperties.Configuration}
        @param conf2: a configuration
        @type conf2: L{MMTK.ParticleProperties.Configuration}
        @returns: the difference vector between the two configurations
                  for each atom, taking into account the universe
                  topology (e.g. minimum-image convention).
        @rtype: L{MMTK.ParticleProperties.ParticleVector}
        """
        d = conf2-conf1
        cell = conf1.cell_parameters
        if cell is not None:
            self._spec.foldCoordinatesIntoBox(d.array)
        return d
            
    def distanceVector(self, p1, p2, conf=None):
        """
        @param p1: a vector or a chemical object whose position is taken
        @param p2: a vector or a chemical object whose position is taken
        @param conf: a configuration (C{None} for the current configuration)
        @returns: the distance vector between p1 and p2 (i.e. the
                  vector from p1 to p2) in the configuration conf,
                  taking into account the universe's topology.
        """
        p1 = self.position(p1, conf)
        p2 = self.position(p2, conf)
        if conf is None:
            return Vector(self._spec.distanceVector(p1.array, p2.array))
        else:
            cell = self._fixCellParameters(conf.cell_parameters)
            if cell is None:
                return Vector(self._spec.distanceVector(p1.array, p2.array))
            else:
                return Vector(self._spec.distanceVector(p1.array, p2.array,
                                                        cell))
            
    def distance(self, p1, p2, conf = None):
        """
        @param p1: a vector or a chemical object whose position is taken
        @param p2: a vector or a chemical object whose position is taken
        @param conf: a configuration (C{None} for the current configuration)
        @returns: the distance between p1 and p2, i.e. the length
                  of the distance vector
        @rtype: C{float}
        """
        return self.distanceVector(p1, p2, conf).length()

    def angle(self, p1, p2, p3, conf = None):
        """
        @param p1: a vector or a chemical object whose position is taken
        @param p2: a vector or a chemical object whose position is taken
        @param p3: a vector or a chemical object whose position is taken
        @param conf: a configuration (C{None} for the current configuration)
        @returns: the angle between the distance vectors p1-p2 and p3-p2
        @rtype: C{float}
        """
        v1 = self.distanceVector(p2, p1, conf)
        v2 = self.distanceVector(p2, p3, conf)
        return v1.angle(v2)

    def dihedral(self, p1, p2, p3, p4, conf = None):
        """
        @param p1: a vector or a chemical object whose position is taken
        @param p2: a vector or a chemical object whose position is taken
        @param p3: a vector or a chemical object whose position is taken
        @param p4: a vector or a chemical object whose position is taken
        @param conf: a configuration (C{None} for the current configuration)
        @returns: the dihedral angle between the plane containing the
                  distance vectors p1-p2 and p3-p2 and the plane containing
                  the distance vectors p2-p3 and p4-p3
        @rtype: C{float}
        """
        v1 = self.distanceVector(p2, p1, conf)
        v2 = self.distanceVector(p3, p2, conf)
        v3 = self.distanceVector(p3, p4, conf)
        a = v1.cross(v2).normal()
        b = v3.cross(v2).normal()
        cos = a*b
        sin = b.cross(a)*v2/v2.length()
        return Transformation.angleFromSineAndCosine(sin, cos)

    def _deleteAtom(self, atom):
        pass

    def basisVectors(self):
        """
        @returns: the basis vectors of the elementary cell of a periodic
                  universe, or C{None} for a non-periodic universe
        """
        return None

    def reciprocalBasisVectors(self):
        """
        @returns: the reciprocal basis vectors of the elementary cell of
                  a periodic universe, or C{None} for a non-periodic universe
        """
        return None

    def cellParameters(self):
        return None

    def setCellParameters(self, parameters):
        if parameters is not None:
            raise ValueError('incompatible cell parameters')

    def _fixCellParameters(self, cell_parameters):
        return cell_parameters

    def cellVolume(self):
        """
        @returns: the volume of the elementary cell of a periodic
                  universe, C{None} for a non-periodic universe
        """
        return None

    def largestDistance(self):
        """
        @returns: the largest possible distance between any two points
                  that can be represented independent of orientation, i.e. the
                  radius of the largest sphere that fits into the simulation
                  cell. Returns C{None} if no such upper limit exists.
        """
        return None

    def contiguousObjectOffset(self, objects = None, conf = None,
                               box_coordinates = False):
        """
        @param objects: a list of chemical objects, or C{None} for all
                        objects in the universe
        @type objects: C{list}
        @param conf: a configuration (C{None} for the current configuration)
        @param box_coordinates: use box coordinates rather than real ones
        @type box_coordinates: C{bool}
        @returns: a set of displacement vectors relative to
                  the conf which, when added to the configuration,
                  create a configuration in which none of the objects
                  is split across the edge of the elementary cell.
                  For nonperiodic universes the return value is C{None}.
        @rtype: L{MMTK.ParticleProperties.ParticleVector}
        """
        return None

    def contiguousObjectConfiguration(self, objects = None, conf = None):
        """
        @param objects: a list of chemical objects, or C{None} for all
                        objects in the universe
        @type objects: C{list}
        @param conf: a configuration (C{None} for the current configuration)
        @returns: configuration conf (default: current configuration)
                  corrected by the contiguous object offsets for that
                  configuration.
        @rtype: L{MMTK.ParticleProperties.Configuration}
        """
        if conf is None:
            conf = self.configuration()
        offset = self.contiguousObjectOffset(objects, conf)
        if offset is not None:
            return conf + offset
        else:
            return copy.copy(conf)

    def realToBoxCoordinates(self, vector):
        """
        Box coordinates are defined only for periodic universes;
        their components have values between -0.5 and 0.5; these
        extreme values correspond to the walls of the simulation box.
        @param vector: a point in the universe
        @returns: the box coordinate equivalent of vector, or the original
                  vector if no box coordinate system exists
        @rtype: C{Scientific.Geometry.Vector}
        """
        return vector
    
    def boxToRealCoordinates(self, vector):
        """
        @param vector: a point in the universe expressed in box coordinates
        @returns: the real-space equivalent of vector
        @rtype: C{Scientific.Geometry.Vector}
        """
        return vector

    def _realToBoxPointArray(self, array, parameters=None):
        return array

    def _boxToRealPointArray(self, array, parameters=None):
        return array

    def cartesianToFractional(self, vector):
        """
        Fractional coordinates are defined only for periodic universes;
        their components have values between 0. and 1.

        @param vector: a point in the universe
        @type vector: C{Scientific.Geometry.Vector}
        @returns: the fractional coordinate equivalent of vector
        @rtype: C{Scientific.N.array_type}
        """
        raise ValueError("Universe is not periodic")

    def cartesianToFractionalMatrix(self):
        raise ValueError("Universe is not periodic")

    def fractionalToCartesian(self, array):
        """
        Fractional coordinates are defined only for periodic universes;
        their components have values between 0. and 1.

        @param array: an array of fractional coordinates
        @type array: C{Scientific.N.array_type}
        @returns: the real-space equivalent of vector
        @rtype: C{Scientific.Geometry.Vector}
        """
        raise ValueError("Universe is not periodic")

    def fractionalToCartesianMatrix(self):
        raise ValueError("Universe is not periodic")

    def foldCoordinatesIntoBox(self):
        return

    def randomPoint(self):
        """
        @returns: a random point from a uniform distribution within
                  the universe. This operation is defined only for
                  finite-volume universes, e.g. periodic universes.
        @rtype: C{Scientific.Geometry.Vector}
        """
        raise TypeError("undefined operation")

    def map(self, function):
        """
        Apply a function to all objects in the universe and
        return the list of the results. If the results are chemical
        objects, a Collection object is returned instead of a list.

        @param function: the function to be applied
        @type function: callable
        @returns: the list or collection of the results
        """
        return self._objects.map(function)

    def description(self, objects = None, index_map = None):
        if objects is None:
            objects = self
        attributes = {}
        for attr in dir(self):
            if attr[0] != '_':
                object = getattr(self, attr)
                if ChemicalObjects.isChemicalObject(object) \
                   or Environment.isEnvironmentObject(object):
                    attributes[object] = attr
        items = []
        for o in objects.objectList():
            attr = attributes.get(o, None)
            if attr is not None:
                items.append(repr(attr))
            items.append(o.description(index_map))
        for o in self._environment:
            attr = attributes.get(o, None)
            if attr is not None:
                items.append(repr(attr))
            items.append(o.description())
        try:
            classname = self.classname_for_trajectories
        except AttributeError:
            classname = self.__class__.__name__
        s = 'c(%s,[%s])' % \
            (`classname + self._descriptionArguments()`,
             ','.join(items))
        return s

    def _graphics(self, conf, distance_fn, model, module, options):
        return self._objects._graphics(conf, distance_fn, model,
                                       module, options)

    def setFromTrajectory(self, trajectory, step = None):
        """
        Set the state of the universe to the one stored in a trajectory.
        This operation is thread-safe; it blocks other threads that
        want to access the configuration or velocities while the data is
        being updated.

        @param trajectory: a trajectory object for this universe
        @type trajectory: L{MMTK.Trajectory.Trajectory}
        @param step: a step number, or C{None} for the default step
                     (0 for a standard trajectory, the last written
                     step for a restart trajectory)
        @type step: C{int}
        """
        if step is None:
            step = trajectory.defaultStep()
        self.acquireWriteStateLock()
        try:
            self.setConfiguration(trajectory.configuration[step], False)
            vel = self.velocities()
            try:
                vel_tr = trajectory.velocities[step]
            except AttributeError:
                if vel is not None:
                    Utility.warning("velocities were not modified because " +
                                    "the trajectory does not contain " +
                                    "velocity data.")
                return
            if vel is None:
                self._atom_properties['velocity'] = vel_tr
            else:
                vel.assign(vel_tr)
        finally:
            self.releaseWriteStateLock()

    #
    # More efficient reimplementations of methods in Collections.GroupOfAtoms
    #
    def numberOfFixedAtoms(self):
        return self.getParticleBoolean('fixed').sumOverParticles()

    def degreesOfFreedom(self):
        return 3*(self.numberOfAtoms()-self.numberOfFixedAtoms()) \
               - self.numberOfDistanceConstraints()

    def mass(self):
        return self.masses().sumOverParticles()

    def centerOfMass(self, conf = None):
        m = self.masses()
        if conf is None:
            conf = self.configuration()
        return (m*conf).sumOverParticles()/m.sumOverParticles()

    def kineticEnergy(self, velocities = None):
        if velocities is None:
            velocities = self.velocities()
        return 0.5*velocities.massWeightedDotProduct(velocities)

    def momentum(self, velocities = None):
        if velocities is None:
            velocities = self.velocities()
        return (self.masses()*velocities).sumOverParticles()

    def translateBy(self, vector):
        conf = self.configuration().array
        N.add(conf, vector.array[N.NewAxis, :], conf)

    def applyTransformation(self, t):
        conf = self.configuration().array
        rot = t.rotation().tensor.array
        conf[:] = N.dot(conf, N.transpose(rot))
        N.add(conf, t.translation().vector.array[N.NewAxis, :], conf)

    def writeXML(self, file):
        file.write('<?xml version="1.0" encoding="ISO-8859-1" ' +
                   'standalone="yes"?>\n\n')
        file.write('<molecularsystem>\n\n')
        file.write('<templates>\n\n')
        memo = {'counter': 1}
        instances = []
        atoms = []
        for object in self._objects.objectList():
            instances = instances + object.writeXML(file, memo, 1)
            atoms = atoms + object.getXMLAtomOrder()
        file.write('\n</templates>\n\n')
        file.write('<universe %s>\n' % self.XMLSpec())
        for instance in instances:
            file.write('  ')
            file.write(instance)
            file.write('\n')
        conf = self.configuration()
        if conf.hasValidPositions():
            file.write('  <configuration>\n')
            file.write('    <atomArray units="units:nm"\n')
            file.write('    x3="')
            for atom in atoms:
                file.write(str(conf[atom][0]))
                file.write(' ')
            file.write('"\n')
            file.write('    y3="')
            for atom in atoms:
                file.write(str(conf[atom][1]))
                file.write(' ')
            file.write('"\n')
            file.write('    z3="')
            for atom in atoms:
                file.write(str(conf[atom][2]))
                file.write(' ')
            file.write('"\n')
            file.write('    />\n')
            file.write('  </configuration>\n')
        file.write('</universe>\n\n') 
        file.write('</molecularsystem>\n')
       

#
# Infinite universes
#
class InfiniteUniverse(Universe):

    """
    Infinite (unbounded and nonperiodic) universe.
    """

    def __init__(self, forcefield=None, **properties):
        """
        @param forcefield: a force field, or C{None} for no force field
        @type forcefield: L{MMTK.ForceField.ForceField.ForceField}
        """
        Universe.__init__(self, forcefield, properties)
        self._createSpec()

    def CdistanceFunction(self):
        from MMTK_universe import infinite_universe_distance_function
        return infinite_universe_distance_function, N.array([0.])

    def CcorrectionFunction(self):
        from MMTK_universe import infinite_universe_correction_function
        return infinite_universe_correction_function, N.array([0.])

    def CvolumeFunction(self):
        from MMTK_universe import infinite_universe_volume_function
        return infinite_universe_volume_function, N.array([0.])

    def CboxTransformationFunction(self):
        return None, N.array([0.])

    def _createSpec(self):
        from MMTK_universe import InfiniteUniverseSpec
        self._spec = InfiniteUniverseSpec()

    def _descriptionArguments(self):
        if self._forcefield is None:
            return '()'
        else:
            return '(%s)' % self._forcefield.description()

    def XMLSpec(self):
        return 'topology="infinite"'

#
# 3D periodic universe base class
#
class Periodic3DUniverse(Universe):

    is_periodic = True

    def setVolume(self, volume):
        """
        Multiplies all edge lengths by the same factor such that the cell
        volume becomes equal to the specified value.
        @param volume: the desired volume
        @type volume: C{float}
        """
        factor = (volume/self.cellVolume())**(1./3.)
        self.scaleSize(factor)

    def foldCoordinatesIntoBox(self):
        self._spec.foldCoordinatesIntoBox(self.configuration().array)
    
    def basisVectors(self):
        return [self.boxToRealCoordinates(Vector(1., 0., 0.)),
                self.boxToRealCoordinates(Vector(0., 1., 0.)),
                self.boxToRealCoordinates(Vector(0., 0., 1.))]

    def cartesianToFractional(self, vector):
        r1, r2, r3 = self.reciprocalBasisVectors()
        return N.array([r1*vector, r2*vector, r3*vector])

    def cartesianToFractionalMatrix(self):
        return N.array(self.reciprocalBasisVectors())

    def fractionalToCartesian(self, array):
        e1, e2, e3 = self.basisVectors()
        return array[0]*e1 + array[1]*e2 + array[2]*e3

    def fractionalToCartesianMatrix(self):
        return N.transpose(self.basisVectors())

    def randomPoint(self):
        return self.boxToRealCoordinates(Random.randomPointInBox(1., 1., 1.))

    def contiguousObjectOffset(self, objects = None, conf = None,
                               box_coordinates = 0):
        from MMTK_universe import contiguous_object_offset
        if objects is None or objects == self or objects == [self]:
            default = True
            objects = self._objects.objectList()
            pairs = self._bond_pairs
        else:
            default = False
            pairs = None
        if conf is None:
            conf = self.configuration()
        cell = self._fixCellParameters(conf.cell_parameters)
        offset = ParticleProperties.ParticleVector(self)
        if pairs is None:
            pairs = []
            for o in objects:
                new_object = True
                if ChemicalObjects.isChemicalObject(o):
                    units = o.bondedUnits()
                elif Collections.isCollection(o) or isUniverse(o):
                    units = set([u
                                 for element in o
                                 for u in element.topLevelChemicalObject()
                                                              .bondedUnits()])
                else:
                    raise ValueError(str(o) + " not a chemical object")
                for bu in units:
                    atoms = [a.index for a in bu.atomsWithDefinedPositions()]
                    mpairs = bu.traverseBondTree(lambda a: a.index)
                    mpairs = [(a1, a2) for (a1, a2) in mpairs
                              if a1 in atoms and a2 in atoms]
                    if len(mpairs) == 0:
                        mpairs = Utility.pairs(atoms)
                    new_object = False
                    pairs.extend(mpairs)
            pairs = N.array(pairs)
            if default:
                self._bond_pairs = pairs
        if cell is None:
            contiguous_object_offset(self._spec, pairs, conf.array,
                                     offset.array, box_coordinates)
        else:
            contiguous_object_offset(self._spec, pairs, conf.array,
                                     offset.array, box_coordinates, cell)
        return offset

    def _graphics(self, conf, distance_fn, model, module, options):
        objects = self._objects._graphics(conf, distance_fn, model,
                                          module, options)
        v1, v2, v3 = self.basisVectors()
        p = -0.5*(v1+v2+v3)
        color = options.get('color', 'white')
        material = module.EmissiveMaterial(color)
        objects.append(module.Line(p, p+v1, material=material))
        objects.append(module.Line(p, p+v2, material=material))
        objects.append(module.Line(p+v1, p+v1+v2, material=material))
        objects.append(module.Line(p+v2, p+v1+v2, material=material))
        objects.append(module.Line(p, p+v3, material=material))
        objects.append(module.Line(p+v1, p+v1+v3, material=material))
        objects.append(module.Line(p+v2, p+v2+v3, material=material))
        objects.append(module.Line(p+v1+v2, p+v1+v2+v3, material=material))
        objects.append(module.Line(p+v3, p+v1+v3, material=material))
        objects.append(module.Line(p+v3, p+v2+v3, material=material))
        objects.append(module.Line(p+v1+v3, p+v1+v2+v3, material=material))
        objects.append(module.Line(p+v2+v3, p+v1+v2+v3, material=material))
        return objects

#
# Orthorhombic universe with periodic boundary conditions
#
class OrthorhombicPeriodicUniverse(Periodic3DUniverse):

    """
    Periodic universe with orthorhombic elementary cell.
    """

    def __init__(self, size = None, forcefield = None, **properties):
        """
        @param size: a sequence of length three specifying the edge
                     lengths along the x, y, and z directions
        @param forcefield: a force field, or C{None} for no force field
        @type forcefield: L{MMTK.ForceField.ForceField.ForceField}
        """
        Universe.__init__(self, forcefield, properties)
        self.data = N.zeros((3,), N.Float)
        if size is not None:
            self.setSize(size)
        self._createSpec()

    is_orthogonal = True

    def __setstate__(self, state):
        Universe.__setstate__(self, state)
        if len(self.data.shape) == 2:
            self.data = self.data[0]

    def setSize(self, size):
        self.data[:] = size

    def scaleSize(self, factor):
        """
        Multiplies all edge lengths by a factor.
        @param factor: the scale factor
        @type factor: C{float}
        """
        self.data[:] = factor*self.data
        self._spec.foldCoordinatesIntoBox(self.configuration().array)

    def setCellParameters(self, parameters):
        if parameters is not None:
            self.data[:] = parameters

    def realToBoxCoordinates(self, vector):
        x, y, z = vector
        return Vector(x/self.data[0],
                      y/self.data[1],
                      z/self.data[2])

    def boxToRealCoordinates(self, vector):
        x, y, z = vector
        return Vector(x*self.data[0],
                      y*self.data[1],
                      z*self.data[2])

    def _realToBoxPointArray(self, array, parameters=None):
        if parameters is None:
            parameters = self.data
        if parameters.shape == (3,):
            parameters = parameters[N.NewAxis, :]
        return array/parameters

    def _boxToRealPointArray(self, array, parameters=None):
        if parameters is None:
            parameters = self.data
        if parameters.shape == (3,):
            parameters = parameters[N.NewAxis, :]
        return array*parameters

    def CdistanceFunction(self):
        from MMTK_universe import orthorhombic_universe_distance_function
        return orthorhombic_universe_distance_function, self.data

    def CcorrectionFunction(self):
        from MMTK_universe import orthorhombic_universe_correction_function
        return orthorhombic_universe_correction_function, self.data

    def CvolumeFunction(self):
        from MMTK_universe import orthorhombic_universe_volume_function
        return orthorhombic_universe_volume_function, self.data

    def CboxTransformationFunction(self):
        from MMTK_universe import orthorhombic_universe_box_transformation
        return orthorhombic_universe_box_transformation, self.data

    def cellParameters(self):
        return self.data

    def reciprocalBasisVectors(self):
        return [Vector(1., 0., 0.)/self.data[0],
                Vector(0., 1., 0.)/self.data[1],
                Vector(0., 0., 1.)/self.data[2]]

    def cellVolume(self):
        return N.multiply.reduce(self.data)

    def largestDistance(self):
        return 0.5*N.minimum.reduce(self.data)

    def _createSpec(self):
        from MMTK_universe import OrthorhombicPeriodicUniverseSpec
        self._spec = OrthorhombicPeriodicUniverseSpec(self.data)

    def _descriptionArguments(self):
        if self._forcefield is None:
            return '((0.,0.,0.),)'
        else:
            return '((0.,0.,0.),%s)' % self._forcefield.description()

    def XMLSpec(self):
        return 'topology="periodic3d" ' + \
               'cellshape="orthorhombic" ' + \
               ('cellsize="%f %f %f" ' % tuple(self.data)) + \
               'units="units:nm"'

#
# Cubic universe with periodic boundary conditions
#
class CubicPeriodicUniverse(OrthorhombicPeriodicUniverse):

    """
    Periodic universe with cubic elementary cell.
    """

    def setSize(self, size):
        """
        Set the edge length to a given value.
        @param size: the new size
        @type size: C{float}
        """
        OrthorhombicPeriodicUniverse.setSize(self, 3*(size,))

    def _descriptionArguments(self):
        if self._forcefield is None:
            return '(0.)'
        else:
            return '(0.,%s)' % self._forcefield.description()

#
# Parallelepipedic universe with periodic boundary conditions
#
class ParallelepipedicPeriodicUniverse(Periodic3DUniverse):

    """
    Periodic universe with parallelepipedic elementary cell.
    """

    def __init__(self, shape = None, forcefield = None, **properties):
        """
        @param shape: the basis vectors
        @type shape: sequence of C{Scientific.Geometry.Vector}
        @param forcefield: a force field, or C{None} for no force field
        @type forcefield: L{MMTK.ForceField.ForceField.ForceField}
        """
        Universe.__init__(self, forcefield, properties)
        self.data = N.zeros((19,), N.Float)
        if shape is not None:
            self.setShape(shape)
        self._createSpec()

    is_periodic = True

    def setShape(self, shape):
        self.data[:9] = N.ravel(N.transpose([list(s) for s in shape]))
        from MMTK_universe import parallelepiped_invert
        parallelepiped_invert(self.data)

    def scaleSize(self, factor):
        """
        Multiplies all edge lengths by a factor.
        @param factor: the scale factor
        @type factor: C{float}
        """
        self.data[:9] = factor*self.data[:9]
        from MMTK_universe import parallelepiped_invert
        parallelepiped_invert(self.data)
        self._spec.foldCoordinatesIntoBox(self.configuration().array)

    def setCellParameters(self, parameters):
        if parameters is not None:
            self.data[:9] = parameters
            from MMTK_universe import parallelepiped_invert
            parallelepiped_invert(self.data)

    def _fixCellParameters(self, cell_parameters):
        full_parameters = 0.*self.data
        full_parameters[:9] = cell_parameters
        from MMTK_universe import parallelepiped_invert
        parallelepiped_invert(full_parameters)
        return full_parameters

    def realToBoxCoordinates(self, vector):
        x, y, z = vector
        return Vector(self.data[0+9]*x + self.data[1+9]*y + self.data[2+9]*z,
                      self.data[3+9]*x + self.data[4+9]*y + self.data[5+9]*z,
                      self.data[6+9]*x + self.data[7+9]*y + self.data[8+9]*z)

    def boxToRealCoordinates(self, vector):
        x, y, z = vector
        return Vector(self.data[0]*x + self.data[1]*y + self.data[2]*z,
                      self.data[3]*x + self.data[4]*y + self.data[5]*z,
                      self.data[6]*x + self.data[7]*y + self.data[8]*z)

    def _realToBoxPointArray(self, array, parameters=None):
        if parameters is None:
            matrix = N.reshape(self.data[9:18], (1, 3, 3))
        else:
            parameters = N.concatenate([parameters, N.zeros((10,), N.Float)])
            from MMTK_universe import parallelepiped_invert
            parallelepiped_invert(parameters)
            matrix = N.reshape(parameters[9:18], (1, 3, 3))
        return N.add.reduce(matrix*array[:, N.NewAxis, :], axis=-1)

    def _boxToRealPointArray(self, array, parameters=None):
        if parameters is None:
            parameters = self.data[:9]
        matrix = N.reshape(parameters, (1, 3, 3))
        return N.add.reduce(matrix*array[:, N.NewAxis, :], axis=-1)

    def CdistanceFunction(self):
        from MMTK_universe import parallelepipedic_universe_distance_function
        return parallelepipedic_universe_distance_function, self.data

    def CcorrectionFunction(self):
        from MMTK_universe import parallelepipedic_universe_correction_function
        return parallelepipedic_universe_correction_function, self.data

    def CvolumeFunction(self):
        from MMTK_universe import parallelepipedic_universe_volume_function
        return parallelepipedic_universe_volume_function, self.data

    def CboxTransformationFunction(self):
        from MMTK_universe import parallelepipedic_universe_box_transformation
        return parallelepipedic_universe_box_transformation, self.data

    def cellParameters(self):
        return self.data[:9]

    def reciprocalBasisVectors(self):
        return [Vector(self.data[9:12]),
                Vector(self.data[12:15]),
                Vector(self.data[15:18])]

    def cellVolume(self):
        return abs(self.data[18])

    def largestDistance(self):
        return min([0.5/v.length() for v in self.reciprocalBasisVectors()])

    def _createSpec(self):
        from MMTK_universe import ParallelepipedicPeriodicUniverseSpec
        self._spec = ParallelepipedicPeriodicUniverseSpec(self.data)

    def _descriptionArguments(self):
        if self._forcefield is None:
            return '((Vector(0.,0.,0.),Vector(0.,0.,0.),Vector(0.,0.,0.)))'
        else:
            return '((Vector(0.,0.,0.),Vector(0.,0.,0.),Vector(0.,0.,0.)),%s)'\
                   % self._forcefield.description()

    def XMLSpec(self):
        return 'topology="periodic3d" ' + \
               'cellshape="parallelepipedic" ' + \
               ('cellshape="%f %f %f %f %f %f %f %f %f" '
                % tuple(self.data[:9])) + \
               'units="units:nm"'

#
# Recognition functions
#
def isUniverse(object):
    """
    @param object: any Python object
    @returns: C{True} if object is a universe.
    """
    return isinstance(object, Universe)
