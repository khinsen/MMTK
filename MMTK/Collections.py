# This module defines collections of chemical objects.
#
# Written by Konrad Hinsen
# last revision: 2010-9-13
#

"""
Collections of chemical objects
"""

__docformat__ = 'epytext'

from MMTK import Utility, Units, ParticleProperties, Visualization
from Scientific.Geometry import Vector, Tensor, Objects3D
from Scientific.Geometry import Quaternion, Transformation
from Scientific import N
import copy, itertools

#
# This class defines groups of atoms. It is used as a base class
# for anything containing atoms, including chemical objects, collections,
# universes etc., but it can't be used directly. All its subclasses
# must define a method atomList() that returns a list of all their atoms.
#
class GroupOfAtoms(object):

    """
    Anything that consists of atoms

    A mix-in class that defines a large set of operations which are
    common to all objects that consist of atoms, i.e. any subset of
    a chemical system. Examples are atoms, molecules, collections,
    or universes.
    """
    
    def numberOfAtoms(self):
        """
        @returns: the number of atoms
        @rtype: C{int}
        """
        return len(self.atomList())

    def numberOfPoints(self):
        """
        @returns: the number of geometrical points that define the
                  object. It is currently equal to the number of atoms,
                  but could be different e.g. for quantum systems, in which
                  each atom is described by a wave function or a path integral.
        @rtype: C{int}
        """
        return sum([a.numberOfPoints() for a in self.atomIterator()])

    numberOfCartesianCoordinates = numberOfPoints

    def beadIterator(self):
        for a in self.atomIterator():
            for b in a.beads():
                yield b

    def numberOfFixedAtoms(self):
        """
        @returns: the number of atoms that are fixed, i.e. that cannot move
        @rtype: C{int}
        """
        n = 0
        for a in self.atomIterator():
            try:
                if a.fixed: n = n + 1
            except AttributeError: pass
        return n

    def degreesOfFreedom(self):
        """
        @returns: the number of mechanical degrees of freedom
        @rtype: C{int}
        """
        return 3*(self.numberOfAtoms()-self.numberOfFixedAtoms())

    def atomCollection(self):
        """
        @returns: a collection containing all atoms in the object
        """
        return Collection(self.atomList())

    def atomsWithDefinedPositions(self, conf = None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: a collection of all atoms that have a position in the
                  given configuration
        """
        return Collection([a for a in self.atomIterator()
                           if Utility.isDefinedPosition(a.position(conf))])

    def mass(self):
        """
        @returns: the total mass
        @rtype: C{float}
        """
        return sum(a._mass for a in self.atomIterator())

    def centerOfMass(self, conf = None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: the center of mass in the given configuration
        @rtype: C{Scientific.Geometry.Vector}
        """
        offset = None
        universe = self.universe()
        if universe is not None:
            offset = universe.contiguousObjectOffset([self], conf)
        m = 0.
        mr = Vector(0.,0.,0.)
        if offset is None:
            for a in self.atomIterator():
                m += a._mass
                mr += a._mass * a.position(conf)
        else:
            for a in self.atomIterator():
                m += a._mass
                mr += a._mass * (a.position(conf)+offset[a])
        return mr/m

    position = centerOfMass

    def centerAndMomentOfInertia(self, conf = None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: the center of mass and the moment of inertia tensor
                  in the given configuration
        """
        from Scientific.Geometry import delta
        offset = None
        universe = self.universe()
        if universe is not None:
            offset = universe.contiguousObjectOffset([self], conf)
        m = 0.
        mr = Vector(0.,0.,0.)
        t = Tensor(3*[3*[0.]])
        for a in self.atomIterator():
            ma = a._mass
            if offset is None:
                r = a.position(conf)
            else:
                r = a.position(conf)+offset[a]
            m += ma
            mr += ma*r
            t += ma*r.dyadicProduct(r)
        cm = mr/m
        t -= m*cm.dyadicProduct(cm)
        t = t.trace()*delta - t
        return cm, t

    def rotationalConstants(self, conf=None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: a sorted array of rotational constants A, B, C
                  in internal units
        """
        com, i = self.centerAndMomentOfInertia(conf)
        pmi = i.eigenvalues()
        return N.sort(Units.h / (8.*N.pi*N.pi*pmi))[::-1]

    def boundingBox(self, conf = None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: two opposite corners of a bounding box around the
                  object. The bounding box is the smallest rectangular
                  bounding box with edges parallel to the coordinate axes.
        @rtype: C{tuple} of two C{Scientific.Geometry.Vector}
        """
        atoms = self.atomList()
        min = atoms[0].beads()[0].position(conf).array
        max = min
        for a in atoms:
            for b in a.beads():
                r = b.position(conf).array
                min = N.minimum(min, r)
                max = N.maximum(max, r)
        return Vector(min), Vector(max)

    def boundingSphere(self, conf = None):
        """
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration} or C{NoneType}
        @returns: a sphere that contains all atoms in the object.
                  This is B{not} the minimal bounding sphere, just B{some}
                  bounding sphere.
        @rtype: Scientific.Geometry.Objects3D.Sphere
        """
        atoms = self.atomList()
        center = sum((a.position(conf) for a in atoms),
                     Vector(0., 0., 0.)) / len(atoms)
        r = 0.
        for p in self.beadIterator():
            r = max(r, (p.position(conf)-center).length())
        return Objects3D.Sphere(center, r)

    def rmsDifference(self, conf1, conf2 = None):
        """
        @param conf1: a configuration object
        @type conf1: L{MMTK.Configuration}
        @param conf2: a configuration object, or C{None} for the
                      current configuration
        @type conf2: L{MMTK.Configuration} or C{NoneType}
        @returns: the RMS (root-mean-square) difference between the
                  conformations of the object in two universe configurations,
                  conf1 and conf2
        @rtype: C{float}
        """
        universe = conf1.universe
        m = 0.
        rms = 0.
        for a in self.atomIterator():
            ma = a._mass
            dr = universe.distanceVector(a.position(conf1), a.position(conf2))
            m += ma
            rms += ma*dr*dr
        return N.sqrt(rms/m)

    def findTransformationAsQuaternion(self, conf1, conf2 = None):
        universe = self.universe()
        if conf1.universe != universe:
            raise ValueError("conformation is for a different universe")
        if conf2 is None:
            conf1, conf2 = conf2, conf1
        else:
            if conf2.universe != universe:
                raise ValueError("conformation is for a different universe")
        ref = conf1
        conf = conf2
        weights = universe.masses()
        weights = weights/self.mass()
        ref_cms = self.centerOfMass(ref).array
        pos = N.zeros((3,), N.Float)
        possq = 0.
        cross = N.zeros((3, 3), N.Float)
        for a in self.atomIterator():
            r = a.position(conf).array
            r_ref = a.position(ref).array-ref_cms
            w = weights[a]
            pos = pos + w*r
            possq = possq + w*N.add.reduce(r*r) \
                          + w*N.add.reduce(r_ref*r_ref)
            cross = cross + w*r[:, N.NewAxis]*r_ref[N.NewAxis, :]
        k = N.zeros((4, 4), N.Float)
        k[0, 0] = -cross[0, 0]-cross[1, 1]-cross[2, 2]
        k[0, 1] = cross[1, 2]-cross[2, 1]
        k[0, 2] = cross[2, 0]-cross[0, 2]
        k[0, 3] = cross[0, 1]-cross[1, 0]
        k[1, 1] = -cross[0, 0]+cross[1, 1]+cross[2, 2]
        k[1, 2] = -cross[0, 1]-cross[1, 0]
        k[1, 3] = -cross[0, 2]-cross[2, 0]
        k[2, 2] = cross[0, 0]-cross[1, 1]+cross[2, 2]
        k[2, 3] = -cross[1, 2]-cross[2, 1]
        k[3, 3] = cross[0, 0]+cross[1, 1]-cross[2, 2]
        for i in range(1, 4):
            for j in range(i):
                k[i, j] = k[j, i]
        k = 2.*k
        for i in range(4):
            k[i, i] = k[i, i] + possq - N.add.reduce(pos*pos)
        from Scientific import LA
        e, v = LA.eigenvectors(k)
        i = N.argmin(e)
        v = v[i]
        if v[0] < 0: v = -v
        if e[i] <= 0.:
            rms = 0.
        else:
            rms = N.sqrt(e[i])
        return Quaternion.Quaternion(v), Vector(ref_cms), \
               Vector(pos), rms

    def findTransformation(self, conf1, conf2 = None):
        """
        @param conf1: a configuration object
        @type conf1: L{MMTK.Configuration}
        @param conf2: a configuration object, or C{None} for the
                      current configuration
        @type conf2: L{MMTK.Configuration} or C{NoneType}
        @returns: the linear transformation that, when applied to
                  the object in configuration conf1, minimizes the
                  RMS distance to the conformation in conf2, and the
                  minimal RMS distance.
                  If conf2 is C{None}, returns the transformation from the
                  current configuration to conf1 and the associated
                  RMS distance.
        """
        q, cm1, cm2, rms = self.findTransformationAsQuaternion(conf1, conf2)
        return Transformation.Translation(cm2) * \
               q.asRotation() * \
               Transformation.Translation(-cm1), \
               rms

    def translateBy(self, vector):
        """
        Translate the object by a displacement vector
        
        @param vector: the displacement vector
        @type vector: C{Scientific.Geometry.Vector}
        """
        for b in self.beadIterator():
            b.translateBy(vector)

    def translateTo(self, position):
        """
        Translate the object such that its center of mass is at position
        @param position: the final position
        @type position: C{Scientific.Geometry.Vector}
        """
        self.translateBy(position-self.centerOfMass())

    def normalizePosition(self):
        """
        Translate the center of mass to the coordinate origin
        """
        self.translateTo(Vector(0., 0., 0.))

    def normalizeConfiguration(self, repr=None):
        """
        Apply a linear transformation such that the center of mass of
        the object is translated to the coordinate origin and its
        principal axes of inertia become parallel to the three
        coordinate axes.

        @param repr: the specific representation for axis alignment:
          - Ir    : x y z <--> b c a
          - IIr   : x y z <--> c a b
          - IIIr  : x y z <--> a b c
          - Il    : x y z <--> c b a
          - IIl   : x y z <--> a c b
          - IIIl  : x y z <--> b a c        
        """
        transformation = self.normalizingTransformation(repr)
        self.applyTransformation(transformation)

    def normalizingTransformation(self, repr=None):
        """
        Calculate a linear transformation that shifts the center of mass
        of the object to the coordinate origin and makes its
        principal axes of inertia parallel to the three coordinate
        axes.

        @param repr: the specific representation for axis alignment:
          Ir    : x y z <--> b c a
          IIr   : x y z <--> c a b
          IIIr  : x y z <--> a b c
          Il    : x y z <--> c b a
          IIl   : x y z <--> a c b
          IIIl  : x y z <--> b a c
        @returns: the normalizing transformation
        @rtype: C{Scientific.Geometry.Transformation.RigidBodyTransformation}
        """
        from Scientific.LA import determinant
        cm, inertia = self.centerAndMomentOfInertia()
        ev, diag = inertia.diagonalization()
        if determinant(diag.array) < 0:
            diag.array[0] = -diag.array[0]
        if repr != None:
            seq = N.argsort(ev)
            if repr == 'Ir':
                seq = N.array([seq[1], seq[2], seq[0]])
            elif repr == 'IIr':
                seq = N.array([seq[2], seq[0], seq[1]])
            elif repr == 'Il':
                seq = N.seq[2::-1]
            elif repr == 'IIl':
                seq[1:3] = N.array([seq[2], seq[1]])
            elif repr == 'IIIl':
                seq[0:2] = N.array([seq[1], seq[0]])
            elif repr != 'IIIr':
                print 'unknown representation'
            diag.array = N.take(diag.array, seq)                
        return Transformation.Rotation(diag)*Transformation.Translation(-cm)

    def applyTransformation(self, t):
        """
        Apply a transformation to the object

        @param t: the transformation to be applied
        @type t: C{Scientific.Geometry.Transformation}
        """
        for b in self.beadIterator():
            b.setPosition(t(b.position()))

    def displacementUnderTransformation(self, t):
        """
        @param t: the transformation to be applied
        @type t: C{Scientific.Geometry.Transformation}
        @returns: the displacement vectors for the atoms in the object
                  that correspond to the transformation |t|.
        @rtype: L{MMTK.ParticleVector}
        """
        d = ParticleProperties.ParticleVector(self.universe())
        for b in self.beadIterator():
            r = b.position()
            d[b] = t(r)-r
        return d

    def rotateAroundCenter(self, axis_direction, angle):
        """
        Rotate the object around an axis that passes through its center
        of mass.

        @param axis_direction: the direction of the axis of rotation
        @type axis_direction: C{Scientific.Geometry.Vector}
        @param angle: the rotation angle (in radians)
        @type angle: C{float}
        """
        cm = self.centerOfMass()
        t = Transformation.Translation(cm) * \
            Transformation.Rotation(axis_direction, angle) * \
            Transformation.Translation(-cm)
        self.applyTransformation(t)

    def rotateAroundOrigin(self, axis_direction, angle):
        """
        Rotate the object around an axis that passes through the
        coordinate origin.

        @param axis_direction: the direction of the axis of rotation
        @type axis_direction: C{Scientific.Geometry.Vector}
        @param angle: the rotation angle (in radians)
        @type angle: C{float}
        """
        self.applyTransformation(Transformation.Rotation(axis_direction, angle))

    def rotateAroundAxis(self, point1, point2, angle):
        """
        Rotate the object arond an axis specified by two points

        @param point1: the first point
        @type point1: C{Scientific.Geometry.Vector}
        @param point2: the second point
        @type point2: C{Scientific.Geometry.Vector}
        @param angle: the rotation angle (in radians)
        @type angle: C{float}
        """
        tr1 = Transformation.Translation(-point1)
        tr2 = Transformation.Rotation(point2-point1, angle)
        tr3 = tr1.inverse()
        self.applyTransformation(tr3*tr2*tr1)

    def writeToFile(self, filename, configuration = None, format = None):
        """
        Write a representation of the object in a given
        configuration to a file.

        @param filename: the name of the file
        @type filename: C{str}
        @param configuration: a configuration object, or C{None} for the
                              current configuration
        @type configuration: L{MMTK.Configuration} or C{NoneType}
        @param format: 'pdb' or 'vrml' (default: guess from filename)
                       A subformat specification can be added, separated
                       by a dot. Subformats of 'vrml' are 'wireframe'
                       (default), 'ball_and_stick', 'highlight' (like
                       'wireframe', but with a small sphere for
                       all atoms that have an attribute 'highlight' with a
                       non-zero value), and 'charge' (wireframe plus small
                       spheres for the atoms whose color indicates the
                       charge on a red-to-green color scale)
        @type format: C{str}
        """
        from MMTK import ConfigIO
        universe = self.universe()
        if universe is not None:
            configuration = universe.contiguousObjectConfiguration(
                               [self], configuration)
        file = ConfigIO.OutputFile(filename, format)
        file.write(self, configuration)
        file.close()

    def view(self, configuration = None, format = 'pdb'):
        """
        Start an external viewer for the object in the given
        configuration.

        @param configuration: the configuration to be visualized
        @type configuration: L{MMTK.Configuration}
        @param format: 'pdb' (for running $PDBVIEWER) or 'vrml'
                       (for running $VRMLVIEWER). An optional
                       subformat specification can be added, see
                       L{writeToFile} for the details.
        """
        universe = self.universe()
        if universe is not None:
            configuration = universe.contiguousObjectConfiguration([self],
                                                                configuration)
        Visualization.viewConfiguration(self, configuration, format)

    def kineticEnergy(self, velocities = None):
        """
        @param velocities: a set of velocities for all atoms, or
                           C{None} for the current velocities
        @type velocities: L{MMTK.ParticleVector}
        @returns: the kinetic energy
        @rtype: C{float}
        """
        if velocities is None:
            velocities = self.atomList()[0].universe().velocities()
        energy = 0.
        for b in self.beadIterator():
            v = velocities[b]
            energy += energy + b._mass*(v*v)
        return 0.5*energy

    def temperature(self, velocities = None):
        """
        @param velocities: a set of velocities for all atoms, or
                           C{None} for the current velocities
        @type velocities: L{MMTK.ParticleVector}
        @returns: the temperature
        @rtype: C{float}
        """
        energy = self.kineticEnergy(velocities)
        return 2.*energy/(self.degreesOfFreedom()*Units.k_B)

    def momentum(self, velocities = None):
        """
        @param velocities: a set of velocities for all atoms, or
                           C{None} for the current velocities
        @type velocities: L{MMTK.ParticleVector}
        @returns: the momentum
        @rtype: C{Scientific.Geometry.Vector}
        """
        if velocities is None:
            velocities = self.atomList()[0].universe().velocities()
        return sum((b._mass*velocities[b] for b in self.beadIterator()),
                   Vector(0., 0., 0.))

    def angularMomentum(self, velocities = None, conf = None):
        """
        @param velocities: a set of velocities for all atoms, or
                           C{None} for the current velocities
        @type velocities: L{MMTK.ParticleVector}
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration}
        @returns: the angluar momentum
        @rtype: C{Scientific.Geometry.Vector}
        """
        if velocities is None:
            velocities = self.atomList()[0].universe().velocities()
        cm = self.centerOfMass(conf)
        return sum((b._mass*b.position(conf).cross(velocities[b])
                    for b in self.beadIterator()),
                   Vector(0., 0., 0.))

    def angularVelocity(self, velocities = None, conf = None):
        """
        @param velocities: a set of velocities for all atoms, or
                           C{None} for the current velocities
        @type velocities: L{MMTK.ParticleVector}
        @param conf: a configuration object, or C{None} for the
                     current configuration
        @type conf: L{MMTK.Configuration}
        @returns: the angluar velocity
        @rtype: C{Scientific.Geometry.Vector}
        """
        if velocities is None:
            velocities = self.atomList()[0].universe().velocities()
        cm, inertia = self.centerAndMomentOfInertia(conf)
        l = sum((b._mass*b.position(conf).cross(velocities[b])
                 for b in self.beadIterator()),
                Vector(0., 0., 0.))
        return inertia.inverse()*l
        
    def universe(self):
        """
        @returns: the universe of which the object is part. For an
                  object that is not part of a universe, the result is
                  C{None}
        """
        atoms = self.atomList()
        if not atoms:
            return None
        universe = atoms[0].universe()
        for a in atoms[1:]:
            if a.universe() is not universe:
                return None
        return universe

    def charge(self):
        """
        @returns: the total charge of the object. This is defined only
                  for objects that are part of a universe with a force
                  field that defines charges.
        @rtype: C{float}
        """
        return self.universe().forcefield().charge(self)

    def dipole(self, reference = None):
        """
        @returns: the total dipole moment of the object. This is defined only
                  for objects that are part of a universe with a force field
                  that defines charges.
        @rtype: C{Scientific.Geometry.Vector}
        """
        return self.universe().forcefield().dipole(self, reference)

    def booleanMask(self):
        """
        @returns: a ParticleScalar object that contains a value of 1
                  for each atom that is in the object and a value of 0 for all
                  other atoms in the universe
        @rtype: L{MMTK.ParticleScalar}
        """
        universe = self.universe()
        if universe is None:
            raise ValueError("object not in a universe")
        array = N.zeros((universe.numberOfAtoms(),), N.Int)
        mask = ParticleProperties.ParticleScalar(universe, array)
        for a in self.atomIterator():
            mask[a] = 1
        return mask

#
# This class defines a general collection that can contain
# chemical objects and other collections.
#
class Collection(GroupOfAtoms, Visualization.Viewable):

    """
    Collection of chemical objects

    Collections permit the grouping of arbitrary chemical objects
    (atoms, molecules, etc.) into one object for the purpose of analysis
    or manipulation.

    Collections permit length inquiry, item extraction by indexing,
    and iteration, like any Python sequence object. Two collections
    can be added to yield a collection that contains the combined
    elements.
    """

    def __init__(self, *objects):
        """
        @param objects: a chemical object or a sequence of chemical objects that
                        define the initial content of the collection.
        """
        self.objects = []
        self.addObject(objects)

    is_collection = 1

    def addObject(self, object):
        """
        Add objects to the collection.

        @param object: the object(s) to be added. If it is another collection
                       or a list, all of its elements are added
        """
        from MMTK.ChemicalObjects import isChemicalObject
        if isChemicalObject(object):
            self.addChemicalObject(object)
        elif isCollection(object):
            self.addChemicalObjectList(object.objectList())
        elif Utility.isSequenceObject(object):
            if object and isChemicalObject(object[0]):
                self.addChemicalObjectList(list(object))
            else:
                for o in object:
                    self.addObject(o)
        else:
            raise TypeError('Wrong object type in collection')

    def addChemicalObject(self, object):
        self.objects.append(object)

    def addChemicalObjectList(self, list):
        self.objects.extend(list)

    def removeObject(self, object):
        """
        Remove an object or a list or collection of objects from the
        collection. The object(s) to be removed must be elements of the
        collection.

        @param object: the object to be removed, or a list or collection
                       of objects whose elements are to be removed
        @raises ValueError: if the object is not an element of the collection
        """
        from MMTK.ChemicalObjects import isChemicalObject
        if isChemicalObject(object):
            self.removeChemicalObject(object)
        elif isCollection(object) or Utility.isSequenceObject(object):
            for o in object:
                self.removeObject(o)
        else:
            raise ValueError('Object not in this collection')

    def removeChemicalObject(self, object):
        try:
            self.objects.remove(object)
        except ValueError:
            raise ValueError('Object not in this collection')

    def selectShell(self, point, r1, r2=0.):
        """
        Select objects in a spherical shell around a central point.

        @param point: the center of the spherical shell
        @type point: C{Scientific.Geometry.Vector}
        @param r1: inner or outer radius of the shell
        @type r1: C{float}
        @param r2: inner or outer radius of the shell (default: 0.)
        @type r2: C{float}
        @returns: a collection of all elements whose
                  distance from point is between r1 and r2
        @rtype: L{Collection}
        """
        if r1 > r2:
            r1, r2 = r2, r1
        universe = self.universe()
        in_shell = []
        for o in self.objects:
            for a in o.atomIterator():
                r =  universe.distance(a.position(), point)
                if r >= r1 and r <= r2:
                    in_shell.append(o)
                    break
        return Collection(in_shell)

    def selectBox(self, p1, p2):
        """
        Select objects in a rectangular volume

        @param p1: one corner of the rectangular volume
        @type p1: C{Scientific.Geometry.Vector}
        @param p2: the other corner of the rectangular volume
        @type p2: C{Scientific.Geometry.Vector}
        @returns: a collection of all elements that lie
                  within the rectangular volume
        @rtype: L{Collection}
        """
        x1 = N.minimum(p1.array, p2.array)
        x2 = N.maximum(p1.array, p2.array)
        in_box = []
        for o in self.objects:
            r = o.position().array
            if N.logical_and.reduce( \
                N.logical_and(N.less_equal(x1, r),
                              N.less(r, x2))):
                in_box.append(o)
        return Collection(in_box)

    def objectList(self, type = None):
        """
        Make a list of all objects in the collection that are instances
        of a specific type or of one of its subtypes.

        @param type: the type that serves as a filter. If C{None},
                     all objects are returned
        @returns: the objects that match the given type
        @rtype: C{list}
        """
        if type is None:
            return self.objects
        else:
            return [o for o in self.objects if isinstance(o, type)]

    def atomList(self):
        """
        @returns: a list containing all atoms of all objects in the collection
        @rtype: C{list}
        """
        atoms = []
        for o in self.objectList():
            atoms.extend(o.atomList())
        return atoms

    def atomIterator(self):
        return itertools.chain(*(o.atomIterator() for o in self.objects))
    
    def numberOfAtoms(self):
        """
        @returns: the total number of atoms in the objects of the collection
        @rtype: C{int}
        """
        return sum(o.numberOfAtoms() for o in self.objectList())
    
    def universe(self):
        """
        @returns: the universe of which all objects in the collection
                  are part. If no such universe exists, the return value
                  is C{None}
        """
        if not self.objects:
            return None
        universe = self.objects[0].universe()
        for o in self.objects[1:]:
            if o.universe() is not universe:
                return None
        return universe

    def __len__(self):
        """
        @returns: the number of objects in the collection
        @rtype: C{int}
        """
        return len(self.objects)

    def __getitem__(self, item):
        """
        @param item: an index into the object list
        @type item: C{int}
        @returns: the object with the given index
        """
        return self.objects[item]

    def __iter__(self):
        return self.objects.__iter__()

    def __add__(self, other):
        return Collection(self.objectList(), other.objectList())

    def __str__(self):
        return "Collection of %d objects" % len(self.objects)

    def map(self, function):
        """
        Apply a function to all objects in the collection and
        return the list of the results. If the results are chemical
        objects, a Collection object is returned instead of a list.

        @param function: the function to be applied
        @type function: callable
        @returns: the list or collection of the results
        """
        from MMTK.ChemicalObjects import isChemicalObject
        list = [function(o) for o in self.objectList()]
        if list and isChemicalObject(list[0]):
            return Collection(list)
        else:
            return list

    def bondedUnits(self):
        bu = []
        for o in self.objects:
            bu = bu + o.bondedUnits()
        return bu

    def degreesOfFreedom(self):
        return GroupOfAtoms.degreesOfFreedom(self) \
               - self.numberOfDistanceConstraints()

    def distanceConstraintList(self):
        """
        @returns: the list of distance constraints
        @rtype: C{list}
        """
        dc = []
        for o in self.objects:
            dc.extend(o.distanceConstraintList())
        return dc

    def numberOfDistanceConstraints(self):
        """
        @returns: the number of distance constraints
        """
        return sum(o.numberOfDistanceConstraints() for o in self.objects)

    def setBondConstraints(self, universe=None):
        """
        Set distance constraints for all bonds
        """
        if universe is None:
            universe = self.universe()
        for o in self.objects:
            o.setBondConstraints(universe)

    def removeDistanceConstraints(self, universe=None):
        """
        Remove all distance constraints
        """
        if universe is None:
            universe = self.universe()
        for o in self.objects:
            o.removeDistanceConstraints(universe)

    def _graphics(self, conf, distance_fn, model, module, options):
        lists = []
        for o in self.objects:
            lists.append(o._graphics(conf, distance_fn, model,
                                     module, options))
        return sum(lists)

    def __copy__(self):
        return self.__class__(copy.copy(self.objects))


# type check for collections

def isCollection(object):
    """
    @param object: any Python object
    @returns: C{True} if the object is a L{Collection}
    """
    return hasattr(object, 'is_collection')

#
# This class defines a partitioned collection. Such collections
# divide their objects into cubic boxes according to their positions.
# It is then possible to find potential neighbours much more efficiently.
#
class PartitionedCollection(Collection):

    """
    Collection with cubic partitions

    A PartitionedCollection differs from a plain Collection by
    sorting its elements into small cubic cells. This makes adding
    objects slower, but geometrical operations like 
    selectShell become much faster for a large number of
    objects.
    """

    def __init__(self, partition_size, *objects):
        """
        @param partition_size: the edge length of the cubic cells
        @param objects: a chemical object or a sequence of chemical objects that
                        define the initial content of the collection.
        """
        self.partition_size = 1.*partition_size
        self.undefined = []
        self.partition = {}
        self.addObject(objects)

    def addChemicalObject(self, object):
        p = object.position()
        if p is None:
            self.undefined.append(object)
        else:
            index = self.partitionIndex(p)
            try:
                partition = self.partition[index]
            except KeyError:
                partition = []
                self.partition[index] = partition
            partition.append(object)
        self.all = None

    def addChemicalObjectList(self, list):
        for object in list:
            self.addChemicalObject(object)

    def removeChemicalObject(self, object):
        p = object.position()
        if p is None:
            self.undefined.remove(object)
        else:
            index = self.partitionIndex(p)
            try:
                partition = self.partition[index]
            except KeyError:
                raise ValueError('Object not in this collection')
            try:
                partition.remove(object)
            except ValueError:
                raise ValueError('Object not in this collection')
        self.all = None

    def partitionIndex(self, x):
        return (int(N.floor(x[0]/self.partition_size)),
                int(N.floor(x[1]/self.partition_size)),
                int(N.floor(x[2]/self.partition_size)))

    def objectList(self):
        return sum(self.partition.values(), [self.undefined])

    def __len__(self):
        return sum(len(p) for p in self.partition.values()) + \
               len(self.undefined)

    def __getitem__(self, item):
        if self.all is None:
            self.all = self.objectList()
        if item >= len(self.all):
            self.all = None
            raise IndexError
        return self.all[item]

    def __copy__(self):
        return self.__class__(self.partition_size,
                              copy.copy(self.objectList()))

    def partitions(self):
        """
        @returns: a list of cubic partitions. Each partition is specified
                  by a tuple containing two vectors (describing the diagonally
                  opposite corners) and the list of objects in the partition.
        """
        list = []
        for index, objects in self.partition.items():
            min = Vector(index)*self.partition_size
            max = min + Vector(3*[self.partition_size])
            list.append((min, max, objects))
        return list

    def selectCube(self, point, edge):
        x = int(round(point[0]/self.partition_size))
        y = int(round(point[1]/self.partition_size))
        z = int(round(point[2]/self.partition_size))
        d = (Vector(x, y, z)*self.partition_size-point).length()
        n = int(N.ceil((edge + d)/(2.*self.partition_size)))
        objects = []
        for nx in range(-n, n):
            for ny in range(-n, n):
                for nz in range(-n, n):
                    try:
                        objects.append(self.partition[(nx+x, ny+y, nz+z)])
                    except KeyError: pass
        return Collection(objects)

    def selectShell(self, point, min, max=0):
        if min > max:
            min, max = max, min
        objects = Collection()
        minsq = min**2
        maxsq = max**2
        for index in self.partition.keys():
            d1 = self.partition_size*N.array(index) - point.array
            d2 = d1 + self.partition_size
            dmin = (d1 > 0.)*d1 - (d2 < 0.)*d2
            dminsq = N.add.reduce(dmin**2)
            dmaxsq = N.add.reduce(N.maximum(d1**2, d2**2))
            if dminsq >= minsq and dmaxsq <= maxsq:
                objects.addObject(self.partition[index])
            elif dmaxsq >= minsq and dminsq <= maxsq:
                o = Collection(self.partition[index]).selectShell(point,
                                                                  min, max)
                objects.addObject(o)
        return objects

    def pairsWithinCutoff(self, cutoff):
        """
        @param cutoff: a cutoff for pair distances
        @returns: a list containing all pairs of objects in the
                  collection whose center-of-mass distance is less than
                  the cutoff
        @rtype: C{list}
        """
        pairs = []
        positions = {}
        for index, objects in self.partition.items():
            pos = map(lambda o: o.position(), objects)
            positions[index] = pos
            for o1, o2 in Utility.pairs(zip(objects, pos)):
                if (o2[1]-o1[1]).length() <= cutoff:
                    pairs.append((o1[0], o2[0]))
        partition_cutoff = int(N.floor((cutoff/self.partition_size)**2))
        ones = N.array([1,1,1])
        zeros = N.array([0,0,0])
        keys = self.partition.keys()
        for i in range(len(keys)):
            p1 = keys[i]
            for j in range(i+1, len(keys)):
                p2 = keys[j]
                d = N.maximum(abs(N.array(p2)-N.array(p1)) -
                              ones, zeros)
                if N.add.reduce(d*d) <= partition_cutoff:
                    for o1, pos1 in zip(self.partition[p1],
                                        positions[p1]):
                        for o2, pos2 in zip(self.partition[p2],
                                            positions[p2]):
                            if (pos2-pos1).length() <= cutoff:
                                pairs.append((o1, o2))
        return pairs

#
# A special form of partitioned collection that stores the atoms
# of all objects that are added to it.
#
class PartitionedAtomCollection(PartitionedCollection):

    """
    Partitioned collection of atoms

    PartitionedAtomCollection objects behave like PartitionedCollection
    atoms, except that they store only atoms. When a composite chemical
    object is added, its atoms are stored instead.
    """

    def __init__(*args):
        apply(PartitionedCollection.__init__, args)

    def addChemicalObject(self, object):
        for atom in object.atomIterator():
            PartitionedCollection.addChemicalObject(self, atom)

    def removeChemicalObject(self, object):
        for atom in object.atomIterator():
            PartitionedCollection.removeChemicalObject(self, atom)
