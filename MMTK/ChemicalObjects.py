# This module implements classes that represent atoms, molecules, and
# complexes. They are made as copies from blueprints in the database.
#
# Written by Konrad Hinsen
#

"""
Atoms, groups, molecules and similar objects
"""

__docformat__ = 'restructuredtext'

from MMTK import Bonds, Collections, Database, \
                 Units, Utility, Visualization
from Scientific.Geometry import Vector
from Scientific.Geometry import Objects3D
from Scientific import N
import copy

#
# The base class for all chemical objects.
#
class ChemicalObject(Collections.GroupOfAtoms, Visualization.Viewable):

    """
    General chemical object

    This is an abstract base class that implements methods which
    are applicable to any chemical object (atom, molecule, etc.).
    """

    def __init__(self, blueprint, memo):
        if isinstance(blueprint, basestring):
            blueprint = self.blueprintclass(blueprint)
        self.type = blueprint.type
        if hasattr(blueprint, 'name'):
            self.name = blueprint.name
        if memo is None: memo = {}
        memo[id(blueprint)] = self
        for attr in blueprint.instance:
            setattr(self, attr,
                    Database.instantiate(getattr(blueprint, attr), memo))

    is_chemical_object = True
    is_incomplete = False
    is_modified = False

    __safe_for_unpickling__ = True
    __had_initargs__ = True

    def __hash__(self):
        return id(self)

    def __getattr__(self, attr):
        if attr[:1] == '_' or attr[:3] == 'is_':
            raise AttributeError
        else:
            return getattr(self.type, attr)

    def isSubsetModel(self):
        return False

    def addProperties(self, properties):
        if properties:
            for item in properties.items():
                if hasattr(self, item[0]) and item[0] != 'name':
                    raise TypeError('attribute '+item[0]+' already defined')
                setattr(self, item[0], item[1])

    def binaryProperty(self, properties, name, default):
        value = default
        try:
            value = properties[name]
            del properties[name]
        except KeyError:
            pass
        return value

    def topLevelChemicalObject(self):
        """Returns the highest-level chemical object of which
        the current object is a part."""
        if self.parent is None or not isChemicalObject(self.parent):
            return self
        else:
            return self.parent.topLevelChemicalObject()

    def universe(self):
        """
        :returns: the universe to which the object belongs,
                  or None if the object does not belong to any universe
        :rtype: :class:`~MMTK.Universe.Universe`
        """
        if self.parent is None:
            return None
        else:
            return self.parent.universe()

    def bondedUnits(self):
        """
        :returns: the largest subobjects which can contain bonds.
                  There are no bonds between any of the subobjects
                  in the list.
        :rtype: list
        """
        return [self]

    def atomList(self):
        """
        :returns: a list containing all atoms in the object
        :rtype: list
        """
        pass

    def atomIterator(self):
        """
        :returns: an iterator over all atoms in the object
        :rtype: iterator
        """
        pass

    def fullName(self):
        """
        :returns: the full name of the object. The full name consists
                  of the proper name of the object preceded by
                  the full name of its parent separated by a dot.
        :rtype: str
        """
        if self.parent is None or not isChemicalObject(self.parent):
            return self.name
        else:
            return self.parent.fullName() + '.' + self.name

    def degreesOfFreedom(self):
        """
        :returns: the number of degrees of freedom of the object
        :rtype: int
        """
        return Collections.GroupOfAtoms.degreesOfFreedom(self) \
               - self.numberOfDistanceConstraints()

    def distanceConstraintList(self):
        """
        :returns: the distance constraints of the object
        :rtype: list
        """
        return []

    def _distanceConstraintList(self):
        return []

    def traverseBondTree(self, function = None):
        return []

    def numberOfDistanceConstraints(self):
        """
        :returns: the number of distance constraints of the object
        :rtype: int
        """
        return 0

    def setBondConstraints(self, universe=None):
        """
        Sets distance constraints for all bonds.
        """
        pass

    def removeDistanceConstraints(self, universe=None):
        """
        Removes all distance constraints.
        """
        pass

    def setRigidBodyConstraints(self, universe = None):
        """
        Sets distance constraints that make the object fully rigid.
        """
        if universe is None:
            universe = self.universe()
        if universe is None:
            import Universe
            universe = Universe.InfiniteUniverse()
        atoms = self.atomList()
        if len(atoms) > 1:
            self.addDistanceConstraint(atoms[0], atoms[1],
                                       universe.distance(atoms[0], atoms[1]))
        if len(atoms) > 2:
            self.addDistanceConstraint(atoms[0], atoms[2],
                                       universe.distance(atoms[0], atoms[2]))
            self.addDistanceConstraint(atoms[1], atoms[2],
                                       universe.distance(atoms[1], atoms[2]))
        if len(atoms) > 3:
            for a in atoms[3:]:
                self.addDistanceConstraint(atoms[0], a,
                                           universe.distance(atoms[0], a))
                self.addDistanceConstraint(atoms[1], a,
                                           universe.distance(atoms[1], a))
                self.addDistanceConstraint(atoms[2], a,
                                           universe.distance(atoms[2], a))

    def getAtomProperty(self, atom, property):
        """
        Retrieve atom properties from the chemical database.

        Note: the property is first looked up in the database entry
        for the object on which the method is called. If the lookup
        fails, the complete hierarchy from the atom to the top-level
        object is constructed and traversed starting from the top-level
        object until the property is found. This permits database entries
        for higher-level objects to override property definitions in
        its constituents.

        At the atom level, the property is retrieved from an attribute
        with the same name. This means that properties at the atom
        level can be defined both in the chemical database and for
        each atom individually by assignment to the attribute.
        
        :returns: the value of the specified property for the given
                  atom from the chemical database.
        """

    def writeXML(self, file, memo, toplevel=1):
        if self.type is None:
            name = 'm' + `memo['counter']`
            memo['counter'] = memo['counter'] + 1
            memo[id(self)] = name
            atoms = copy.copy(self.atoms)
            bonds = copy.copy(self.bonds)
            for group in self.groups:
                group.writeXML(file, memo, 0)
                for atom in group.atoms:
                    atoms.remove(atom)
                for bond in group.bonds:
                    bonds.remove(bond)
            file.write('<molecule id="%s">\n' % name)
            for group in self.groups:
                file.write('  <molecule ref="%s" title="%s"/>\n'
                           % (memo[id(group.type)], group.name))
            if atoms:
                file.write('  <atomArray>\n')
                for atom in atoms:
                    file.write('    <atom title="%s" elementType="%s"/>\n'
                               % (atom.name, atom.type.symbol))
                file.write('  </atomArray>\n')
            if bonds:
                file.write('  <bondArray>\n')
                for bond in bonds:
                    a1n = self.relativeName(bond.a1)
                    a2n = self.relativeName(bond.a2)
                    file.write('    <bond atomRefs2="%s %s"/>\n' % (a1n, a2n))
                file.write('  </bondArray>\n')
            file.write('</molecule>\n')
        else:
            name = self.name
            self.type.writeXML(file, memo)
            atom_names = self.type.getXMLAtomOrder()
        if toplevel:
            return ['<molecule ref="%s"/>' % name]
        else:
            return None

    def getXMLAtomOrder(self):
        if self.type is None:
            atoms = []
            for group in self.groups:
                atoms.extend(group.getXMLAtomOrder())
            for atom in self.atoms:
                if atom not in atoms:
                    atoms.append(atom)
        else:
            atom_names = self.type.getXMLAtomOrder()
            atoms = []
            for name in atom_names:
                parts = name.split(':')
                object = self
                for attr in parts:
                    object = getattr(object, attr)
                atoms.append(object)
        return atoms

    def relativeName(self, object):
        name = ''
        while object is not None and object != self:
            name = object.name + ':' + name
            object = object.parent
        return name[:-1]

    def relativeName_unused(self, object):
        name = ''
        while object is not None and object != self:
            for attr, value in object.parent.__dict__.items():
                if value is object:
                    name = attr + ':' + name
                    break
            object = object.parent
        return name[:-1]

    def description(self, index_map = None):
        tag = Utility.uniqueAttribute()
        s = self._description(tag, index_map, 1)
        for a in self.atomList():
            delattr(a, tag)
        return s

    def __repr__(self):
        return self.__class__.__name__ + ' ' + self.fullName()
    __str__ = __repr__

    def __copy__(self):
        return copy.deepcopy(self, {id(self.parent): None})

# Type check

def isChemicalObject(object):
    """
    :returns: True if object is a chemical object
    """
    return hasattr(object, 'is_chemical_object')

#
# The second base class for all composite chemical objects.
#
class CompositeChemicalObject(object):

    """
    Chemical object containing subobjects

    This is an abstract base class that implements methods
    which can be used with any composite chemical object,
    i.e. any chemical object that is not an atom.
    """

    def __init__(self, properties):
        if properties.has_key('configuration'):
            conf = properties['configuration']
            self.configurations[conf].applyTo(self)
            del properties['configuration']
        elif hasattr(self, 'configurations') and \
             self.configurations.has_key('default'):
            self.configurations['default'].applyTo(self)
        if properties.has_key('position'):
            self.translateTo(properties['position'])
            del properties['position']
        self.addProperties(properties)

    def atomList(self):
        return self.atoms

    def atomIterator(self):
        return iter(self.atoms)

    def setPosition(self, atom, position):
        if atom.__class__ is Database.AtomReference:
            atom = self.atoms[atom.number]
        atom.setPosition(position)

    def setIndex(self, atom, index):
        if atom.__class__ is Database.AtomReference:
            atom = self.atoms[atom.number]
        atom.setIndex(index)

    def getAtom(self, atom):
        if atom.__class__ is Database.AtomReference:
            atom = self.atoms[atom.number]
        return atom

    def getReference(self, atom):
        if atom.__class__ is Database.AtomReference:
            return atom
        return Database.AtomReference(self.atoms.index(atom))

    def getAtomProperty(self, atom, property, levels = None):
        try:
            return atom.__dict__[property]
        except KeyError:
            try:
                return getattr(self, property)[self.getReference(atom)]
            except (AttributeError, KeyError):
                if levels is None:
                    object = atom
                    levels = []
                    while object != self:
                        levels.append(object)
                        object = object.parent
                if not levels:
                    raise KeyError('Property ' + property +
                                    ' not defined for  ', `atom`)
                return levels[-1].getAtomProperty(atom, property, levels[:-1])

    def deleteUndefinedAtoms(self):
        delete = [a for a in self.atoms if a.position() is None]
        for a in delete:
            a.delete()

    def _deleteAtom(self, atom):
        self.atoms.remove(atom)
        self.is_modified = True
        self.type = None
        if self.parent is not None:
            self.parent._deleteAtom(atom)

    def distanceConstraintList(self):
        dc = self._distanceConstraintList()
        for o in self._subunits():
            dc = dc + o._distanceConstraintList()
        return dc

    def _distanceConstraintList(self):
        try:
            return self.distance_constraints
        except AttributeError:
            return []

    def numberOfDistanceConstraints(self):
        n = len(self._distanceConstraintList()) + \
            sum(len(o._distanceConstraintList()) for o in self._subunits())
        return n

    def setBondConstraints(self, universe=None):
        if universe is None:
            universe = self.universe()
        bond_database = universe.bondLengthDatabase()
        for o in self.bondedUnits():
            o._setBondConstraints(universe, bond_database)

    def _setBondConstraints(self, universe, bond_database):
        self.distance_constraints = []
        for bond in self.bonds:
            d = bond_database.bondLength(bond)
            if d is None:
                d = universe.distance(bond.a1, bond.a2)
            self.addDistanceConstraint(bond.a1, bond.a2, d)

    def addDistanceConstraint(self, atom1, atom2, distance):
        try:
            self.distance_constraints.append((atom1, atom2, distance))
        except AttributeError:
            self.distance_constraints = [(atom1, atom2, distance)]

    def removeDistanceConstraints(self, universe=None):
        try:
            del self.distance_constraints
        except AttributeError:
            pass
        for o in self._subunits():
            o.removeDistanceConstraints()

    def traverseBondTree(self, function = None):
        self.setBondAttributes()
        todo = [self.atoms[0]]
        done = {todo[0]: True}
        bonds = []
        while todo:
            next_todo = []
            for atom in todo:
                bonded = atom.bondedTo()
                for other in bonded:
                    if not done.get(other, False):
                        if function is None:
                            bonds.append((atom, other))
                        else:
                            bonds.append((function(atom), function(other)))
                        next_todo.append(other)
                        done[other] = True
            todo = next_todo
        self.clearBondAttributes()
        return bonds

    def _description(self, tag, index_map, toplevel):
        letter, kwargs = self._descriptionSpec()
        s = [letter, '(', `self.name`, ',[']
        s.extend([o._description(tag, index_map, 0) + ','
                  for o in self._subunits()])
        s.extend([a._description(tag, index_map, 0) + ','
                  for a in self.atoms if not hasattr(a, tag)])
        s.append(']')
        if toplevel:
            s.extend([',', `self._typeName()`])
        if kwargs is not None:
            s.extend([',', kwargs])
        constraints = self._distanceConstraintList()
        if constraints:
            s.append(',dc=[')
            if index_map is None:
                s.extend(['(%d,%d,%f),' % (c[0].index, c[1].index, c[2])
                          for c in constraints])
            else:
                s.extend(['(%d,%d,%f),' % (index_map[c[0].index],
                                           index_map[c[1].index], c[2])
                          for c in constraints])
            s.append(']')
        s.append(')')
        return ''.join(s)

    def _typeName(self):
        return self.type.name

    def _graphics(self, conf, distance_fn, model, module, options):
        glist = []
        for bu in self.bondedUnits():
            for a in bu.atomList():
                glist.extend(a._graphics(conf, distance_fn, model,
                                         module, options))
            if hasattr(bu, 'bonds'):
                for b in bu.bonds:
                    glist.extend(b._graphics(conf, distance_fn, model,
                                             module, options))
        return glist

#
# The classes for atoms, groups, molecules, and complexes.
#
class Atom(ChemicalObject):

    """
    Atom
    """

    def __init__(self, atom_spec, _memo = None, **properties):
        """
        :param atom_spec: a string (not case sensitive) specifying
                          the chemical element
        :type atom_spec: str
        :keyword position: the position of the atom
        :type position: Scientific.Geometry.Vector
        :keyword name: a name given to the atom
        :type name: str
        """
        Utility.uniqueID.registerObject(self)
        ChemicalObject.__init__(self, atom_spec, _memo)
        self._mass = self.type.average_mass
        self.array = None
        self.index = None
        if properties.has_key('position'):
            self.setPosition(properties['position'])
            del properties['position']
        self.addProperties(properties)

    blueprintclass = Database.BlueprintAtom

    def __getstate__(self):
        state = copy.copy(self.__dict__)
        if self.array is not None:
            state['array'] = None
            state['pos'] = Vector(self.array[self.index,:])
        return state

    def atomList(self):
        return [self]

    def atomIterator(self):
        yield self

    def setPosition(self, position):
        """
        Changes the position of the atom.

        :param position: the new position
        :type position: Scientific.Geometry.Vector
        """
        if position is None:
            if self.array is None:
                try: del self.pos
                except AttributeError: pass
            else:
                self.array[self.index,0] = Utility.undefined
                self.array[self.index,1] = Utility.undefined
                self.array[self.index,2] = Utility.undefined
        else:
            if self.array is None:
                self.pos = position
            else:
                self.array[self.index,0] = position[0]
                self.array[self.index,1] = position[1]
                self.array[self.index,2] = position[2]

    translateTo = setPosition

    def position(self, conf = None):
        """
        :returns: the position in configuration conf. If conf is 
                  'None', use the current configuration. If the atom has
                  not been assigned a position, the return value is None.
        """
        if conf is None:
            if self.array is None:
                try:
                    return self.pos
                except AttributeError:
                    return None
            else:
                if N.logical_or.reduce(
                    N.greater(self.array[self.index, :],
                              Utility.undefined_limit)):
                    return None
                else:
                    return Vector(self.array[self.index,:])
        else:
            return conf[self]

    centerOfMass = position

    def setMass(self, mass):
        """
        Changes the mass of the atom.

        :param mass: the mass
        :type mass: float
        """
        self._mass = mass
        universe = self.universe()
        if universe is not None:
            universe._changed(True)

    def getAtom(self, atom):
        return self

    def translateBy(self, vector):
        if self.array is None:
            self.pos = self.pos + vector
        else:
            self.array[self.index,0] = self.array[self.index,0] + vector[0]
            self.array[self.index,1] = self.array[self.index,1] + vector[1]
            self.array[self.index,2] = self.array[self.index,2] + vector[2]

    def numberOfPoints(self):
        return 1

    numberOfCartesianCoordinates = numberOfPoints

    def setIndex(self, index):
        if self.index is not None and self.index != index:
            raise ValueError('Wrong atom index')
        self.index = index

    def setArray(self, index):
        self.index = index
        self.array = None

    def getArray(self):
        return self.array

    def unsetArray(self):
        self.pos = self.position()
        self.array = None

    def setBondAttribute(self, atom):
        try:
            self.bonded_to__.append(atom)
        except AttributeError:
            self.bonded_to__ = [atom]

    def clearBondAttribute(self):
        try:
            del self.bonded_to__
        except AttributeError:
            pass

    def bondedTo(self):
        """
        :returns: a list of all atoms to which a chemical bond exists.
        :rtype: list
        """
        try:
            return self.bonded_to__
        except AttributeError:
            if self.parent is None or not isChemicalObject(self.parent):
                return []
            else:
                return self.parent.bondedTo(self)

    def delete(self):
        if self.parent is not None:
            self.parent._deleteAtom(self)           

    def getAtomProperty(self, atom, property, levels = None):
        if self != atom:
            raise ValueError("Wrong atom")
        return getattr(self, property)

    def _description(self, tag, index_map, toplevel):
        setattr(self, tag, None)
        if index_map is None:
            index = self.index
        else:
            index = index_map[self.index]
        if toplevel:
            return 'A(' + `self.name` + ',' + `index` + ',' + \
                   `self.symbol` + ')'
        else:
            return 'A(' + `self.name` + ',' + `index` + ')'

    def _graphics(self, conf, distance_fn, model, module, options):
        #PJC change:
        if model == 'ball_and_stick':
            color = self._atomColor(self, options)
            material = module.DiffuseMaterial(color)
            radius = options.get('ball_radius', 0.03)
            return [module.Sphere(self.position(), radius, material=material)]
        elif model == 'vdw' or model == 'vdw_and_stick':
            color = self._atomColor(self, options)
            material = module.DiffuseMaterial(color)
            try:
                radius = self.vdW_radius
            except:
                radius = options.get('ball_radius', 0.03)
            return [module.Sphere(self.position(), radius, material=material)]
        else:
            return []

        if model != 'ball_and_stick':
            return []
        color = self._atomColor(self, options)
        material = module.DiffuseMaterial(color)
        radius = options.get('ball_radius', 0.03)
        return [module.Sphere(self.position(), radius, material=material)]

    def writeXML(self, file, memo, toplevel=1):
        return ['<atom title="%s" elementType="%s"/>'
                % (self.name, self.type.symbol)]

    def getXMLAtomOrder(self):
        return [self]

class Group(CompositeChemicalObject, ChemicalObject):

    """
    Group of bonded atoms

    Groups can contain atoms and other groups, and link them by chemical
    bonds. They are used to represent functional groups or any other
    part of a molecule that has a well-defined identity.

    Groups cannot be created in application programs, but only in
    database definitions for molecules or through a MoleculeFactory.
    """

    def __init__(self, group_spec, _memo = None, **properties):
        """
        :param group_spec: a string (not case sensitive) that specifies
                           the group name in the chemical database
        :type group_spec: str
        :keyword position: the position of the center of mass of the group
        :type position: Scientific.Geometry.Vector
        :keyword name: a name given to the group
        :type name: str
        """
        if group_spec is not None:
            # group_spec is None when called from MoleculeFactory
            ChemicalObject.__init__(self, group_spec, _memo)
            self.addProperties(properties)

    blueprintclass = Database.BlueprintGroup
    is_incomplete = True

    def bondedTo(self, atom):
        if self.parent is None or not isChemicalObject(self.parent):
            return []
        else:
            return self.parent.bondedTo(atom)

    def setBondAttributes(self):
        pass

    def clearBondAttributes(self):
        pass

    def _subunits(self):
        return self.groups

    def _descriptionSpec(self):
        return "G", None

class Molecule(CompositeChemicalObject, ChemicalObject):

    """Molecule

    Molecules consist of atoms and groups linked by bonds.
    """

    def __init__(self, molecule_spec, _memo = None, **properties):
        """
        :param molecule_spec: a string (not case sensitive) that specifies
                              the molecule name in the chemical database
        :type molecule_spec: str
        :keyword position: the position of the center of mass of the molecule
        :type position: Scientific.Geometry.Vector
        :keyword name: a name given to the molecule
        :type name: str
        :keyword configuration: the name of a configuration listed in the
                                database definition of the molecule, which
                                is used to initialize the atom positions.
                                If no configuration is specified, the
                                configuration named "default" will be used,
                                if it exists. Otherwise the atom positions
                                are undefined.
        :type configuration: str
        """
        if molecule_spec is not None:
            # molecule_spec is None when called from MoleculeFactory
            ChemicalObject.__init__(self, molecule_spec, _memo)
            properties = copy.copy(properties)
            CompositeChemicalObject.__init__(self, properties)
            self.bonds = Bonds.BondList(self.bonds)

    blueprintclass = Database.BlueprintMolecule

    def bondedTo(self, atom):
        return self.bonds.bondedTo(atom)

    def setBondAttributes(self):
        self.bonds.setBondAttributes()

    def clearBondAttributes(self):
        for a in self.atoms:
            a.clearBondAttribute()

    def _subunits(self):
        return self.groups

    def _descriptionSpec(self):
        return "M", None

    def addGroup(self, group, bond_atom_pairs):
        for a1, a2 in bond_atom_pairs:
            o1 = a1.topLevelChemicalObject()
            o2 = a2.topLevelChemicalObject()
            if set([o1, o2]) != set([self, group]):
                raise ValueError("bond %s-%s outside object" %
                                  (str(a1), str(a2)))
        self.groups.append(group)
        self.atoms = self.atoms + group.atoms
        group.parent = self
        self.clearBondAttributes()
        for a1, a2 in bond_atom_pairs:
            self.bonds.append(Bonds.Bond((a1, a2)))
        for b in group.bonds:
            self.bonds.append(b)

    # construct positions of missing hydrogens
    def findHydrogenPositions(self):
        """
        Find reasonable positions for hydrogen atoms that have no
        position assigned.

        This method uses a heuristic approach based on standard geometry
        data. It was developed for proteins and DNA and may not give
        good results for other molecules. It raises an exception
        if presented with a topology it cannot handle.
        """
        self.setBondAttributes()
        try:
            unknown = {}
            for a in self.atoms:
                if a.position() is None:
                    if a.symbol != 'H':
                        raise ValueError('position of ' + a.fullName() + \
                                          ' is undefined')
                    bonded = a.bondedTo()[0]
                    unknown.setdefault(bonded, []).append(a)
            for a, list in unknown.items():
                bonded = a.bondedTo()
                n = len(bonded)
                known = [b for b in bonded if b.position() is not None]
                nb = len(list)
                try:
                    method = self._h_methods[a.symbol][n][nb]
                except KeyError:
                    raise ValueError("Can't handle this yet: " +
                                      a.symbol + ' with ' + `n` + ' bonds (' +
                                      a.fullName() + ').')
                method(self, a, known, list)
        finally:
            self.clearBondAttributes()

    # default C-H bond length and X-C-H angle
    _ch_bond = 1.09*Units.Ang
    _hch_angle = N.arccos(-1./3.)*Units.rad
    _nh_bond = 1.03*Units.Ang
    _hnh_angle = 120.*Units.deg
    _oh_bond = 0.95*Units.Ang
    _coh_angle = 114.9*Units.deg
    _sh_bond = 1.007*Units.Ang
    _csh_angle = 96.5*Units.deg

    def _C4oneH(self, atom, known, unknown):
        r = atom.position()
        n0 = (known[0].position()-r).normal()
        n1 = (known[1].position()-r).normal()
        n2 = (known[2].position()-r).normal()
        n3 = (n0 + n1 + n2).normal()
        unknown[0].setPosition(r-self._ch_bond*n3)

    def _C4twoH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        r2 = known[1].position()
        plane = Objects3D.Plane(r, r1, r2)
        axis = -((r1-r)+(r2-r)).normal()
        plane = plane.rotate(Objects3D.Line(r, axis), 90.*Units.deg)
        cone = Objects3D.Cone(r, axis, 0.5*self._hch_angle)
        sphere = Objects3D.Sphere(r, self._ch_bond)
        circle = sphere.intersectWith(cone)
        points = circle.intersectWith(plane)
        unknown[0].setPosition(points[0])
        unknown[1].setPosition(points[1])

    def _C4threeH(self, atom, known, unknown):
        self._tetrahedralH(atom, known, unknown, self._ch_bond)

    def _C3oneH(self, atom, known, unknown):
        r = atom.position()
        n1 = (known[0].position()-r).normal()
        n2 = (known[1].position()-r).normal()
        n3 = -(n1 + n2).normal()
        unknown[0].setPosition(r+self._ch_bond*n3)

    def _C3twoH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        others = filter(lambda a: a.symbol != 'H', known[0].bondedTo())
        r2 = others[0].position()
        try:
            plane = Objects3D.Plane(r, r1, r2)
        except ZeroDivisionError:
            # We get here if all three points are colinear.
            # Add a small random displacement as a fix.
            from MMTK.Random import randomPointInSphere
            plane = Objects3D.Plane(r, r1, r2 + randomPointInSphere(0.001))
        axis = (r-r1).normal()
        cone = Objects3D.Cone(r, axis, 0.5*self._hch_angle)
        sphere = Objects3D.Sphere(r, self._ch_bond)
        circle = sphere.intersectWith(cone)
        points = circle.intersectWith(plane)
        unknown[0].setPosition(points[0])
        unknown[1].setPosition(points[1])

    def _C2oneH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        x = r + self._ch_bond * (r - r1).normal()
        unknown[0].setPosition(x)

    def _N2oneH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        others = filter(lambda a: a.symbol != 'H', known[0].bondedTo())
        r2 = others[0].position()
        try:
            plane = Objects3D.Plane(r, r1, r2)
        except ZeroDivisionError:
            # We get here when all three points are colinear.
            # Add a small random displacement as a fix.
            from MMTK.Random import randomPointInSphere
            plane = Objects3D.Plane(r, r1, r2 + randomPointInSphere(0.001))
        axis = (r-r1).normal()
        cone = Objects3D.Cone(r, axis, 0.5*self._hch_angle)
        sphere = Objects3D.Sphere(r, self._nh_bond)
        circle = sphere.intersectWith(cone)
        points = circle.intersectWith(plane)
        unknown[0].setPosition(points[0])

    def _N3oneH(self, atom, known, unknown):
        r = atom.position()
        n1 = (known[0].position()-r).normal()
        n2 = (known[1].position()-r).normal()
        n3 = -(n1 + n2).normal()
        unknown[0].setPosition(r+self._nh_bond*n3)

    def _N3twoH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        others = filter(lambda a: a.symbol != 'H', known[0].bondedTo())
        r2 = others[0].position()
        plane = Objects3D.Plane(r, r1, r2)
        axis = (r-r1).normal()
        cone = Objects3D.Cone(r, axis, 0.5*self._hnh_angle)
        sphere = Objects3D.Sphere(r, self._nh_bond)
        circle = sphere.intersectWith(cone)
        points = circle.intersectWith(plane)
        unknown[0].setPosition(points[0])
        unknown[1].setPosition(points[1])

    def _N4threeH(self, atom, known, unknown):
        self._tetrahedralH(atom, known, unknown, self._nh_bond)

    def _N4twoH(self, atom, known, unknown):
        r = atom.position()
        r1 = known[0].position()
        r2 = known[1].position()
        plane = Objects3D.Plane(r, r1, r2)
        axis = -((r1-r)+(r2-r)).normal()
        plane = plane.rotate(Objects3D.Line(r, axis), 90.*Units.deg)
        cone = Objects3D.Cone(r, axis, 0.5*self._hnh_angle)
        sphere = Objects3D.Sphere(r, self._nh_bond)
        circle = sphere.intersectWith(cone)
        points = circle.intersectWith(plane)
        unknown[0].setPosition(points[0])
        unknown[1].setPosition(points[1])

    def _N4oneH(self, atom, known, unknown):
        r = atom.position()
        n0 = (known[0].position()-r).normal()
        n1 = (known[1].position()-r).normal()
        n2 = (known[2].position()-r).normal()
        n3 = (n0 + n1 + n2).normal()
        unknown[0].setPosition(r-self._nh_bond*n3)

    def _O2(self, atom, known, unknown):
        others = known[0].bondedTo()
        for a in others:
            r = a.position()
            if a != atom and r is not None: break
        dihedral = 180.*Units.deg
        self._findPosition(unknown[0], atom.position(), known[0].position(), r,
                           self._oh_bond, self._coh_angle, dihedral)

    def _S2(self, atom, known, unknown):
        c2 = filter(lambda a: a.symbol == 'C', known[0].bondedTo())[0]
        self._findPosition(unknown[0], atom.position(), known[0].position(),
                           c2.position(),
                           self._sh_bond, self._csh_angle,
                           180.*Units.deg)

    def _tetrahedralH(self, atom, known, unknown, bond):
        r = atom.position()
        n = (known[0].position()-r).normal()
        cone = Objects3D.Cone(r, n, N.arccos(-1./3.))
        sphere = Objects3D.Sphere(r, bond)
        circle = sphere.intersectWith(cone)
        others = filter(lambda a: a.symbol != 'H', known[0].bondedTo())
        others.remove(atom)
        other = others[0]
        ref = (Objects3D.Plane(circle.center, circle.normal) \
               .projectionOf(other.position())-circle.center).normal()
        p0 = circle.center + ref*circle.radius
        p0 = Objects3D.rotatePoint(p0,
                                  Objects3D.Line(circle.center, circle.normal),
                                  60.*Units.deg)
        p1 = Objects3D.rotatePoint(p0,
                                  Objects3D.Line(circle.center, circle.normal),
                                  120.*Units.deg)
        p2 = Objects3D.rotatePoint(p1,
                                  Objects3D.Line(circle.center, circle.normal),
                                  120.*Units.deg)
        unknown[0].setPosition(p0)
        unknown[1].setPosition(p1)
        unknown[2].setPosition(p2)

    def _findPosition(self, unknown, a1, a2, a3, bond, angle, dihedral):
        sphere = Objects3D.Sphere(a1, bond)
        cone = Objects3D.Cone(a1, a2-a1, angle)
        plane = Objects3D.Plane(a3, a2, a1)
        plane = plane.rotate(Objects3D.Line(a1, a2-a1), dihedral)
        points = sphere.intersectWith(cone).intersectWith(plane)
        for p in points:
            if (a1-a2).cross(p-a1)*(plane.normal) > 0:
                unknown.setPosition(p)
                break

    _h_methods = {'C': {4: {3: _C4threeH,
                            2: _C4twoH,
                            1: _C4oneH},
                        3: {2: _C3twoH,
                            1: _C3oneH},
                        2: {1: _C2oneH}},
                  'N': {4: {3: _N4threeH,
                            2: _N4twoH,
                            1: _N4oneH},
                        3: {2: _N3twoH,
                            1: _N3oneH},
                        2: {1: _N2oneH}},
                  'O': {2: {1: _O2}},
                  'S': {2: {1: _S2}},
                  }


class ChainMolecule(Molecule):

    """
    ChainMolecule

    ChainMolecules are generated by a MoleculeFactory from
    templates that have a sequence attribute.
    """

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, item):
        return getattr(self, self.sequence[item])


class Crystal(CompositeChemicalObject, ChemicalObject):

    def __init__(self, blueprint, _memo = None, **properties):
        ChemicalObject.__init__(self, blueprint, _memo)
        properties = copy.copy(properties)
        CompositeChemicalObject.__init__(self, properties)
        self.bonds = Bonds.BondList(self.bonds)

    blueprintclass = Database.BlueprintCrystal

    def _subunits(self):
        return self.groups

    def _descriptionSpec(self):
        return "X", None

class Complex(CompositeChemicalObject, ChemicalObject):

    """
    Complex

    A complex is an assembly of molecules that are not connected by
    chemical bonds.
    """

    def __init__(self, complex_spec, _memo = None, **properties):
        """
        :param complex_spec: a string (not case sensitive) that specifies
                              the complex name in the chemical database
        :type complex_spec: str
        :keyword position: the position of the center of mass of the complex
        :type position: Scientific.Geometry.Vector
        :keyword name: a name given to the complex
        :type name: str
        :keyword configuration: the name of a configuration listed in the
                                database definition of the complex, which
                                is used to initialize the atom positions.
                                If no configuration is specified, the
                                configuration named "default" will be used,
                                if it exists. Otherwise the atom positions
                                are undefined.
        :type configuration: str
        """
        ChemicalObject.__init__(self, complex_spec, _memo)
        properties = copy.copy(properties)
        CompositeChemicalObject.__init__(self, properties)

    blueprintclass = Database.BlueprintComplex

    def recreateAtomList(self):
        self.atoms = []
        for m in self.molecules:
            self.atoms.extend(m.atoms)

    def bondedUnits(self):
        return self.molecules

    def setBondAttributes(self):
        for m in self.molecules:
            m.setBondAttributes()

    def clearBondAttributes(self):
        for m in self.molecules:
            m.clearBondAttributes()

    def _subunits(self):
        return self.molecules

    def _descriptionSpec(self):
        return "C", None

    def writeXML(self, file, memo, toplevel=1):
        if self.type is None:
            name = 'm' + `memo['counter']`
            memo['counter'] = memo['counter'] + 1
            memo[id(self)] = name
            for molecule in self.molecules:
                molecule.writeXML(file, memo, 0)
            file.write('<molecule id="%s">\n' % name)
            for molecule in self.molecules:
                file.write('  <molecule ref="%s"/>\n' % memo[id(molecule)])
            file.write('</molecule>\n')
        else:
            ChemicalObject.writeXML(self, file, memo, toplevel)
        if toplevel:
            return ['<molecule ref="%s"/>' % name]
        else:
            return None

    def getXMLAtomOrder(self):
        if self.type is None:
            atoms = []
            for molecule in self.molecules:
                atoms.extend(molecule.getXMLAtomOrder())
            return atoms
        else:
            return ChemicalObject.getXMLAtomOrder(self)

Database.registerInstanceClass(Atom.blueprintclass, Atom)
Database.registerInstanceClass(Group.blueprintclass, Group)
Database.registerInstanceClass(Molecule.blueprintclass, Molecule)
Database.registerInstanceClass(Crystal.blueprintclass, Crystal)
Database.registerInstanceClass(Complex.blueprintclass, Complex)


class AtomCluster(CompositeChemicalObject, ChemicalObject):

    """
    An agglomeration of atoms

    An atom cluster acts like a molecule without any bonds or atom
    properties. It can be used to represent a group of atoms that
    are known to form a chemical unit but whose chemical properties
    are not sufficiently known to define a molecule.
    """

    def __init__(self, atoms, **properties):
        """
        :param atoms: a list of atoms in the cluster
        :type atoms: list
        :keyword position: the position of the center of mass of the cluster
        :type position: Scientific.Geometry.Vector
        :keyword name: a name given to the cluster
        :type name: str
        """
        self.atoms = list(atoms)
        self.parent = None
        self.name = ''
        self.type = None
        for a in self.atoms:
            if a.parent is not None:
                raise ValueError(repr(a)+' is part of ' + repr(a.parent))
            a.parent = self
            if a.name != '':
                setattr(self, a.name, a)
        properties = copy.copy(properties)
        CompositeChemicalObject.__init__(self, properties)
        self.bonds = Bonds.BondList([])

    def bondedTo(self, atom):
        return []

    def setBondAttributes(self):
        pass

    def clearBondAttributes(self):
        pass

    def _subunits(self):
        return []

    def _description(self, tag, index_map, toplevel):
        s = 'AC(' + `self.name` + ',['
        for a in self.atoms:
            s = s + a._description(tag, index_map, 1) + ','
        s = s + ']'
        constraints = self._distanceConstraintList()
        if constraints:
            s = s + ',dc=['
            if index_map is None:
                for c in constraints:
                    s = s + '(%d,%d,%f),' % (c[0].index, c[1].index, c[2])
            else:
                for c in constraints:
                    s = s + '(%d,%d,%f),' % (index_map[c[0].index],
                                             index_map[c[1].index], c[2])
            s = s + ']'
        return s + ')'
