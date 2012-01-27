# This module manages the chemical database.
#
# Written by Konrad Hinsen
#

"""
Management of the chemical database

The database contains definitions for atoms, groups, molecules, and
complexes. Each definition is a Python file that is executed in the
global environment of the module 'XEnvironment' (X standing for
'Atom', 'Group', 'Molecule', 'Complex', 'Protein', or 'Crystal').
Definitions made in that file will end up as attributes of an object
that is later used as a blueprint to create chemical objects.
"""

__docformat__ = 'restructuredtext'

from MMTK import Utility
import copy
import os
import sys

#
# Find database path
#
try:
    path = os.environ['MMTKDATABASE'].split()
except KeyError:
    path = ['~/.mmtk/Database',
            os.path.join(os.path.split(__file__)[0], 'Database')]

for i in range(len(path)):
    if not Utility.isURL(path[i]):
        path[i] = os.path.expanduser(path[i])

#
# Some miscellaneous functions for use by other modules
#
def databasePath(filename, directory, try_direct = False):
    if Utility.isURL(filename):
        return filename
    filename = os.path.expanduser(filename)
    if try_direct and os.path.exists(filename):
        return os.path.normcase(filename)
    entries = []
    if os.path.split(filename)[0] == '':
        for p in path:
            if Utility.isURL(p):
                url = Utility.joinURL(p, directory+'/'+filename)
                if Utility.checkURL(url):
                    entries.append(url)
            else:
                full_name = os.path.join(os.path.join(p, directory), filename)
                if os.path.exists(full_name):
                    entries.append(os.path.normcase(full_name))
    if len(entries) == 0:
        raise IOError("Database entry %s/%s not found" % (directory, filename))
    else:
        if len(entries) > 1:
            Utility.warning("multiple database entries for %s/%s, using first one"
                            % (directory, filename))
            for e in entries:
                sys.stderr.write(e+'\n')
        return entries[0]
        
def PDBPath(filename):
    return databasePath(filename, 'PDB', True)

def addDatabaseDirectory(directory):
    "Add a directory to the database search path"
    if not Utility.isURL(directory):
        directory = os.path.expanduser(directory)
        path.append(directory)

#
# The class that represents a database. There will be one instance
# for atoms, one for groups etc.
# 
class Database(object):

    def __init__(self, directory, type_constructor):
        self.directory = directory
        self.type_constructor = type_constructor
        self.types = {}

    def findType(self, name):
        name = name.lower()
        if not self.types.has_key(name):
            filename = databasePath(name, self.directory, False)
            self.types[name] = self.type_constructor(filename, name)
        return self.types[name]

#
# The base class for all the type classes. It defines how definitions
# are loaded from the database.
#
class ChemicalObjectType(object):

    def __init__(self, filename, database_name, module, instancevars):
        self.filename = filename
        self.database_name = database_name
        file_text = Utility.readURL(filename)
        newvars = {}
        exec file_text in vars(module), newvars
        for name, value in newvars.items():
            setattr(self, name, value)
        self.parent = None
        if not hasattr(self, 'instance'): self.instance = []
        for attr in instancevars+('parent',):
            if not hasattr(self, attr): setattr(self, attr, [])
            if attr not in self.instance: self.instance.append(attr)
        attributes = vars(self).items()
        attributes.sort(lambda a, b: cmp(a[0], b[0]))
        for name, object in attributes:
            if hasattr(object, 'is_instance_var'):
                if name not in self.instance:
                    self.instance.append(name)
                object.parent = self
                object.name = name
            if hasattr(object, 'object_list'):
                getattr(self, object.object_list).append(object)

    is_chemical_object_type = True

    def setReferences(self):
        atom_refs = []
        for i in range(len(self.atoms)):
            atom_refs.append(AtomReference(i))
        for attr in vars(self).items():
            if attr[0] not in self.instance:
                setattr(self, attr[0],
                        Utility.substitute(getattr(self, attr[0]),
                                           self.atoms, atom_refs))

    # Type objects are singletons, they are never copied
    def __copy__(self, memo = None):
        return self
    __deepcopy__ = __copy__

    # Pickle support. When pickled/unpickled through
    # MMTK.Utility.save and MMTK.Utility.load, type objects
    # are treated as external persistent objects. This is handled
    # by a specialized Pickler and Unpickler (see MMTK.Utility)
    # and the _restoreId() methods in the subclasses of
    # ChemicalObjectType. This mechanism maintains the singleton
    # nature of type objects.
    # When pickling using the unmodified pickle/cPickle routines,
    # there is no way to prevent type objects from being duplicated
    # at unpickle time. The best we can do is to make the newly created
    # type object a clone of the singleton object. This is handled
    # by the __setstate__ and __getstate__ methods that follow

    def __getstate__(self):
        return self._restoreId()

    def __setstate__(self, state):
        import sys
        singleton_type = eval("MMTK."+state, sys.modules)
        self.__dict__.update(singleton_type.__dict__)


    def writeXML(self, file, memo):
        if memo.get(id(self), None) is not None:
            return
        try:
            name = self.name
        except AttributeError:
            name = 'm' + `memo['counter']`
            memo['counter'] = memo['counter'] + 1
        memo[id(self)] = name
        atoms = copy.copy(self.atoms)
        bonds = copy.copy(self.bonds)
        for group in self.groups:
            group.type.writeXML(file, memo)
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

    def getXMLAtomOrder(self):
        atoms = copy.copy(self.atoms)
        atom_names = []
        for group in self.groups:
            for atom in group.atoms:
                atoms.remove(atom)
            group_name = self.relativeName(group)
            for atom in group.type.getXMLAtomOrder():
                atom_names.append(group_name + ':' + atom)
        for atom in atoms:
            atom_names.append(self.relativeName(atom))
        return atom_names

    def relativeName(self, object):
        name = ''
        while object is not None and object != self:
            for attr, value in object.parent.__dict__.items():
                if value is object:
                    name = attr + ':' + name
                    break
            object = object.parent
        return name[:-1]

#
# Atom type class
#
class AtomType(ChemicalObjectType):

    error = 'AtomTypeError'

    def __init__(self, filename, database_name):
        from MMTK import AtomEnvironment
        ChemicalObjectType.__init__(self, filename, database_name,
                                    AtomEnvironment, ())
        if not isinstance(self.mass, list):
            self.mass = [(self.mass, 100.)]
        total_probability = sum([m[1] for m in self.mass])
        if abs(total_probability-100.) > 1.e-4:
            raise self.error('Inconsistent mass specification: ' +
                              `total_probability-100.` + ' percent missing')
        self.average_mass = sum([m[0]*m[1] for m in self.mass])/100
        if not hasattr(self, 'pdbmap'):
            name = self.symbol.upper()
            self.pdbmap = [(name, {name: None})]

    def _restoreId(self):
        return 'Database.atom_types.findType("' + \
               self.database_name + '")'

    def writeXML(self, file, memo):
        file.write('<atom/>\n')

    def getXMLAtomOrder(self):
        return [self.name]

#
# Group type class
#
class GroupType(ChemicalObjectType):

    error = 'GroupTypeError'

    def __init__(self, filename, database_name):
        from MMTK import GroupEnvironment
        ChemicalObjectType.__init__(self, filename, database_name,
                                    GroupEnvironment,
                                    ('atoms', 'groups', 'bonds',
                                     'chain_links'))
        for g in self.groups:
            self.atoms = self.atoms + g.atoms
            self.bonds = self.bonds + g.bonds
        self.setReferences()

    def _restoreId(self):
        return 'Database.group_types.findType("' + \
               self.database_name + '")'

#
# Molecule type class
#
class MoleculeType(ChemicalObjectType):

    error = 'MoleculeTypeError'

    def __init__(self, filename, database_name):
        from MMTK import MoleculeEnvironment
        ChemicalObjectType.__init__(self, filename, database_name,
                                    MoleculeEnvironment,
                                    ('atoms', 'groups', 'bonds'))
        for g in self.groups:
            self.atoms = self.atoms + g.atoms
            self.bonds = self.bonds + g.bonds
        self.setReferences()

    def _restoreId(self):
        return 'Database.molecule_types.findType("' + \
               self.database_name + '")'

#
# Crystal type class
#
class CrystalType(ChemicalObjectType):

    error = 'CrystalTypeError'

    def __init__(self, filename, database_name):
        from MMTK import CrystalEnvironment
        ChemicalObjectType.__init__(self, filename, database_name,
                                    CrystalEnvironment,
                                    ('atoms', 'groups', 'molecules', 'bonds'))
        for g in self.groups:
            self.atoms = self.atoms + g.atoms
            self.bonds = self.bonds + g.bonds
        for m in self.molecules:
            self.atoms = self.atoms + m.atoms
            self.bonds = self.bonds + m.bonds
        self.setReferences()

    def _restoreId(self):
        return 'Database.crystal_types.findType("' + \
               self.database_name + '")'

#
# Complex type class
#
class ComplexType(ChemicalObjectType):

    error = 'ComplexTypeError'

    def __init__(self, filename, database_name):
        import ComplexEnvironment
        ChemicalObjectType.__init__(self, filename, database_name,
                                    ComplexEnvironment,
                                    ('atoms', 'molecules'))
        for m in self.molecules:
            self.atoms = self.atoms + m.atoms
        self.setReferences()

    def _restoreId(self):
        return 'Database.complex_types.findType("' + \
               self.database_name + '")'

#
# An atom reference object is substituted for all references to
# atom objects that are not in instance variables. A reference
# object contains only the number of the atom in the list of
# atoms of its parent object.
#
class AtomReference(object):

    def __init__(self, number):
        self.number = number

    def increaseBy(self, offset):
        self.number = self.number + offset

    def __repr__(self):
        return '<Atom number ' + `self.number` + '>'
    __str__ = __repr__

    def __eq__(self, other):
        if isinstance(other, AtomReference):
            return self.number == other.number
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, AtomReference):
            return self.number != other.number
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self.number)

#
# The base class for all types that just contain references to
# the file names. These are for objects that tend to be big
# and used only in small numbers.
#
class ReferenceType(object):

    def __init__(self, filename, database_name, environment):
        self.filename = filename
        self.database_name = database_name
        self.environment = environment

    def createObject(self, newvars):
        file_text = Utility.readURL(self.filename)
        exec file_text in vars(self.environment), newvars

class ProteinType(ReferenceType):

    def __init__(self, filename, database_name):
        from MMTK import ProteinEnvironment
        ReferenceType.__init__(self, filename, database_name,
                               ProteinEnvironment)

    def _restoreId(self):
        return 'Database.protein_types.findType("' + \
               self.database_name + '")'

#
# The five databases.
#
atom_types = Database('Atoms', AtomType)
group_types = Database('Groups', GroupType)
molecule_types = Database('Molecules', MoleculeType)
crystal_types = Database('Crystals', CrystalType)
complex_types = Database('Complexes', ComplexType)
protein_types = Database('Proteins', ProteinType)

#
# The following classes represent the chemical objects
# in the type blueprints. They contain no information
# in addition to references to an object type. They
# are needed only to establish a test of identity;
# e.g. each hydrogen atom in a molecule must be represented
# by a different object.
#
class BlueprintObject(object):

    def __init__(self, original, database, memo):
        if isinstance(original, basestring):
            original = database.findType(original)
            self.type = original
        elif hasattr(original, 'is_blueprint'):
            self.type = original.type
            if hasattr(original, 'name'):
                self.name = original.name
        else:
            self.type = original
        if memo is None: memo = {}
        memo[id(original)] = self
        for attr in self.type.instance:
            setattr(self, attr, _blueprintCopy(getattr(original, attr), memo))

    is_instance_var = True
    is_blueprint = True

    def __getattr__(self, attr):
        return getattr(self.type, attr)

    def __copy__(self, memo = None):
        return self
    __deepcopy__ = __copy__

class BlueprintAtom(BlueprintObject):
    def __init__(self, type, memo = None):
        BlueprintObject.__init__(self, type, atom_types, memo)
    object_list = 'atoms'

class BlueprintGroup(BlueprintObject):
    def __init__(self, type, memo = None):
        BlueprintObject.__init__(self, type, group_types, memo)
    object_list = 'groups'

class BlueprintMolecule(BlueprintObject):
    def __init__(self, type, memo = None):
        BlueprintObject.__init__(self, type, molecule_types, memo)
    object_list = 'molecules'

class BlueprintCrystal(BlueprintObject):
    def __init__(self, type, memo = None):
        BlueprintObject.__init__(self, type, crystal_types, memo)
    object_list = 'crystals'

class BlueprintComplex(BlueprintObject):
    def __init__(self, type, memo = None):
        BlueprintObject.__init__(self, type, complex_types, memo)
    object_list = 'complexes'

#
# Blueprint class corresponding to ReferenceType
#
class ReferenceBlueprint(object):

    def __init__(self, original, database):
        if isinstance(original, basestring):
            self.type = database.findType(original)
        elif hasattr(original, 'is_blueprint'):
            self.type = original.type
        else:
            self.type = type
        self.type.createObject(self.__dict__)

    is_blueprint = True

class BlueprintProtein(ReferenceBlueprint):

    def __init__(self, type):
        ReferenceBlueprint.__init__(self, type, protein_types)

#
# This function copies the appropriate attributes of
# a blueprint object.
#
def _blueprintCopy(object, memo):
    key = id(object)
    try:
        return memo[key]
    except KeyError:
        pass
    if isinstance(object, list):
        return [_blueprintCopy(o, memo) for o in object]
    if hasattr(object, 'is_blueprint'):
        new = object.__class__(object, memo)
    elif hasattr(object, '_blueprintCopy'):
        new = object._blueprintCopy(memo)
    else:
        new = object
    memo[key] = new
    return new

#
# The bond class just keeps track of the two atoms involved.
#
class BlueprintBond(object):

    def __init__(self, a1, a2):
        self.a1 = a1
        self.a2 = a2

    is_instance_var = True

    def _blueprintCopy(self, memo):
        return BlueprintBond(_blueprintCopy(self.a1, memo),
                             _blueprintCopy(self.a2, memo))

    object_list = 'bonds'

    def __copy__(self, memo = None):
        return self
    __deepcopy__ = __copy__

#
# The function "instantiate" returns the real object corresponding
# to a given blueprint object. A list of previously instantiated
# objects is kept to make sure that references to the same blueprint
# object don't create several real objects.
#
def instantiate(blueprint, memo):
    key = id(blueprint)
    try:
        return memo[key]
    except KeyError:
        pass
    if isinstance(blueprint, list):
        newobject = [instantiate(e, memo) for e in blueprint]
    else:
        try:
            newobject = _instanceclass[blueprint.__class__](blueprint, memo)
        except KeyError:
            newobject = blueprint
    memo[key] = newobject
    return newobject

#
# Register the instanceclasses for each blueprintclass
#
def registerInstanceClass(blueprint, instance):
    _instanceclass[blueprint] = instance

_instanceclass = {}
