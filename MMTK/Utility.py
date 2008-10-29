# This module contains useful functions that are needed in various
# places.
#
# Written by Konrad Hinsen
# last revision: 2008-10-29
#

_undocumented = 1

import Database
import os, string, sys, types
from Scientific import N

# Constants

undefined_limit = 1.e30
undefined = 10.*undefined_limit

# Error

class MMTKError(Exception):
    pass

#
# Unique ID store.
# Certain objects must have a unique ID to allow unique ordering,
# e.g. of atoms in a bond. Usually the id() function is sufficient,
# except when running in a parallel environment, where the IDs
# must be identical across all processors as long as the code is
# identical.
#
class StandardUniqueIDGenerator:

    def __call__(self, object):
        return id(object)

    def registerObject(self, object):
        pass

class DeterministicUniqueIDGenerator:

    def __init__(self):
        self.number = 0
        self.id = {}

    def __call__(self, object):
        return self.id[object]

    def registerObject(self, object):
        self.id[object] = self.number
        self.number = self.number + 1

parallel = 1
try:
    from Scientific.MPI import world
    if world is None:
        parallel = 0
    elif world.size == 1:
        parallel = 0
    del world
except ImportError:
    parallel = 0
if parallel:
    uniqueID = DeterministicUniqueIDGenerator()
else:
    uniqueID = StandardUniqueIDGenerator()
del parallel

#
# This function substitutes all references to one object by references to
# another object in a given object recursively. The  substitution is
# specified by a dictionary.
# Attention: For some object types, substitution is destructive!
#
def substitute(obj, *exchange):
    if len(exchange) == 1:
        exchange = exchange[0]
    else:
        dict = {}
        map(lambda k, v, d=dict: _put(d, k, v), exchange[0], exchange[1])
        exchange = dict
    if type(obj) == types.ListType:
        return map(lambda e, x=exchange: substitute(e, x), obj)
    if type(obj) == types.TupleType:
        return tuple(map(lambda e, x=exchange: substitute(e, x), obj))
    elif type(obj) == types.DictionaryType:
        newdict = {}
        for key, value in obj.items():
            newdict[substitute(key, exchange)] = substitute(value, exchange)
        return newdict
    elif not isinstance(obj, type) and hasattr(obj, '_substitute'):
        for attr in dir(obj):
            if attr[:2] != '__':
                setattr(obj, attr, substitute(getattr(obj, attr), exchange))
        return obj
    elif exchange.has_key(obj):
        return exchange[obj]
    else:
        return obj

def _put(dict, key, value):
    dict[key] = value

#
# Return a unique attribute name
#
_unique_attributes = 0
def uniqueAttribute():
    global _unique_attributes
    _unique_attributes = (_unique_attributes + 1) % 10000
    return '_' + `_unique_attributes` + '__'

#
# Return a list of all pairs of objects in a given list
#
def pairs(list):
    p = []
    for i in range(len(list)):
        for j in range(i+1,len(list)):
            p.append((list[i], list[j]))
    return p

#
# Type check for sequence objects
#
def isSequenceObject(obj):
    t = type(obj)
    return t == types.ListType or t == types.TupleType \
           or (t == types.InstanceType and hasattr(obj, '__getitem__') \
               and hasattr(obj, '__len__'))
#
# Check if an object represents a well-defined position
#
def isDefinedPosition(p):
    if p is None:
        return 0
    if N.add.reduce(N.greater(p.array, undefined_limit)) > 0:
        return 0
    return 1

#
# Print a warning with reasonable line breaks.
#
def warning(text):
    words = string.split(text)
    text = 'Warning:'
    l = len(text)
    while words:
        lw = len(words[0])
        if l + lw + 1 < 60:
            text = text + ' ' + words[0]
            l = l + lw
        else:
            text = text + '\n' + 9*' ' + words[0]
            l = lw + 9
        words = words[1:]
    sys.stderr.write(text+"\n")

#
# Pickler and unpickler taking care of non-pickled objects
#

try:
    array_package = N.package
except AttributeError:
    array_package = 'Numeric'
if array_package == 'Numeric':
    BasePickler = N.Pickler
    BaseUnpickler = N.Unpickler
else:
    from pickle import Pickler as BasePickler
    from pickle import Unpickler as BaseUnpickler
del array_package

class Pickler(BasePickler):

    def persistent_id(self, obj):
        if hasattr(obj, 'is_chemical_object_type'):
            id = obj._restoreId()
            return id
        else:
            return None

class _EmptyClass:
    pass

class Unpickler(BaseUnpickler):

    def __init__(self, *args):
        BaseUnpickler.__init__(self, *args)
        self.dispatch['i'] = Unpickler.load_inst

    def persistent_load(self, id):
        return eval(id)

    def find_class(self, module, name):
        env = {}
        try:
            exec 'from %s import %s' % (module, name) in env
            klass = env[name]
        except ImportError:
            from NewModuleNames import new_name
            nmodule, nname = new_name.get((module, name),
                                          ("MMTK."+module, name))
            try:
                exec 'from %s import %s' % (nmodule, nname) in env
                klass = env[nname]
            except ImportError:
                raise SystemError("Failed to import class %s from module %s" %
                                   (name, module))
        return klass

    # Modified load_inst removes argument lists for classes that used
    # to have __getinitargs__ but don't have it any more.
    # This makes it possible to read pickle files from older MMTK versions.
    def load_inst(self):
        k = self.marker()
        args = tuple(self.stack[k+1:])
        del self.stack[k:]
        module = self.readline()[:-1]
        name = self.readline()[:-1]
        klass = self.find_class(module, name)
        instantiated = 0
        if ((not args or hasattr(klass, "__had_initargs__"))
            and type(klass) is types.ClassType
            and not hasattr(klass, "__getinitargs__")):
            try:
                value = _EmptyClass()
                value.__class__ = klass
                instantiated = 1
            except RuntimeError:
                # In restricted execution, assignment to inst.__class__ is
                # prohibited
                pass
        if not instantiated:
            try:
                if not hasattr(klass, '__safe_for_unpickling__'):
                    raise UnpicklingError('%s is not safe for unpickling' %
                                          klass)
                value = apply(klass, args)
            except TypeError, err:
                raise TypeError, "in constructor for %s: %s" % (
                    klass.__name__, str(err)), sys.exc_info()[2]
        self.append(value)


#
# General routines for writing objects to files and reading them back
#
def save(obj, filename):
    """Writes |obj| to a newly created file with the name |filename|,
    for later retrieval by 'load()'."""
    import ChemicalObjects
    filename = os.path.expanduser(filename)
    file = open(filename, 'wb')
    if ChemicalObjects.isChemicalObject(obj):
        parent = obj.parent
        obj.parent = None
        Pickler(file).dump(obj)
        obj.parent = parent
    else:
        Pickler(file).dump(obj)
    file.close()

def load(filename):
    """Loads the file indicated by |filename|, which must have been produced
    by 'save()', and returns the object stored in that file."""
    filename = os.path.expanduser(filename)
    file = open(filename, 'rb')
    obj = Unpickler(file).load()
    file.close()
    return obj

#
# URL related functions
#
def isURL(filename):
    return string.find(filename, ':/') > 1

def joinURL(url, filename):
    if url[-1] == '/':
        return url+filename
    else:
        return url+'/'+filename

def checkURL(filename):
    if isURL(filename):
        import urllib
        try:
            urllib.urlopen(filename)
            return 1
        except IOError:
            return 0
    else:
        return os.path.exists(filename)

def readURL(filename):
    try:
        if isURL(filename):
            import urllib
            file = urllib.urlopen(filename)
        else:
            file = open(filename)
    except IOError, details:
        if details[0] == 2:
            print "File " + filename + " not found."
        raise IOError(details)
    data = file.read()
    file.close()
    return data
