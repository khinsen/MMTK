# This module contains useful functions that are needed in various
# places.
#
# Written by Konrad Hinsen
#

import os, sys
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
        self.number += self.number

parallel = True
try:
    from Scientific.MPI import world
    if world is None:
        parallel = False
    elif world.size == 1:
        parallel = False
    del world
except ImportError:
    parallel = False
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
        exchange = dict(zip(exchange[0], exchange[1]))
    if isinstance(obj, list):
        return [substitute(e, exchange) for e in obj]
    if isinstance(obj, tuple):
        return tuple([substitute(e, exchange) for e in obj])
    if isinstance(obj, dict):
        newdict = {}
        for key, value in obj.items():
            newdict[substitute(key, exchange)] = substitute(value, exchange)
        return newdict
    if hasattr(obj, '_substitute'):
        for attr in vars(obj).keys():
            setattr(obj, attr, substitute(getattr(obj, attr), exchange))
        return obj
    try:
        return exchange[obj]
    except KeyError:
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
# Ensure that items in a pair are ordered
#
def normalizePair(pair):
    i, j = pair
    if i > j:
        return j, i
    else:
        return i, j

#
# Return an iterator over all pairs of objects in a given sequence
#
def pairs(seq):
    n = len(seq)
    for i in range(n):
        a = seq[i]
        for j in range(i+1, n):
            b = seq[j]
            yield (a, b)

#
# Return an iterator over all ordered pairs of objects in a given sequence
#
def orderedPairs(seq):
    n = len(seq)
    for i in range(n):
        a = seq[i]
        for j in range(i+1, n):
            b = seq[j]
            if a > b:
                yield (b, a)
            else:
                yield (a, b)

#
# Type check for sequence objects
#
def isSequenceObject(obj):
    try:
        it = iter(obj)
        return True
    except:
        return False

#
# Check if an object represents a well-defined position
#
def isDefinedPosition(p):
    if p is None:
        return False
    if N.add.reduce(N.greater(p.array, undefined_limit)) > 0:
        return False
    return True

#
# Print a warning with reasonable line breaks.
#
def warning(text):
    words = text.split()
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
        from MMTK import Database
        return eval(id)

    def find_class(self, module, name):
        env = {}
        exec 'from %s import %s' % (module, name) in env
        return env[name]

    # Modified load_inst removes argument lists for classes that used
    # to have __getinitargs__ but don't have it any more.
    # This makes it possible to read pickle files from older MMTK versions.
    def load_inst(self):
        import types
        k = self.marker()
        args = tuple(self.stack[k+1:])
        del self.stack[k:]
        module = self.readline()[:-1]
        name = self.readline()[:-1]
        klass = self.find_class(module, name)
        instantiated = False
        if ((not args or hasattr(klass, "__had_initargs__"))
            and type(klass) is types.ClassType
            and not hasattr(klass, "__getinitargs__")):
            try:
                value = _EmptyClass()
                value.__class__ = klass
                instantiated = True
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
    return filename.find(':/') > 1

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
            return True
        except IOError:
            return False
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
