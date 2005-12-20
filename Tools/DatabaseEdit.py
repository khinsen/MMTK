# Functions that help to create/modify the database
#
# Written by Konrad Hinsen
# last revision: 1999-9-2
#

_undocumented = 1

import Database
import copy, os, string

class TypeClone(Database.ChemicalObjectType):

    def __init__(self, filename, module, instancevars):
	self._filename = os.path.split(filename)[-1]
	module.Atom = BlueprintAtom
	module.Group = BlueprintGroup
	module.Molecule = BlueprintMolecule
	module.Bond = BlueprintBond
	Database.ChemicalObjectType.__init__(self, filename, module,
					     instancevars)
	self._written = 0

    def attributes(self):
	dict = copy.copy(self.__dict__)
	del dict['_filename']
	del dict['_written']
	del dict['instance']
	del dict['parent']
	del dict['filename']
	try:
	    del dict['amber_atom_type']
	except KeyError: pass
	try:
	    del dict['amber_charge']
	except KeyError: pass
	return dict

    def writeToFile(self, filename_converter):
	if not self._written:
	    filename = filename_converter(self._filename, self.class_name)
	    filename = Database.DatabasePath(filename, 'Output')
	    file = open(filename, 'w')
	    print '*** File ' + filename
	    dict = self.attributes()
	    try:
		atoms = copy.copy(dict['atoms'])
		del dict['atoms']
	    except KeyError:
		atoms = []
	    try:
		groups = copy.copy(dict['groups'])
		del dict['groups']
	    except KeyError:
		groups = []
	    try:
		molecules = copy.copy(dict['molecules'])
		del dict['molecules']
	    except KeyError:
		molecules = []
	    try:
		bonds = copy.copy(dict['bonds'])
		del dict['bonds']
	    except KeyError:
		bonds = []
	    values = []
	    for attr, value in dict.items():
		if hasattr(value, 'is_blueprint'):
		    file.write(attr + ' = ')
		    self.writeValue(file, value, filename_converter)
		    value.setString(attr)
		    file.write('\n')
		    eval(value.object_list).remove(value)
		    del dict[attr]
	    if atoms:
		file.write('atoms = ')
		self.writeValue(file, atoms, filename_converter)
		file.write('\n')
	    if groups:
		file.write('groups = ')
		self.writeValue(file, groups, filename_converter)
		file.write('\n')
	    if molecules:
		file.write('molecules = ')
		self.writeValue(file, molecules, filename_converter)
		file.write('\n')
	    if bonds:
		file.write('bonds = ')
		self.writeValue(file, bonds, filename_converter)
		file.write('\n')
	    for attr, value in dict.items():
		file.write(attr + ' = ')
		self.writeValue(file, value, filename_converter)
		file.write('\n')
		del dict[attr]
	    file.close()
	    self._written = 1

    def writeValue(self, file, value, filename_converter):
	if hasattr(value, 'is_blueprint'):
	    value.writeToFile(self, file, filename_converter)
	elif type(value) == type([]):
	    file.write('[')
	    for element in value:
		self.writeValue(file, element, filename_converter)
		file.write(', ')
	    file.write(']')
	elif type(value) == type(()):
	    file.write('(')
	    for element in value:
		self.writeValue(file, element, filename_converter)
		file.write(', ')
	    file.write(')')
	elif type(value) == type({}):
	    file.write('{')
	    for k, v in value.items():
		self.writeValue(file, k, filename_converter)
		file.write(': ')
		self.writeValue(file, v, filename_converter)
		file.write(', ')
	    file.write('}')
	else:
	    file.write(repr(value))

    def parentName(self, parent):
	if parent != self:
	    raise TypeError
	return ''

    def deleteAtoms(self, identifier):
	dict = self.attributes()
	try: del dict['atoms']
	except KeyError: pass
	try: del dict['groups']
	except KeyError: pass
	try: del dict['molecules']
	except KeyError: pass
	try: del dict['bonds']
	except KeyError: pass
	for attr, value in dict.items():
	    if hasattr(value, 'is_blueprint'):
		if value.__class__ == BlueprintAtom and identifier(value):
		    delattr(self, attr)
	    else:
		setattr(self, attr, _deleteAtoms(value, identifier))
	self.atoms = filter(lambda a, i=identifier: not i(a), self.atoms)
	self.bonds = filter(lambda b, i=identifier: not (i(b.a1) or i(b.a2)),
			    self.bonds)

def _deleteAtoms(value, identifier):
    if type(value) == type([]):
	return map(_deleteAtoms, value, len(value)*[identifier])
    if type(value) == type(()):
	return tuple(map(_deleteAtoms, value, len(value)*[identifier]))
    if type(value) == type({}):
	new = {}
	for k, v in value.items():
	    if not identifier(k) and not identifier(v):
		new[k] = _deleteAtoms(v, identifier)
	return new
    else:
	return value

class AtomTypeClone(TypeClone):

    error = 'AtomTypeError'

    def __init__(self, filename):
	import AtomEnvironment
	TypeClone.__init__(self, filename, AtomEnvironment, ())

    class_name = 'Atom'

class GroupTypeClone(TypeClone):

    error = 'GroupTypeError'

    def __init__(self, filename):
	import GroupEnvironment
	TypeClone.__init__(self, filename, GroupEnvironment,
			   ('atoms', 'groups', 'bonds'))

    class_name = 'Group'

class MoleculeTypeClone(TypeClone):

    error = 'MoleculeTypeError'

    def __init__(self, filename):
	import MoleculeEnvironment
	TypeClone.__init__(self, filename, MoleculeEnvironment,
			   ('atoms', 'groups', 'bonds'))

    class_name = 'Molecule'


atom_types = Database.Database('Atoms', AtomTypeClone)
group_types = Database.Database('Groups', GroupTypeClone)
molecule_types = Database.Database('Molecules', MoleculeTypeClone)
#crystal_types = Database.Database('Crystals', CrystalTypeClone)
#complex_types = Database.Database('Complexes', ComplexTypeClone)

class BlueprintObjectClone:

    def parentName(self, parent):
	if hasattr(self, 'parent') and self.parent != parent:
	    return self.parent.parentName() + '.' + self._string
	else:
	    return self._string

    def setString(self, string):
	self._string = string

    def writeToFile(self, parent, file, filename_converter):
	if hasattr(self, 'parent') and self.parent != parent:
	    file.write(self.parent.parentName(parent) + '.')
	if self._string is not None:
	    file.write(self._string)
	else:
	    self.writeRepresentation(parent, file, filename_converter)

class BlueprintObject(Database.BlueprintObject, BlueprintObjectClone):

    def __init__(self, original, database, memo):
	Database.BlueprintObject.__init__(self, original, database, memo)
	self._string = None

    def writeRepresentation(self, parent, file, filename_converter):
	self.type.writeToFile(filename_converter)
	if hasattr(self.parent, 'is_blueprint'):
	    file.write(self.name)
	else:
	    filename = filename_converter(self.type._filename, self.class_name)
	    file.write(self.class_name + '(' + repr(filename) + ')')


class BlueprintAtom(BlueprintObject):
    def __init__(self, type, memo = None):
	BlueprintObject.__init__(self, type, atom_types, memo)
    object_list = 'atoms'
    class_name = 'Atom'

class BlueprintGroup(BlueprintObject):
    def __init__(self, type, memo = None):
	BlueprintObject.__init__(self, type, group_types, memo)
    object_list = 'groups'
    class_name = 'Group'

class BlueprintMolecule(BlueprintObject):
    def __init__(self, type, memo = None):
	BlueprintObject.__init__(self, type, molecule_types, memo)
    object_list = 'molecules'
    class_name = 'Molecule'

class BlueprintCrystal(BlueprintObject):
    def __init__(self, type, memo = None):
	BlueprintObject.__init__(self, type, crystal_types, memo)
    object_list = 'crystals'

class BlueprintComplex(BlueprintObject):
    def __init__(self, type, memo = None):
	BlueprintObject.__init__(self, type, complex_types, memo)
    object_list = 'complexes'

class BlueprintBond(BlueprintObjectClone):

    def __init__(self, a1, a2):
	self.a1 = a1
	self.a2 = a2
	self._string = None

    is_blueprint = 1
    object_list = 'bonds'

    def writeRepresentation(self, parent, file, filename_converter):
	file.write('Bond(')
	self.a1.writeToFile(parent, file, filename_converter)
	file.write(', ')
	self.a2.writeToFile(parent, file, filename_converter)
	file.write(')')

#
# Delete hydrogens in all amino acids
#

residue_names = ['alanine',
		 'arginine',
		 'asparagine',
		 'aspartic_acid',
		 'cysteine',
		 'cystine_ss',
		 'glutamine',
		 'glutamic_acid',
		 'glycine', 
		 'histidine',
		 'isoleucine',
		 'leucine',
		 'lysine',
		 'methionine',
		 'phenylalanine',
		 'proline', 
		 'serine', 
		 'threonine',
		 'tryptophan',
		 'tyrosine',
		 'valine',
		 'histidine_deltah',
		 'histidine_epsilonh',
		 'histidine_plus']
nt = map(lambda n: n + '_nt', residue_names)
ct = map(lambda n: n + '_ct', residue_names)
residue_names = residue_names + nt + ct

def noh(x, class_name):
    if class_name == 'Atom':
	return x
    else:
	return x + '_noh'

def uni(x, class_name):
    if class_name == 'Atom':
	return x
    else:
	return x + '_uni'

def uni2(x, class_name):
    if class_name == 'Atom':
	return x
    else:
	return x + '_uni2'

def hydrogenIdentifier(atom):
    return hasattr(atom, 'is_blueprint') and atom.__class__ == BlueprintAtom \
	   and atom.type.symbol == 'H'

def amberToCHARMM19(atom):
    return hasattr(atom, 'is_blueprint') and atom.__class__ == BlueprintAtom \
	   and atom.type.symbol == 'LP'

def markedIdentifier(atom):
    return hasattr(atom, 'is_blueprint') and atom.__class__ == BlueprintAtom \
	   and hasattr(atom, 'delete')

def markNonpolar(group):
    for b in group.bonds:
	if b.a1.symbol == 'C' and b.a2.symbol == 'H':
	    b.a2.delete = 1
	    print b.a2.name
	if b.a1.symbol == 'H' and b.a2.symbol == 'C':
	    b.a1.delete = 1
	    print b.a1.name

for name in residue_names:
    group_types.findType(name)
    group_types.findType(name+'_uni')
    group_types.findType(name+'_uni2')
    group_types.findType(name+'_noh')

def numberedHydrogen(name):
    return name[0] == 'H' and name[-1] in string.digits

def newname(name):
    return name[-1]+name[:-1]

for t in group_types.types.values():
    print t.name
    if hasattr(t, 'pdbmap'):
	try:
	    pdbalt = t.pdb_alternative
	except AttributeError:
	    pdbalt = {}
	for n1, n2 in pdbalt.items():
	    if numberedHydrogen(n2):
		pdbalt[n1] = newname(n2)
	pdbdict = t.pdbmap[0][1]
	for name, object in pdbdict.items():
	    print name
	    if numberedHydrogen(name):
		new = newname(name)
		del pdbdict[name]
		pdbdict[new] = object
		pdbalt[name] = new
	if len(pdbalt) > 0:
	    t.pdb_alternative = pdbalt

for t in group_types.types.values():
    t.writeToFile(lambda a, b: a)

## for t in group_types.types.values():
##     t.deleteAtoms(amberToCHARMM19)
## for name in residue_names:
##     t = group_types.findType(name)
##     t.writeToFile(uni2)

## for t in group_types.types.values():
##     t.deleteAtoms(hydrogenIdentifier)
## for name in residue_names:
##     t = group_types.findType(name)
##     t.writeToFile(noh)

## for name, t in group_types.types.items():
##     if string.find(name, 'peptide') != -1 or \
##        string.find(name, 'sidechain') != -1:
## 	print name
## 	markNonpolar(t)
## group_types.types['gly_sidechain'].H_alpha_3.delete = 1
## group_types.types['gly_nt_sidechain'].H_alpha_3.delete = 1

## for t in group_types.types.values():
##     t.deleteAtoms(markedIdentifier)

## for name in residue_names:
##     t = group_types.findType(name)
##     t.writeToFile(uni2)
