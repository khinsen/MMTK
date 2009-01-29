# This module implements classes for nucleotide chains.
#
# Written by Konrad Hinsen
# last revision: 2009-1-29
#

"""
Nucleic acid chains
"""

__docformat__ = 'epytext'

from MMTK import Biopolymers, Bonds, ChemicalObjects, Collections, \
                 ConfigIO, Database, Universe, Utility
from Scientific.Geometry import Vector

from MMTK.Biopolymers import defineNucleicAcidResidue

#
# Nucleotides are special groups
#
class Nucleotide(Biopolymers.Residue):

    """
    Nucleic acid residue

    Nucleotides are a special kind of group. Like any other
    group, they are defined in the chemical database. Each residue
    has two or three subgroups ('sugar' and 'base', plus 'phosphate'
    except for 5'-terminal residues) and is usually
    connected to other residues to form a nucleotide chain. The database
    contains three variants of each residue (5'-terminal, 3'-terminal,
    non-terminal).
    """

    def __init__(self, name = None, model = 'all'):
        """
        @param name: the name of the nucleotide in the chemical database. This
                     is the full name of the residue plus the suffix
                     "_5ter" or "_3ter" for the terminal variants.
        @type name: C{str}
        @param model: one of "all" (all-atom), "none" (no hydrogens),
                      "polar" (united-atom with only polar hydrogens),
                      "polar_charmm" (like "polar", but defining
                      polar hydrogens like in the CHARMM force field).
                      Currently the database has definitions only for "all".
        @type model: C{str}
        """
	if name is not None:
	    blueprint = _residueBlueprint(name, model)
	    ChemicalObjects.Group.__init__(self, blueprint)
	    self.model = model
	    self._init()

    def backbone(self):
        """
        @returns: the sugar and phosphate groups
        @rtype: L{MMTK.ChemicalObjects.Group}
        """
        bb = self.sugar
        if hasattr(self, 'phosphate'):
            bb = Collections.Collection([bb, self.phosphate])
        return bb

    def bases(self):
        """
        @returns: the base group
        @rtype: L{MMTK.ChemicalObjects.Group}
        """
        return self.base


def _residueBlueprint(name, model):
    try:
	blueprint = _residue_blueprints[(name, model)]
    except KeyError:
	if model == 'polar':
	    name = name + '_uni'
	elif model == 'polar_charmm':
	    name = name + '_uni2'
	elif model == 'none':
	    name = name + '_noh'
	blueprint = Database.BlueprintGroup(name)
	_residue_blueprints[(name, model)] = blueprint
    return blueprint

_residue_blueprints = {}

#
# Nucleotide chains are molecules with added features.
#
class NucleotideChain(Biopolymers.ResidueChain):

    """
    Nucleotide chain

    Nucleotide chains consist of nucleotides that are linked together.
    They are a special kind of molecule, i.e. all molecule operations
    are available.

    Nucleotide chains act as sequences of residues. If C{n} is a NucleotideChain
    object, then

     - C{len(n)} yields the number of nucleotides

     - C{n[i]} yields nucleotide number C{i}

     - C{n[i:j]} yields the subchain from nucleotide number C{i} up to but
                 excluding nucleotide number C{j}
    """

    def __init__(self, sequence, **properties):
        """
        @param sequence: the nucleotide sequence. This can be a string
                         containing the one-letter codes, or a list
                         of two-letter codes (a "d" or "r" for the type of
                         sugar, and the one-letter base code), or a
                         L{MMTK.PDB.PDBNucleotideChain} object.
                         If a PDBNucleotideChain object is supplied, the
                         atomic positions it contains are assigned to the
                         atoms of the newly generated nucleotide chain,
                         otherwise the positions of all atoms are undefined.
        @keyword model: one of "all" (all-atom), "no_hydrogens" or "none"
                        (no hydrogens), "polar_hydrogens" or "polar"
                        (united-atom with only polar hydrogens),
                        "polar_charmm" (like "polar", but defining
                        polar hydrogens like in the CHARMM force field),
                        "polar_opls" (like "polar", but defining
                        polar hydrogens like in the latest OPLS force field).
                        Default is "all". Currently the database contains
                        definitions only for "all".
        @type model: C{str}
        @keyword terminus_5: if C{True}, the first residue is constructed
                             using the 5'-terminal variant, if C{False} the
                             non-terminal version is used. Default is C{True}.
        @type terminus_5: C{bool}
        @keyword terminus_3: if C{True}, the last residue is constructed
                             using the 3'-terminal variant, if C{False} the
                             non-terminal version is used. Default is C{True}.
        @type terminus_3: C{bool}
        @keyword circular: if C{True}, a bond is constructed
                           between the first and the last residue.
                           Default is C{False}.
        @type circular: C{bool}
        @keyword name: a name for the chain (a string)
        @type name: C{str}
        """
	if sequence is not None:
	    hydrogens = self.binaryProperty(properties, 'hydrogens', 'all')
	    if hydrogens == 1:
		hydrogens = 'all'
	    elif hydrogens == 0:
		hydrogens = 'none'
	    term5 = self.binaryProperty(properties, 'terminus_5', True)
	    term3 = self.binaryProperty(properties, 'terminus_3', True)
	    circular = self.binaryProperty(properties, 'circular', False)
	    try:
		model = properties['model'].lower()
	    except KeyError:
		model = hydrogens
	    self.version_spec = {'hydrogens': hydrogens,
				 'terminus_5': term5,
				 'terminus_3': term3,
				 'model': model,
                                 'circular': circular}
            if isinstance(sequence[0], str):
		conf = None
	    else:
		conf = sequence
		sequence = [r.name for r in sequence]
	    sequence = [Biopolymers._fullName(r) for r in sequence]
            if term5:
                if sequence[0].find('5ter') == -1:
                    sequence[0] += '_5ter'
            if term3:
                if sequence[-1].find('3ter') == -1:
                    sequence[-1] += '_3ter'

            self.groups = []
            n = 0
            for residue in sequence:
                n += 1
                r = Nucleotide(residue, model)
                r.name = r.symbol + '_' + `n`
                r.sequence_number = n
                r.parent = self
                self.groups.append(r)

            self._setupChain(circular, properties, conf)

    is_nucleotide_chain = True

    def __getslice__(self, first, last):
	return NucleotideSubChain(self, self.groups[first:last])

    def backbone(self):
        """
        @returns: the sugar and phosphate groups of all nucleotides
        @rtype: L{MMTK.Collections.Collection}
        """
        bb = Collections.Collection([])
        for residue in self.groups:
            try:
                bb.addObject(residue.phosphate)
            except AttributeError:
                pass
            bb.addObject(residue.sugar)
        return bb

    def bases(self):
        """
        @returns: the base groups of all nucleotides
        @rtype: L{MMTK.Collections.Collection}
        """
	return Collections.Collection([r.base for r in self.groups])

    def _descriptionSpec(self):
        kwargs = ','.join([name + '=' + `value`
                           for name, value in self.version_spec.items()])
	return "N", kwargs

    def _typeName(self):
        return ''.join([s.ljust(3) for s in self.sequence()])

    def _graphics(self, conf, distance_fn, model, module, options):
	if model != 'backbone':
	    return ChemicalObjects.Molecule._graphics(self, conf,
						      distance_fn, model,
						      module, options)
	color = options.get('color', 'black')
	material = module.EmissiveMaterial(color)
	objects = []
	for i in range(1, len(self.groups)-1):
	    a1 = self.groups[i].phosphate.P
	    a2 = self.groups[i+1].phosphate.P
	    p1 = a1.position(conf)
	    p2 = a2.position(conf)
	    if p1 is not None and p2 is not None:
		bond_vector = 0.5*distance_fn(a1, a2, conf)
		cut = bond_vector != 0.5*(p2-p1)
		if not cut:
		    objects.append(module.Line(p1, p2, material = material))
		else:
		    objects.append(module.Line(p1, p1+bond_vector,
					       material = material))
		    objects.append(module.Line(p2, p2-bond_vector,
					       material = material))
	return objects

#
# Subchains are created by slicing chains or extracting a chain from
# a group of connected chains.
#
class NucleotideSubChain(NucleotideChain):

    """
    A contiguous part of a nucleotide chain

    NucleotideSubChain objects are the result of slicing operations on
    NucleotideChain objects. They cannot be created directly.
    NucleotideSubChain objects permit all operations of NucleotideChain
    objects, but cannot be added to a universe.
    """

    def __init__(self, chain, groups, name = ''):
	self.groups = groups
	self.atoms = []
	self.bonds = []
	for g in self.groups:
	    self.atoms = self.atoms + g.atoms
	    self.bonds = self.bonds + g.bonds
        for i in range(len(self.groups)-1):
            self.bonds.append(Bonds.Bond((self.groups[i].sugar.O_3,
                                          self.groups[i+1].phosphate.P)))
	self.bonds = Bonds.BondList(self.bonds)
	self.name = name
	self.parent = chain.parent
	self.type = None
	self.configurations = {}
	self.part_of = chain

    is_incomplete = True

    def __repr__(self):
	if self.name == '':
	    return 'SubChain of ' + repr(self.part_of)
	else:
	    return ChemicalObjects.Molecule.__repr__(self)
    __str__ = __repr__


#
# Type check functions
#
def isNucleotideChain(x):
    "Returns C{True} if x is a NucleotideChain."
    return hasattr(x, 'is_nucleotide_chain')
