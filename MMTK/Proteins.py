# This module implements classes for peptide chains and proteins.
#
# Written by Konrad Hinsen
#

"""
Peptide chains and proteins
"""

__docformat__ = 'restructuredtext'

from MMTK import Biopolymers, Bonds, ChemicalObjects, Collections, \
                 ConfigIO, Database, Units, Universe, Utility
from Scientific.Geometry import Vector

from MMTK.Biopolymers import defineAminoAcidResidue

#
# Residues are special groups
#
class Residue(Biopolymers.Residue):

    """
    Amino acid residue

    Amino acid residues are a special kind of group. They are defined
    in the chemical database. Each residue has two subgroups
    ('peptide' and 'sidechain') and is usually connected to other
    residues to form a peptide chain. The database contains three
    variants of each residue (N-terminal, C-terminal,
    non-terminal) and various models (all-atom, united-atom,
    |C_alpha|).
    """

    def __init__(self, name = None, model = 'all'):
        """
        :param name: the name of the residue in the chemical database. This
                     is the full name of the residue plus the suffix
                     "_nt" or "_ct" for the terminal variants.
        :type name: str
        :param model: one of "all" (all-atom), "none" (no hydrogens),
                      "polar" (united-atom with only polar hydrogens),
                      "polar_charmm" (like "polar", but defining
                      polar hydrogens like in the CHARMM force field),
                      "polar_opls" (like "polar", but defining
                      polar hydrogens like in the latest OPLS force field),
                      "calpha" (only the |C_alpha| atom).
        :type model: str
        """
        if name is not None:
            blueprint = _residueBlueprint(name, model)
            ChemicalObjects.Group.__init__(self, blueprint)
            self.model = model
            self._init()

    def _init(self):
        Biopolymers.Residue._init(self)
        # create peptide attribute for calpha model
        if self.model == 'calpha':
            self.peptide = self

    def isNTerminus(self):
        return hasattr(self.peptide, 'H_3')

    def isCTerminus(self):
        return hasattr(self.peptide, 'O_2')

    def _makeCystine(self):
        if self.model == 'calpha':
            return self
        if self.symbol.lower() != 'cys':
            raise ValueError(`self` + " is not cysteine.")
        new_residue = 'cystine_ss'
        if self.isNTerminus():
            new_residue = new_residue + '_nt'
        elif self.isCTerminus():
            new_residue = new_residue + '_ct'
        new_residue = Residue(new_residue, self.model)
        for g in ['peptide', 'sidechain']:
            g_old = getattr(self, g)
            g_new = getattr(new_residue, g)
            for a in getattr(g_new, 'atoms'):
                set_method = getattr(getattr(g_new, a.name), 'setPosition')
                set_method(getattr(getattr(g_old, a.name), 'position')())
        return new_residue

    def isSubsetModel(self):
        return self.model == 'calpha'

    def backbone(self):
        """
        :returns: the peptide group
        :rtype: :class:`~MMTK.ChemicalObjects.Group`
        """
        return self.peptide

    def sidechains(self):
        """
        :returns: the sidechain group
        :rtype: :class:`~MMTK.ChemicalObjects.Group`
        """
        return self.sidechain

    def phiPsi(self, conf = None):
        """
        :returns: the values of the backbone dihedral angles phi and psi.
        :rtype: tuple (float, float)
        """
        universe = self.universe()
        if universe is None:
            universe = Universe.InfiniteUniverse()
        C = None
        for a in self.peptide.N.bondedTo():
            if a.parent.parent != self:
                C = a
                break
        if C is None:
            phi = None
        else:
            phi = universe.dihedral(self.peptide.C, self.peptide.C_alpha,
                                    self.peptide.N, C, conf)
        N = None
        for a in self.peptide.C.bondedTo():
            if a.parent.parent != self:
                N = a
                break
        if N is None:
            psi = None
        else:
            psi = universe.dihedral(N, self.peptide.C, self.peptide.C_alpha,
                                    self.peptide.N, conf)
        return phi, psi

    def phiAngle(self):
        """
        :returns: an object representing the phi angle and allowing to modify it
        :rtype: MMTK.InternalCoordinates.DihedralAngle
        """
        from MMTK.InternalCoordinates import DihedralAngle
        C = None
        for a in self.peptide.N.bondedTo():
            if a.parent.parent != self:
                C = a
                break
        if C is None:
            raise ValueError("residue is N-terminus")
        return DihedralAngle(self.peptide.C, self.peptide.C_alpha,
                             self.peptide.N, C)

    def psiAngle(self):
        """
        :returns: an object representing the psi angle and allowing to modify it
        :rtype: MMTK.InternalCoordinates.DihedralAngle
        """
        from MMTK.InternalCoordinates import DihedralAngle
        N = None
        for a in self.peptide.C.bondedTo():
            if a.parent.parent != self:
                N = a
                break
        if N is None:
            raise ValueError("residue is C-terminus")
        return DihedralAngle(N, self.peptide.C, self.peptide.C_alpha,
                             self.peptide.N)

    def chiAngle(self):
        """
        :returns: an object representing the chi angle and allowing to modify it
        :rtype: MMTK.InternalCoordinates.DihedralAngle
        """
        from MMTK.InternalCoordinates import DihedralAngle
        try:
            C_beta = self.sidechain.C_beta
        except AttributeError:
            raise ValueError("no C_beta in sidechain")
        X = None
        for atom_name in ['C_gamma', 'C_gamma_1', 'S_gamma',
                          'O_gamma', 'O_gamma_1', 'H_beta_1']:
            try:
                X = getattr(self.sidechain, atom_name)
                break
            except AttributeError:
                pass
        if X is None:
            raise ValueError("no sidechain reference atom found")
        return DihedralAngle(self.peptide.N, self.peptide.C_alpha,
                             C_beta, X)


def _residueBlueprint(name, model):
    try:
        blueprint = _residue_blueprints[(name, model)]
    except KeyError:
        if model == 'polar':
            name = name + '_uni'
        elif model == 'polar_charmm':
            name = name + '_uni2'
        elif model == 'polar_oldopls':
            name = name + '_uni3'
        elif model == 'none':
            name = name + '_noh'
        elif model == 'calpha':
            name = name + '_calpha'
        blueprint = Database.BlueprintGroup(name)
        _residue_blueprints[(name, model)] = blueprint
    return blueprint

_residue_blueprints = {}

#
# Peptide chains are molecules with added features.
#
class PeptideChain(Biopolymers.ResidueChain):

    """
    Peptide chain

    Peptide chains consist of amino acid residues that are linked
    by peptide bonds. They are a special kind of molecule, i.e.
    all molecule operations are available.

    Peptide chains act as sequences of residues. If p is a PeptideChain
    object, then

     * len(p) yields the number of residues
     * p[i] yields residue number i
     * p[i:j] yields the subchain from residue number i up to
                 but excluding residue number j

    :param sequence: the amino acid sequence. This can be a string
                     containing the one-letter codes, or a list
                     of three-letter codes, or a
                     :class:`~MMTK.PDB.PDBPeptideChain` object.
                     If a PDBPeptideChain object is supplied, the atomic
                     positions it contains are assigned to the atoms
                     of the newly generated peptide chain, otherwise the
                     positions of all atoms are undefined.
    :keyword model: one of "all" (all-atom), "no_hydrogens" or "none"
                    (no hydrogens), "polar_hydrogens" or "polar"
                    (united-atom with only polar hydrogens),
                    "polar_charmm" (like "polar", but defining
                    polar hydrogens like in the CHARMM force field),
                    "polar_opls" (like "polar", but defining
                    polar hydrogens like in the latest OPLS force field),
                    "calpha" (only the |C_alpha| atom of each residue).
                    Default is "all".
    :type model: str
    :keyword n_terminus: if True, the first residue is constructed
                         using the N-terminal variant, if False the
                         non-terminal version is used. Default is True.
    :type n_terminus: bool
    :keyword c_terminus: if True, the last residue is constructed
                         using the C-terminal variant, if False the
                         non-terminal version is used. Default is True.
    :type c_terminus: bool
    :keyword circular: if True, a peptide bond is constructed
                       between the first and the last residues.
                       Default is False.
    :type circular: bool
    :keyword name: a name for the chain (a string)
    :type name: str

    """

    def __init__(self, sequence, **properties):
        if sequence is not None:
            model = 'all'
            if properties.has_key('model'):
                model = properties['model'].lower()
            elif properties.has_key('hydrogens'):
                model = properties['hydrogens']
                if model == 1: model = 'all'
                elif model == 0: model = 'none'
                else: model = model.lower()
            if model == 'no_hydrogens':
                model = 'none'
            elif model == 'polar_hydrogens':
                model = 'polar'
            n_term = self.binaryProperty(properties, 'n_terminus', True)
            c_term = self.binaryProperty(properties, 'c_terminus', True)
            circular = self.binaryProperty(properties, 'circular', False)
            self.version_spec = {'n_terminus': n_term,
                                 'c_terminus': c_term,
                                 'model': model,
                                 'circular': circular}
            if isinstance(sequence[0], basestring):
                conf = None
                numbers = range(len(sequence))
            else:
                conf = sequence
                sequence = conf.sequence()
                numbers = [r.number for r in conf]
            sequence = map(Biopolymers._fullName, sequence)
            if model != 'calpha':
                if n_term:
                    sequence[0] = sequence[0] + '_nt'
                if c_term:
                    sequence[-1] = sequence[-1] + '_ct'

            self.groups = []
            n = 0
            for residue, number in zip(sequence, numbers):
                n = n + 1
                r = Residue(residue, model)
                r.name = r.symbol + str(number)
                r.sequence_number = n
                r.parent = self
                self.groups.append(r)

            self._setupChain(circular, properties, conf)

    is_peptide_chain = True

    def __getslice__(self, first, last):
        return SubChain(self, self.groups[first:last])

    def sequence(self):
        """
        :returns: the primary sequence as a list of three-letter
                  residue codes.
        :rtype: list
        """
        return [r.symbol for r in self.groups]

    def backbone(self):
        """
        :returns: the peptide groups of all residues
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        backbone = Collections.Collection()
        for r in self.groups:
            try:
                backbone.addObject(r.peptide)
            except AttributeError:
                pass
        return backbone
    
    def sidechains(self):
        """
        :returns: the sidechain groups of all residues
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        sidechains = Collections.Collection()
        for r in self.groups:
            try:
                sidechains.addObject(r.sidechain)
            except AttributeError:
                pass
        return sidechains

    def phiPsi(self, conf = None):
        """
        :returns: a list of the (phi, psi) backbone angles for each residue
        :rtype: list of tuple of float
        """
        universe = self.universe()
        if universe is None:
            universe = Universe.InfiniteUniverse()
        angles = []
        for i in range(len(self)):
            r = self[i]
            if i == 0:
                phi = None
            else:
                phi = universe.dihedral(r.peptide.C, r.peptide.C_alpha,
                                        r.peptide.N,
                                        self[i-1].peptide.C, conf)
            if i == len(self)-1:
                psi = None
            else:
                psi = universe.dihedral(self[i+1].peptide.N,
                                        r.peptide.C, r.peptide.C_alpha,
                                        r.peptide.N, conf)
            angles.append((phi, psi))
        return angles

    def replaceResidue(self, r_old, r_new):
        """
        :param r_old: the residue to be replaced (must be part of the chain)
        :type r_old: Residue
        :param r_new: the residue that replaces r_old
        :type r_new: Residue
        """
        n = self.groups.index(r_old)
        for a in r_old.atoms:
            self.atoms.remove(a)
        obsolete_bonds = []
        for b in self.bonds:
            if b.a1 in r_old.atoms or b.a2 in r_old.atoms:
                obsolete_bonds.append(b)
        for b in obsolete_bonds:
            self.bonds.remove(b)
        r_old.parent = None
        self.atoms.extend(r_new.atoms)
        self.bonds.extend(r_new.bonds)
        r_new.sequence_number = n+1
        if r_old.name.startswith(r_old.symbol):
            r_new.name = r_new.symbol+r_old.name[len(r_old.symbol):]
        else:
            r_new.name = r_new.symbol+`n+1`
        r_new.parent = self
        self.groups[n] = r_new
        if n > 0:
            peptide_old = self.bonds.bondsOf(r_old.peptide.N)
            if peptide_old:
                self.bonds.remove(peptide_old[0])
            if not (self.groups[n-1].isCTerminus()
                    or self.groups[n].isNTerminus()):
                # ConnectedChain objects can have N/C-terminal
                # residues inside the (virtual) chain, so the
                # test is necessary.
                self.bonds.append(Bonds.Bond((self.groups[n-1].peptide.C,
                                              self.groups[n].peptide.N)))
        if n < len(self.groups)-1:
            peptide_old = self.bonds.bondsOf(r_old.peptide.C)
            if peptide_old:
                self.bonds.remove(peptide_old[0])
            if not (self.groups[n].isCTerminus()
                    or self.groups[n+1].isNTerminus()):
                self.bonds.append(Bonds.Bond((self.groups[n].peptide.C,
                                              self.groups[n+1].peptide.N)))
        if isinstance(self.parent, ChemicalObjects.Complex):
            self.parent.recreateAtomList()
        universe = self.universe()
        if universe is not None:
            universe._changed(True)

    # add sulfur bridges between cysteine residues
    def _addSSBridges(self, bonds):
        for b in bonds:
            cys1 = b[0]
            if cys1.symbol.lower() == 'cyx':
                cys_ss1 = cys1
            else:
                cys_ss1 = cys1._makeCystine()
                self.replaceResidue(cys1, cys_ss1)
            cys2 = b[1]
            if cys2.symbol.lower() == 'cyx':
                cys_ss2 = cys2
            else:
                cys_ss2 = cys2._makeCystine()
                self.replaceResidue(cys2, cys_ss2)
            self.bonds.append(Bonds.Bond((cys_ss1.sidechain.S_gamma,
                                          cys_ss2.sidechain.S_gamma)))

    def _descriptionSpec(self):
        kwargs = ','.join([name + '=' + `self.version_spec[name]`
                           for name in sorted(self.version_spec.keys())])
	return "S", kwargs

    def _typeName(self):
        return ''.join(self.sequence())

    def _graphics(self, conf, distance_fn, model, module, options):
        if model != 'backbone':
            return ChemicalObjects.Molecule._graphics(self, conf,
                                                      distance_fn, model,
                                                      module, options)
        color = options.get('color', 'black')
        material = module.EmissiveMaterial(color)
        objects = []
        for i in range(len(self.groups)-1):
            a1 = self.groups[i].peptide.C_alpha
            a2 = self.groups[i+1].peptide.C_alpha
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
class SubChain(PeptideChain):

    """
    A contiguous part of a peptide chain

    SubChain objects are the result of slicing operations on
    PeptideChain objects. They cannot be created directly.
    SubChain objects permit all operations of PeptideChain
    objects, but cannot be added to a universe.
    """

    def __init__(self, chain=None, groups=None, name = ''):
        if chain is not None:
            self.groups = groups
            self.atoms = []
            self.bonds = []
            for g in self.groups:
                self.atoms.extend(g.atoms)
                self.bonds.extend(g.bonds)
            for i in range(len(self.groups)-1):
                link1 = self.groups[i].chain_links[1]
                link2 = self.groups[i+1].chain_links[0]
                self.bonds.append(Bonds.Bond((link1, link2)))
            self.bonds = Bonds.BondList(self.bonds)
            self.name = name
            self.model = chain.model
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

    def replaceResidue(self, r_old, r_new):
        for a in r_old.atoms:
            self.atoms.remove(a)
        obsolete_bonds = []
        for b in self.bonds:
            if b.a1 in r_old.atoms or b.a2 in r_old.atoms:
                obsolete_bonds.append(b)
        for b in obsolete_bonds:
            self.bonds.remove(b)
        n = self.groups.index(r_old)
        if n > 0:
            for b in self.bonds.bondsOf(r_old.peptide.N):
                self.bonds.remove(b)
        if n < len(self.groups)-1:
            for b in self.bonds.bondsOf(r_old.peptide.C):
                self.bonds.remove(b)
        PeptideChain.replaceResidue(self.part_of, r_old, r_new)
        self.groups[n] = r_new
        self.atoms.extend(r_new.atoms)
        self.bonds.extend(r_new.bonds)
        if n > 0:
            self.bonds.append(Bonds.Bond((self.groups[n-1].peptide.C,
                                          self.groups[n].peptide.N)))
        if n < len(self.groups)-1:
            self.bonds.append(Bonds.Bond((self.groups[n].peptide.C,
                                          self.groups[n+1].peptide.N)))

    def _distanceConstraintList(self):
        atoms = self.atomList()
        return [(a1, a2, d)
                for a1, a2, d in self.part_of._distanceConstraintList()
                if a1 in atoms and a2 in atoms]

    def addDistanceConstraint(self, atom1, atom2, distance):
        chain = self
        while True:
            try:
                chain = chain.part_of
            except AttributeError:
                break
        try:
            chain.distance_constraints.append((atom1, atom2, distance))
        except AttributeError:
            chain.distance_constraints = [(atom1, atom2, distance)]

    def removeDistanceConstraints(self, universe=None):
        raise NotImplementedError

#
# Connected chains are collections of peptide chains connected by s-s bridges.
#
class ConnectedChains(PeptideChain):

    """
    Peptide chains connected by disulfide bridges
    
    A group of peptide chains connected by disulfide bridges must be considered
    a single molecule due to the presence of chemical bonds. Such a molecule
    is represented by a ConnectedChains object. These objects are created
    automatically when a Protein object is assembled. They are normally
    not used directly by application programs. When a chain with disulfide
    bridges to other chains is extracted from a Protein object, the
    return value is a SubChain object that indirectly refers to a
    ConnectedChains object.
    """

    def __init__(self, chains=None):
        if chains is not None:
            self.chains = []
            self.groups = []
            self.atoms = []
            self.bonds = Bonds.BondList([])
            self.chain_names = []
            self.model = chains[0].model
            version_spec = chains[0].version_spec
            for c in chains:
                if c.version_spec['model'] != version_spec['model']:
                    raise ValueError("mixing chains of different model: " +
                                      c.version_spec['model'] + "/" +
                                      version_spec['model'])
                ng = len(self.groups)
                self.chains.append((c.name, ng, ng+len(c.groups),
                                    c.version_spec))
                self.groups.extend(c.groups)
                self.atoms.extend(c.atoms)
                self.bonds.extend(c.bonds)
                try: name = c.name
                except AttributeError: name = ''
                self.chain_names.append(name)
            for g in self.groups:
                g.parent = self
            self.name = ''
            self.parent = None
            self.type = None
            self.configurations = {}
    is_connected_chains = True

    def _finalize(self):
        for i in range(len(self.chains)):
            c = self.chains[i]
            sub_chain = SubChain(self, self.groups[c[1]:c[2]], c[0])
            sub_chain.version_spec = c[3]
            for g in sub_chain.groups:
                g.parent = sub_chain
            self.chains[i] = sub_chain

    def __len__(self):
        return len(self.chains)

    def __getitem__(self, item):
        return self.chains[item]

    def __getslice__(self, first, last):
        raise TypeError("Can't slice connected chains")

    def _graphics(self, conf, distance_fn, model, module, options):
        if model != 'backbone':
            return ChemicalObjects.Molecule._graphics(self, conf,
                                                      distance_fn, model,
                                                      module, options)
        objects = []
        for chain in self:
            objects = objects + chain._graphics(conf, distance_fn,
                                                model, module, options)
        return objects

#
# Proteins are complexes of peptide chains, connected peptide chains,
# and possibly other things.
#
class Protein(ChemicalObjects.Complex):

    """
    Protein

    A Protein object is a special kind of :class:`~MMTK.ChemicalObjects.Complex`
    object which is made up of peptide chains and possibly ligands.

    If the atoms in the peptide chains that make up a protein have
    defined positions, sulfur bridges within chains and between
    chains will be constructed automatically during protein generation
    based on a distance criterion between cystein sidechains.


    Proteins act as sequences of chains. If p is a Protein object, then

    * len(p) yields the number of chains
    * p[i] yields chain number i

    """

    def __init__(self, *items, **properties):
        """
        :param items: either a sequence of peptide chain objects, or
                      a string, which is interpreted as the name of a
                      database definition for a protein.
                      If that definition does not exist, the string
                      is taken to be the name of a PDB file, from which
                      all peptide chains are constructed and
                      assembled into a protein.
        :keyword model: one of "all" (all-atom), "no_hydrogens" or "none"
                        (no hydrogens),"polar_hydrogens" or "polar"
                        (united-atom with only polar hydrogens),
                        "polar_charmm" (like "polar", but defining
                        polar hydrogens like in the CHARMM force field),
                        "polar_opls" (like "polar", but defining
                        polar hydrogens like in the latest OPLS force field),
                        "calpha" (only the |C_alpha| atom of each residue).
                        Default is "all".
        :type model: str
        :keyword position: the center-of-mass position of the protein
        :type position: Scientific.Geometry.Vector
        :keyword name: a name for the protein
        :type name: str
        """
        if items == (None,):
            return
        self.name = ''
        if len(items) == 1 and type(items[0]) == type(''):
            try:
                filename = Database.databasePath(items[0], 'Proteins')
                found = 1
            except IOError:
                found = 0
            if found:
                blueprint = Database.BlueprintProtein(items[0])
                items = blueprint.chains
                for attr, value in vars(blueprint).items():
                    if attr not in ['type', 'chains']:
                        setattr(self, attr, value)
            else:
                import PDB
                conf = PDB.PDBConfiguration(items[0])
                model = properties.get('model', 'all')
                items = conf.createPeptideChains(model)
        molecules = []
        for i in items:
            if ChemicalObjects.isChemicalObject(i):
                molecules.append(i)
            else:
                molecules = molecules + list(i)
        for m, i in zip(molecules, range(len(molecules))):
            m._numbers = [i]
            if not m.name:
                m.name = 'chain'+`i`
        ss = self._findSSBridges(molecules)
        new_mol = {}
        for m in molecules:
            new_mol[m] = ([m],[])
        for bond in ss:
            m1 = new_mol[bond[0].topLevelChemicalObject()]
            m2 = new_mol[bond[1].topLevelChemicalObject()]
            if m1 == m2:
                m1[1].append(bond)
            else:
                combined = (m1[0] + m2[0], m1[1] + m2[1] + [bond])
                for m in combined[0]:
                    new_mol[m] = combined
        self.molecules = []
        while new_mol:
            m = new_mol.values()[0]
            for i in m[0]:
                del new_mol[i]
            bonds = m[1]
            if len(m[0]) == 1:
                m = m[0][0]
                m._addSSBridges(bonds)
            else:
                numbers = sum((i._numbers for i in m[0]), [])
                m = ConnectedChains(m[0])
                m._numbers = numbers
                m._addSSBridges(bonds)
                m._finalize()
                for c in m:
                    c.parent = self
            m.parent = self
            self.molecules.append(m)

        self.atoms = []
        self.chains = []
        for m in self.molecules:
            self.atoms.extend(m.atoms)
            if hasattr(m, 'is_connected_chains'):
                for c, name, i in zip(range(len(m)),
                                   m.chain_names, m._numbers):
                    self.chains.append((m, c, name, i))
            else:
                try: name = m.name
                except AttributeError: name = ''
                self.chains.append((m, None, name, m._numbers[0]))
        self.chains.sort(lambda c1, c2: cmp(c1[3], c2[3]))
        self.chains = map(lambda c: c[:3], self.chains)

        self.parent = None
        self.type = None
        self.configurations = {}
        try:
            self.name = properties['name']
            del properties['name']
        except KeyError: pass
        if properties.has_key('position'):
            self.translateTo(properties['position'])
            del properties['position']
        self.addProperties(properties)

        undefined = 0
        for a in self.atoms:
            if a.position() is None:
                undefined += 1
        if undefined > 0 and undefined != len(self.atoms):
            Utility.warning('Some atoms in a protein ' +
                            'have undefined positions.')

    is_protein = True

    def __len__(self):
        return len(self.chains)

    def __getitem__(self, item):
        if isinstance(item, int):
            m, c, name = self.chains[item]
        else:
            for m, c, name in self.chains:
                if name == item:
                    break
            if name != item:
                raise ValueError('No chain with name ' + item)
        if c is None:
            return m
        else:
            return m[c]

    def residuesOfType(self, *types):
        """
        :param types: a sequence of residue codes (one- or three-letter)
        :type types: sequence of str
        :returns: all residues whose type (one- or three-letter code)
                  is contained in types
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        rlist = Collections.Collection([])
        for m in self.molecules:
            if isPeptideChain(m):
                rlist = rlist + apply(m.residuesOfType, types)
        return rlist

    def backbone(self):
        """
        :returns: the peptide groups of all residues in all chains
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        rlist = Collections.Collection([])
        for m in self.molecules:
            if isPeptideChain(m):
                rlist = rlist + m.backbone()
        return rlist

    def sidechains(self):
        """
        :returns: the sidechain groups of all residues in all chains
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        rlist = Collections.Collection([])
        for m in self.molecules:
            if isPeptideChain(m):
                rlist = rlist + m.sidechains()
        return rlist

    def residues(self):
        """
        :returns: all residues in all chains
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        rlist = Collections.Collection([])
        for m in self.molecules:
            if isPeptideChain(m):
                rlist = rlist + m.residues()
        return rlist

    def phiPsi(self, conf = None):
        """
        :returns: a list of the (phi, psi) backbone angles for all residue
                  in all chains
        :rtype: list of list of tuple of float
        """
        return [chain.phiPsi(conf) for chain in self]

    _ss_bond_max = 0.25*Units.nm

    def _findSSBridges(self, molecules):
        molecules = filter(lambda m: hasattr(m, 'is_peptide_chain'), molecules)
        cys = Collections.Collection([])
        for m in molecules:
            if m.version_spec['model'] != 'calpha':
                cys = cys + m.residuesOfType('cys') + m.residuesOfType('cyx')
        s = cys.map(lambda r: r.sidechain.S_gamma)
        ns = len(s)
        ss = []
        for i in xrange(ns-1):
            for j in xrange(i+1,ns):
                r1 = s[i].position()
                r2 = s[j].position()
                if r1 and r2 and (r1-r2).length() < self._ss_bond_max:
                    ss.append((cys[i], cys[j]))
        return ss

    def _subunits(self):
        return list(self)

    def _description(self, tag, index_map, toplevel):
        if not toplevel:
            raise ValueError
        return 'l(' + `self.__class__.__name__` + ',' + `self.name` + ',[' + \
               ','.join(o._description(tag, index_map, True) for o in self) + \
               '])'

    def _graphics(self, conf, distance_fn, model, module, options):
        if model != 'backbone':
            return ChemicalObjects.Complex._graphics(self, conf, distance_fn,
                                                     model, module, options)
        objects = []
        for chain in self:
            objects.extend(chain._graphics(conf, distance_fn,
                                           model, module, options))
        return objects

#
# Type check functions
#
def isPeptideChain(x):
    """
    :param x: any object
    :returns: True if x is a peptide chain
    :rtype: bool
    """
    return hasattr(x, 'is_peptide_chain')

def isProtein(x):
    """
    :param x: any object
    :returns: True if x is a protein
    :rtype: bool
    """
    return hasattr(x, 'is_protein')
