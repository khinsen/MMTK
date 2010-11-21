# This module deals with input and output of configurations in PDB format.
#
# Written by Konrad Hinsen
#

"""
PDB files

This module provides classes that represent molecules in a PDB file.
They permit various manipulations and the creation of MMTK objects.
Note that the classes defined in this module are specializations
of classed defined in Scientific.IO.PDB; the methods defined in
that module are also available.
"""

__docformat__ = 'restructuredtext'

from MMTK import ChemicalObjects, Collections, Database, Units, \
                 Universe, Utility
from Scientific.Geometry import Vector
import Scientific.IO.PDB
import copy, string

#
# The chain classes from Scientific.IO.PDB are extended by methods
# that construct MMTK objects and set configurations.
#
class PDBChain(object):

    def applyTo(self, chain, map = 'pdbmap', alt = 'pdb_alternative',
                atom_map = None):
        if len(chain) != len(self):
            raise ValueError("chain lengths do not match")
        for i in range(len(chain)):
            residue = chain[i]
            pdbmap = getattr(residue, map)
            try: altmap = getattr(residue, alt)
            except AttributeError: altmap = {}
            setResidueConfiguration(residue, self[i], pdbmap[0], altmap,
                                    atom_map)
    
class PDBPeptideChain(Scientific.IO.PDB.PeptideChain, PDBChain):

    """
    Peptide chain in a PDB file

    See the description of Scientific.IO.PDB.PeptideChain for the
    constructor and additional methods. In MMTK, PDBPeptideChain
    objects are usually obtained from a PDBConfiguration object via
    its attribute peptide_chains (see the documentation of
    Scientific.IO.PDB.Structure).
    """

    def createPeptideChain(self, model = 'all',
                           n_terminus=None, c_terminus=None):
        """
        :returns: a :class:`~MMTK.Proteins.PeptideChain` object corresponding to the
                  peptide chain in the PDB file. The parameter model
                  has the same meaning as for the PeptideChain constructor.
        :rtype: :class:`~MMTK.Proteins.PeptideChain`
        """
        self.identifyProtonation()
        from MMTK import Proteins
        properties = {'model': model}
        if self.segment_id != '':
            properties['name'] = self.segment_id
        elif self.chain_id != '':
            properties['name'] = self.chain_id
        if c_terminus is None:
            properties['c_terminus'] = self.isTerminated()
        else:
            properties['c_terminus'] = c_terminus
        if n_terminus is not None:
            properties['n_terminus'] = n_terminus
        chain = apply(Proteins.PeptideChain, (self,), properties)
        if model != 'no_hydrogens':
            chain.findHydrogenPositions()
        return chain

    def identifyProtonation(self):
        for residue in self.residues:
            if residue.name == 'HIS':
                count_hd = 0
                count_he = 0
                for atom in residue:
                    if 'HD' in atom.name:
                        count_hd += 1
                    if 'HE' in atom.name:
                        count_he += 1
                if count_hd + count_he == 0:
                    # default for crystallographic structures without hydrogens
                    residue.name = 'HIE'
                elif count_he == 2:
                    if count_hd == 2:
                        residue.name = 'HIP'
                    else:
                        residue.name = 'HIE'
                else:
                    residue.name = 'HID'
            elif residue.name == 'GLU':
                for atom in residue:
                    if 'HE' in atom.name:
                        residue.name = 'GLP'
                        break
            elif residue.name == 'ASP':
                for atom in residue:
                    if 'HD' in atom.name:
                        residue.name = 'APP'
                        break
            elif residue.name == 'LYS':
                count_hz = 0
                for atom in residue:
                    if 'HZ' in atom.name:
                        count_hz += 1
                if count_hz > 0 and count_hz < 3:
                    # for count_hz == 0 (most probably a crystallographic
                    # structure), keep LYS which is the most frequent one.
                    residue.name = 'LYP'

class PDBNucleotideChain(Scientific.IO.PDB.NucleotideChain, PDBChain):

    """
    Nucleotide chain in a PDB file

    See the description of Scientific.IO.PDB.NucleotideChain for the
    constructor and additional methods. In MMTK, PDBNucleotideChain
    objects are usually obtained from a PDBConfiguration object via
    its attribute nucleotide_chains (see the documentation of
    Scientific.IO.PDB.Structure).
    """

    def createNucleotideChain(self, model='all'):
        """
        :returns: a :class:`~MMTK.NucleicAcids.NucleotideChain` object corresponding 
                  to the nucleotide chain in the PDB file. The parameter model
                  has the same meaning as for the NucleotideChain constructor.
        :rtype: :class:`~MMTK.NucleicAcids.NucleotideChain`
        """
        from MMTK import NucleicAcids
        properties = {'model': model}
        if self.segment_id != '':
            properties['name'] = self.segment_id
        if self[0].hasPhosphate():
            properties['terminus_5'] = 0
        chain = apply(NucleicAcids.NucleotideChain, (self,), properties)
        if model != 'no_hydrogens':
            chain.findHydrogenPositions()
        return chain

class PDBMolecule(Scientific.IO.PDB.Molecule):

    """
    Molecule in a PDB file

    See the description of Scientific.IO.PDB.Molecule for the
    constructor and additional methods. In MMTK, PDBMolecule objects
    are usually obtained from a PDBConfiguration object via its
    attribute molecules (see the documentation of
    Scientific.IO.PDB.Structure). A molecule is by definition any
    residue in a PDB file that is not an amino acid or nucleotide
    residue.
    """

    def applyTo(self, molecule, map = 'pdbmap', alt = 'pdb_alternative',
                atom_map = None):
        pdbmap = getattr(molecule, map)
        try: altmap = getattr(molecule, alt)
        except AttributeError: altmap = {}
        setResidueConfiguration(molecule, self, pdbmap[0], altmap, atom_map)

    def createMolecule(self, name=None):
        """
        :returns: a :class:`~MMTK.ChemicalObjects.Molecule` object corresponding 
                  to the molecule in the PDB file. The parameter name 
                  specifies the molecule name as defined in the chemical database.
                  It can be left out for known molecules (currently
                  only water).
        :rtype: :class:`~MMTK.ChemicalObjects.Molecule`
        """
        if name is None:
            name = molecule_names[self.name]
        m = ChemicalObjects.Molecule(name)
        setConfiguration(m, [self])
        return m

#
# The structure class from Scientific.IO.PDB is extended by methods
# that construct MMTK objects and set configurations.
#
class PDBConfiguration(Scientific.IO.PDB.Structure):

    """
    Everything in a PDB file

    A PDBConfiguration object represents the full contents of a PDB
    file. It can be used to create MMTK objects for all or part
    of the molecules, or to change the configuration of an existing
    system.

    :param file_or_filename: the name of a PDB file, or a file object
    :param model: the number of the model to be used from a
                  multiple model file
    :type model: int
    :param alternate_code: the alternate code to be used for atoms that
                           have multiple positions
    :type alternate_code: str
    """
    
    def __init__(self, file_or_filename, model = 0, alternate_code = 'A'):
        if isinstance(file_or_filename, str):
            file_or_filename = Database.PDBPath(file_or_filename)
        Scientific.IO.PDB.Structure.__init__(self, file_or_filename,
                                             model, alternate_code)
        self._numberAtoms()
        self._convertUnits()

    peptide_chain_constructor = PDBPeptideChain
    nucleotide_chain_constructor = PDBNucleotideChain
    molecule_constructor = PDBMolecule

    def _numberAtoms(self):
        n = 1
        for residue in self.residues:
            for atom in residue:
                atom.number = n
                n += 1

    def _convertUnits(self):
        for residue in self.residues:
            for atom in residue:
                atom.position = atom.position*Units.Ang
                try:
                    b = atom.properties['temperature_factor']
                    atom.properties['temperature_factor'] = b*Units.Ang**2
                except KeyError:
                    pass
                try:
                    u = atom.properties['u']
                    atom.properties['u'] = u*Units.Ang**2
                except KeyError:
                    pass

        # All these attributes exist only if ScientificPython >= 2.7.5 is used.
        # The Scaling transformation was introduced with the same version,
        # so if it exists, the rest should work as well.
        try:
            from Scientific.Geometry.Transformation import Scaling
        except ImportError:
            return
        for attribute in ['a', 'b', 'c']:
            value = getattr(self, attribute)
            if value is not None:
                setattr(self, attribute, value*Units.Ang)
        for attribute in ['alpha', 'beta', 'gamma']:
            value = getattr(self, attribute)
            if value is not None:
                setattr(self, attribute, value*Units.deg)
        if self.to_fractional is not None:
            self.to_fractional = self.to_fractional*Scaling(1./Units.Ang)
            v1 = self.to_fractional(Vector(1., 0., 0.))
            v2 = self.to_fractional(Vector(0., 1., 0.))
            v3 = self.to_fractional(Vector(0., 0., 1.))
            self.reciprocal_basis = (Vector(v1[0], v2[0], v3[0]),
                                     Vector(v1[1], v2[1], v3[1]),
                                     Vector(v1[2], v2[2], v3[2]))
        else:
            self.reciprocal_basis = None
        if self.from_fractional is not None:
            self.from_fractional = Scaling(Units.Ang)*self.from_fractional
            self.basis = (self.from_fractional(Vector(1., 0., 0.)),
                          self.from_fractional(Vector(0., 1., 0.)),
                          self.from_fractional(Vector(0., 0., 1.)))
        else:
            self.basis = None
        for i in range(len(self.ncs_transformations)):
            tr = self.ncs_transformations[i]
            tr_new = Scaling(Units.Ang)*tr*Scaling(1./Units.Ang)
            tr_new.given = tr.given
            tr_new.serial = tr.serial
            self.ncs_transformations[i] = tr_new
        for i in range(len(self.cs_transformations)):
            tr = self.cs_transformations[i]
            tr_new = Scaling(Units.Ang)*tr*Scaling(1./Units.Ang)
            self.cs_transformations[i] = tr_new

    def createUnitCellUniverse(self):
        """
        Constructs an empty universe (OrthrhombicPeriodicUniverse or
        ParallelepipedicPeriodicUniverse) representing the
        unit cell of the crystal. If the PDB file does not define
        a unit cell at all, an InfiniteUniverse is returned.
        
        :returns: a universe
        :rtype: :class:`~MMTK.Universe.Universe`
        """
        if self.from_fractional is None:
            return Universe.InfiniteUniverse()
        e1 = self.from_fractional(Vector(1., 0., 0.))
        e2 = self.from_fractional(Vector(0., 1., 0.))
        e3 = self.from_fractional(Vector(0., 0., 1.))
        if abs(e1.normal()*Vector(1., 0., 0.)-1.) < 1.e-15 \
               and abs(e2.normal()*Vector(0., 1., 0.)-1.) < 1.e-15 \
               and abs(e3.normal()*Vector(0., 0., 1.)-1.) < 1.e-15:
            return \
               Universe.OrthorhombicPeriodicUniverse((e1.length(),
                                                      e2.length(),
                                                      e3.length()))
        return Universe.ParallelepipedicPeriodicUniverse((e1, e2, e3))

    def createPeptideChains(self, model='all'):
        """
        :returns: a list of :class:`~MMTK.Proteins.PeptideChain` objects, one for each
                  peptide chain in the PDB file. The parameter model
                  has the same meaning as for the PeptideChain constructor.
        """
        return [chain.createPeptideChain(model)
                for chain in self.peptide_chains]

    def createNucleotideChains(self, model='all'):
        """
        :returns: a list of :class:`~MMTK.NucleicAcids.NucleotideChain` objects, one for each
                  nucleotide chain in the PDB file. The parameter model
                  has the same meaning as for the NucleotideChain constructor.
        """
        return [chain.createNucleotideChain(model)
                for chain in self.nucleotide_chains]

    def createMolecules(self, names = None, permit_undefined=True):
        """
        :param names: If a list of molecule names (as defined in the
                      chemical database) and/or PDB residue names,
                      only molecules mentioned in this list will be
                      constructed. If a dictionary, it is used to map
                      PDB residue names to molecule names. With the
                      default (None), only water molecules are
                      built.
        :type names: list
        :param permit_undefined: If False, an exception is raised
                                 when a PDB residue is encountered for
                                 which no molecule name is supplied
                                 in names. If True, an AtomCluster
                                 object is constructed for each unknown
                                 molecule.
        :returns: a collection of :class:`~MMTK.ChemicalObjects.Molecule` objects,
                  one for each molecule in the PDB file. Each PDB residue not 
                  describing an amino acid or nucleotide residue is considered a
                  molecule.
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        collection = Collections.Collection()
        mol_dicts = [molecule_names]
        if type(names) == type({}):
            mol_dicts.append(names)
            names = None
        for name in self.molecules.keys():
            full_name = None
            for dict in mol_dicts:
                full_name = dict.get(name, None)
            if names is None or name in names or full_name in names:
                if full_name is None and not permit_undefined:
                    raise ValueError("no definition for molecule " + name)
                for molecule in self.molecules[name]:
                    if full_name:
                        m = ChemicalObjects.Molecule(full_name)
                        setConfiguration(m, [molecule])
                    else:
                        pdbdict = {}
                        atoms = []
                        i = 0
                        for atom in molecule:
                            aname = atom.name
                            while aname[0] in string.digits:
                                aname = aname[1:] + aname[0]
                            try:
                                element = atom['element']
                                a = ChemicalObjects.Atom(element, name = aname)
                            except KeyError:
                                try:
                                    a = ChemicalObjects.Atom(aname[:2],
                                                             name = aname)
                                except IOError:
                                    a = ChemicalObjects.Atom(aname[:1],
                                                             name = aname)
                            a.setPosition(atom.position)
                            atoms.append(a)
                            pdbdict[atom.name] = Database.AtomReference(i)
                            i += 1
                        m = ChemicalObjects.AtomCluster(atoms, name = name)
                        if len(pdbdict) == len(molecule):
                            # pdbmap is correct only if the AtomCluster has
                            # unique atom names
                            m.pdbmap = [(name, pdbdict)]
                        setConfiguration(m, [molecule])
                    collection.addObject(m)
        return collection

    def createGroups(self, mapping):
        groups = []
        for name in self.molecules.keys():
            full_name = mapping.get(name, None)
            if full_name is not None:
                for molecule in self.molecules[name]:
                    g = ChemicalObjects.Group(full_name)
                    setConfiguration(g, [molecule], toplevel=0)
                    groups.append(g)
        return groups

    def createAll(self, molecule_names = None, permit_undefined=True):
        """
        :returns: a collection containing all objects contained in the
                  PDB file, i.e. the combination of the objects
                  returned by :func:`~MMTK.PDB.PDBConfiguration.createPeptideChains`,
                  :func:`~MMTK.PDB.PDBConfiguration.createNucleotideChains`,
                  and :func:`~MMTK.PDB.PDBConfiguration.createMolecules`.
                  The parameters have the same meaning as for
                  :func:`~MMTK.PDB.PDBConfiguration.createMolecules`.
        """
        collection = Collections.Collection()
        peptide_chains = self.createPeptideChains()
        if peptide_chains:
            import Proteins
            collection.addObject(Proteins.Protein(peptide_chains))
        nucleotide_chains = self.createNucleotideChains()
        collection.addObject(nucleotide_chains)
        molecules = self.createMolecules(molecule_names, permit_undefined)
        collection.addObject(molecules)
        return collection

    def asuToUnitCell(self, asu_contents, compact=True):
        """
        :param asu_contents: the molecules in the asymmetric unit, usually
                             obtained from :func:`~MMTK.PDB.PDBConfiguration.createAll()`.
        :param compact: if True, all molecules images are shifted such that
                        their centers of mass lie inside the unit cell.
        :type compact: bool
        :returns: a collection containing all molecules in the unit cell,
                  obtained by copying and moving the molecules from the
                  asymmetric unit according to the crystallographic
                  symmetry operations.
        :rtype: :class:`~MMTK.Collections.Collection`
        """
        unit_cell_contents = Collections.Collection()
        for symop in self.cs_transformations:
            transformation = symop.asLinearTransformation()
            rotation = transformation.tensor
            translation = transformation.vector
            image = copy.deepcopy(asu_contents)
            for atom in image.atomList():
                atom.setPosition(symop(atom.position()))
            if compact:
                cm = image.centerOfMass()
                cm_fr = self.to_fractional(cm)
                cm_fr = Vector(cm_fr[0] % 1., cm_fr[1] % 1., cm_fr[2] % 1.) \
                        - Vector(0.5, 0.5, 0.5)
                cm = self.from_fractional(cm_fr)
                image.translateTo(cm)
            unit_cell_contents.addObject(image)
        return unit_cell_contents

    def applyTo(self, object, atom_map=None):
        """
        Sets the configuration of object from the coordinates in the
        PDB file. The object must be compatible with the PDB file, i.e.
        contain the same subobjects and in the same order. This is usually
        only guaranteed if the object was created by the method

        :func:`~MMTK.PDB.PDBConfiguration.createAll` from a PDB file with the same layout.
        :param object: a chemical object or collection of chemical objects
        """
        setConfiguration(object, self.residues, atom_map=atom_map)

#
# An alternative name for compatibility in Database files.
#
PDBFile = PDBConfiguration

#
# Set atom coordinates from PDB configuration.
#
def setResidueConfiguration(object, pdb_residue, pdbmap, altmap,
                            atom_map = None):
    defined = 0
    for atom in pdb_residue:
        name = atom.name
        try: name = altmap[name]
        except KeyError: pass
        try:
            pdbname = pdbmap[1][name]
        except KeyError:
            pdbname = None
            if not object.isSubsetModel():
                raise IOError('Atom '+atom.name+' of PDB residue ' +
                               pdb_residue.name+' not found in residue ' +
                               pdbmap[0] + ' of object ' + object.fullName())
        if pdbname:
            object.setPosition(pdbname, atom.position)
            try:
                object.setIndex(pdbname, atom.number-1)
            except ValueError:
                pass
            if atom_map is not None:
                atom_map[object.getAtom(pdbname)] = atom
            defined += 1
    return defined

def setConfiguration(object, pdb_residues,
                     map = 'pdbmap', alt = 'pdb_alternative',
                     atom_map = None, toplevel = True):
    defined = 0
    if hasattr(object, 'is_protein'):
        i = 0
        for chain in object:
            l = len(chain)
            defined += setConfiguration(chain, pdb_residues[i:i+l],
                                        map, alt, atom_map, False)
            i = i + l
    elif hasattr(object, 'is_chain'):
        for i in range(len(object)):
            defined += setConfiguration(object[i], pdb_residues[i:i+1],
                                        map, alt, atom_map, False)
    elif hasattr(object, map):
        pdbmap = getattr(object, map)
        try: altmap = getattr(object, alt)
        except AttributeError: altmap = {}
        nres = len(pdb_residues)
        if len(pdbmap) != nres:
            raise IOError('PDB configuration does not match object ' +
                           object.fullName())
        for i in range(nres):
            defined += setResidueConfiguration(object, pdb_residues[i],
                                               pdbmap[i], altmap, atom_map)
    elif Collections.isCollection(object):
        nres = len(pdb_residues)
        if len(object) != nres:
            raise IOError('PDB configuration does not match object ' +
                           object.fullName())
        for i in range(nres):
            defined += setConfiguration(object[i], [pdb_residues[i]],
                                        map, alt, atom_map, False)
    else:
        try:
            name = object.fullName()
        except AttributeError:
            try:
                name = object.name
            except AttributeError:
                name = '???'
        raise IOError('PDB configuration does not match object ' + name)
              
    if toplevel and defined < object.numberOfAtoms():
        name = '[unnamed object]'
        try:
            name = object.fullName()
        except: pass
        if name: name = ' in ' + name
        Utility.warning(`object.numberOfAtoms()-defined` + ' atom(s)' + name +
                        ' were not assigned (new) positions.')
    return defined


#
# Create objects from a PDB configuration.
#
molecule_names = {'HOH': 'water', 'TIP': 'water', 'TIP3': 'water',
                  'WAT': 'water', 'HEM': 'heme'}

def defineMolecule(code, name):
    if molecule_names.has_key(code):
        raise ValueError("PDB code " + code + " already used")
    molecule_names[code] = name

#
# This object represents a PDB file for output.
#
class PDBOutputFile(object):

    """
    PDB file for output
    """

    def __init__(self, filename, subformat= None):
        """
        :param filename: the name of the PDB file that is created
        :type filename: str
        :param subformat: a variant of the PDB format; these subformats
                          are defined in module Scientific.IO.PDB. The
                          default is the standard PDB format.
        :type subformat: str
        """
        self.file = Scientific.IO.PDB.PDBFile(filename, 'w', subformat)
        self.warning = False
        self.atom_sequence = []
        self.model_number = None

    def nextModel(self):
        """
        Start a new model
        """
        if self.model_number is None:
            self.model_number = 1
        else:
            self.file.writeLine('ENDMDL', '')
            self.model_number = self.model_number + 1
        self.file.writeLine('MODEL', {'serial_number': self.model_number})

    def writeModel(self, object, configuration = None, slice_index = None, tag = None):
        """
        Write  an object to the file
        :param object: the object to be written
        :type object: :class:`~MMTK.Collections.GroupOfAtoms`
        :param configuration: the configuration from which the coordinates 
                              are taken (default: current configuration)
        :type configuration: :class:`~MMTK.ParticleProperties.Configuration`
        """
        if not ChemicalObjects.isChemicalObject(object):
            for o in object:
                self.writeModel(o, configuration, slice_index)
        else:
            toplevel = tag is None
            if toplevel:
                tag = Utility.uniqueAttribute()
            if hasattr(object, 'pdbmap'):
                for residue in object.pdbmap:
                    self.file.nextResidue(residue[0], )
                    sorted_atoms = residue[1].items()
                    sorted_atoms.sort(lambda x, y:
                                      cmp(x[1].number, y[1].number))
                    for atom_name, atom in sorted_atoms:
                        atom = object.getAtom(atom)
                        self.writeAtom(atom, atom_name, configuration, slice_index)
                        setattr(atom, tag, None)
            else:
                if hasattr(object, 'is_protein'):
                    for chain in object:                    
                        self.writeModel(chain, configuration, slice_index, tag)
                elif hasattr(object, 'is_chain'):
                        self.file.nextChain(None, object.name)
                        for residue in object:
                            self.writeModel(residue, configuration, slice_index, tag)
                        self.file.terminateChain()
                elif hasattr(object, 'molecules'):
                    for m in object.molecules:
                        self.writeModel(m, configuration, slice_index, tag)
                elif hasattr(object, 'groups'):
                    for g in object.groups:
                        self.writeModel(g, configuration, slice_index, tag)
            if toplevel:
                for a in object.atomList():
                    if not hasattr(a, tag):
                        self.writeModel(a, configuration, slice_index, tag)
                    delattr(a, tag)
    
    write = writeModel

    def writeAtom(self, atom, atom_name, configuration = None, slice_index = None):
        """
        Writes the position of the atom.
        """
        if slice_index is not None:
            bead_number = slice_index / (self.nbeads / atom.numberOfBeads())
            p = atom.beadPositions()[bead_number]
        else:
            p = atom.position(configuration)
        if Utility.isDefinedPosition(p):
            try: occ = atom.occupancy
            except AttributeError: occ = 0.
            try: temp = atom.temperature_factor
            except AttributeError: temp = 0.
            self.file.writeAtom(atom_name, p/Units.Ang,
                          occ, temp, atom.type.symbol)
            self.atom_sequence.append(atom)
        else:
            self.warning = True

    def close(self):
        """
        Closes the file. Must be called in order to prevent data loss.
        """
        if self.model_number is not None:
            self.file.writeLine('ENDMDL', '')
        self.file.close()
        if self.warning:
            Utility.warning('Some atoms are missing in the output file ' + \
                            'because their positions are undefined.')
            self.warning = False


#
# This object represents a PDB file for outputting path integrals to PDB.
#
class PDBOutputFileModel(PDBOutputFile):

    """
    PDB file for output
    """

    def write(self, object, configuration = None, tag = None):
        """
        Writes path integral PDBs with bead positions specified by model number. 
        """
        nbead_values = list(set(a.numberOfBeads() for a in object.atomList()))
        nbead_values.sort()
        self.nbeads = nbead_values[-1]
        for slice_index in range(self.nbeads):
            self.nextModel()
            self.writeModel(object, configuration, slice_index, tag)

#
# This object represents a PDB file for output.
#
class PDBOutputFileAltLoc(PDBOutputFile):

    """
    PDB file for output
    """

    def write(self, object, configuration = None, tag = None):
        """
        Writes path integral PDBs with bead positions specified by atom alternate locations.
        """
        nbead_values = list(set(a.numberOfBeads() for a in object.atomList()))
        nbead_values.sort()
        self.nbeads = nbead_values[-1]
        if self.nbeads > 36:
            raise IOError("At least one atom with greater than 36 beads, alternate location PDB write not supported, use pi_model")

        self.writeModel(object, configuration, tag = tag)

    def writeAtom(self, atom, atom_name, configuration = None, slice_index = None):
        """
        Writes the position of the atom, looping over beads
        """
        altloc_charmap = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
        for bead in range(atom.numberOfBeads()):
            p = atom.beadPositions()[bead]
            if Utility.isDefinedPosition(p):
                try: occ = atom.occupancy
                except AttributeError: occ = 0.
                try: temp = atom.temperature_factor
                except AttributeError: temp = 0.
                self.file.writeAtom(atom_name, p/Units.Ang,
                              occ, temp, atom.type.symbol, altloc_charmap[bead])
                self.atom_sequence.append(atom)
            else:
                self.warning = True

