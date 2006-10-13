# Interface to the CCPN Data Model
#
# Warning: this is a first draft and thus likely to change!
# At the moment, it is a one-way interface that creates MMTK objects
# from CCPN MolSystems and MolStructures.
#
# Written by Konrad Hinsen
# last revision: 2006-10-13
#

"""Interface to the CCPN Data Model

The CCPN Data Model was created to facilitate exchange between programs
working with macromolecules and related experimental data. For more information,
see http://www.ccpn.ac.uk/datamodel/datamodel.html

Both the CCPN Data Model and MMTK's internal data model are
object-oriented representations of molecular systems and they share many
common features. However, there are important differences in terminology 
that can easily create confusion. Here is a rough comparison chart between
the central classes in both models. The equivalences are of course not perfect.


CCPN                    MMTK
====                    ====

Molecule                MoleculeType
MolResidue              BlueprintGroup
ChemComp/ChemCompVar    GroupType

MolSystem               set of Molecules without coordinates
MolSystem.Chain         Molecule without coordinates

MolStructure            set of Molecules with coordinates
MolStructure.CoordChain Molecule with coordinates

An important difference is that in the CCPN model, all molecules are chains of
residues. Non-polymeric molecules are represented as chains containing a single
residue.
"""

from MMTK.MoleculeFactory import MoleculeFactory
from MMTK import Units, Vector, Collection
from ccp.api.molecule.MolSystem import MolSystem, MolStructure
import sys

class CCPNMoleculeFactory(MoleculeFactory):

    """A MoleculeFactory built from the contents of a CCPN MolSystem.
    """

    def __init__(self, mol_system):
        MoleculeFactory.__init__(self)
        self.mol_system = mol_system
        self.group_name_mapping = {}
        self.makeAllGroups()

    def retrieveMoleculeForChain(self, mol_system_chain):
        return self.retrieveMolecule(self.group_name_mapping[mol_system_chain])
        
    def makeAllGroups(self):
        for chain in self.mol_system.chains:
            residues = chain.residues
            for residue in residues:
                chem_comp_var = residue.molResidue.chemCompVar
                self.makeGroupFromChemCompVar(chem_comp_var)
            mol_type = chain.molecule.molType
            if len(residues) > 1:
                self.makeGroupFromChain(chain)
                self.group_name_mapping[chain] = chain.molecule.name
            else:
                group_name = self.groupNameFromChemCompVar(chem_comp_var)
                self.group_name_mapping[chain] = group_name

    def groupNameFromChemCompVar(self, chem_comp_var):
        group_name = chem_comp_var.name
        descriptor = chem_comp_var.descriptor
        if descriptor != "neutral":
            group_name = group_name + '/' + descriptor
        linking = chem_comp_var.linking
        if linking != "none":
            group_name = group_name + '/' + linking
        return group_name

    def makeGroupFromChemCompVar(self, chem_comp_var):
        group_name = self.groupNameFromChemCompVar(chem_comp_var)
        if self.groups.has_key(group_name):
            return
        self.createGroup(group_name)
        for atom in chem_comp_var.findAllChemAtoms(className="ChemAtom"):
            self.addAtom(group_name, atom.name, atom.elementSymbol)
        for bond in chem_comp_var.chemBonds:
            atoms = bond.findAllChemAtoms(className="ChemAtom")
            if len(atoms) < 2:
                continue
            atom1 = atoms[0].name
            atom2 = atoms[1].name
            self.addBond(group_name, atom1, atom2)

    def makeGroupFromChain(self, chain):
        group_name = chain.molecule.name
        if self.groups.has_key(group_name):
            return
        self.createGroup(group_name)
        residue_ids = []
        for residue in chain.residues:
            chem_comp_var = residue.molResidue.chemCompVar
            residue_name = self.groupNameFromChemCompVar(chem_comp_var)
            residue_id = residue.ccpCode + str(residue.seqId)
            self.addSubgroup(group_name, residue_id, residue_name)
            residue_ids.append(residue_id)
        for link in chain.molecule.molResLinks:
            le1, le2 = link.molResLinkEnds
            la1 = le1.linkEnd.boundChemAtom
            la2 = le2.linkEnd.boundChemAtom
            s1 = le1.parent.serial
            s2 = le2.parent.serial
            self.addBond(group_name,
                         "%s.%s" % (residue_ids[s1-1], la1.name),
                         "%s.%s" % (residue_ids[s2-1], la2.name))


class CoordMapping(dict):

    """Mapoing from MolSystem.Atom objects to tuples of MolStructure.Coord
    objects taken from a corresponding MolStructure.
    """

    def __init__(self, mol_structure):
        for chain in mol_structure.coordChains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    self[atom.atom] = atom.coords

    def __getitem__(self, key):
        return self.get(key, ())


class CCPNMolSystemAdapter(object):
    
    """Interface between a CCPN MolSystem and corresponding MMTK
    objects. At the moment, this is one-way only (MMTK objects
    can be generated from CCPN MolSystems), but ultimately it will
    be two-way.
    """

    def __init__(self, mol_system):
        self.mol_system = mol_system
        self.mol_structures = mol_system.molStructures
        self.factory = CCPNMoleculeFactory(self.mol_system)
        self.coord_mappings = [CoordMapping(ms) for ms in self.mol_structures]
        self.ccpn_to_mmtk_atom_map = {}
        self.mmtk_to_ccpn_atom_map = {}
        self.ccpn_to_mmtk_molecule_map = {}
        self.mmtk_to_ccpn_molecule_map = {}

    def makeMolecule(self, ccpn_chain):
        molecule = self.factory.retrieveMoleculeForChain(ccpn_chain)
        self.ccpn_to_mmtk_molecule_map[ccpn_chain] = molecule
        self.mmtk_to_ccpn_molecule_map[molecule] = ccpn_chain
        residues = ccpn_chain.residues
        for residue in residues:
            if len(residues) == 1:
                my_residue = molecule
            else:
                residue_id = residue.ccpCode + str(residue.seqId)
                my_residue = getattr(molecule, residue_id)
            for atom in residue.atoms:
                my_atom = getattr(my_residue, atom.name)
                self.ccpn_to_mmtk_atom_map[atom] =  my_atom
                self.mmtk_to_ccpn_atom_map[my_atom] = atom
        if len(self.coord_mappings) == 1:
            self.assignPositions([molecule])
        return molecule

    def makeAllMolecules(self):
        molecules = Collection()
        for ccpn_chain in self.mol_system.chains:
            molecules.addObject(self.makeMolecule(ccpn_chain))
        return molecules

    def assignPositions(self, molecules, mol_structure_index=0):
        mapping = self.coord_mappings[mol_structure_index]
        for molecule in molecules:
            for atom in molecule.atomList():
                coords = mapping[self.mmtk_to_ccpn_atom_map[atom]]
                if coords:
                    position = Vector(coords[0].x, coords[0].y, coords[0].z)
                    atom.setPosition(position*Units.Ang)
        