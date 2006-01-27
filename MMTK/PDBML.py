# Read PDBML files
#
# Written by Konrad Hinsen
# last revision: 2006-1-27
#

#
# This is experimental code, not yet a real MMTK module. Patience...
#

import elementtree.ElementTree as ET
#import cElementTree as ET
from Scientific.Geometry import Vector

class Structure:

    def __init__(self, filename, model = 1, alternate_code = 'A'):
        self.model_num = model
        self.alt_id = alternate_code
        self.prefix = None
        self.chem_comp = None
        self.entity_names = None
        self.entities = None
        self.chain_entities = None
        self.chains = None
        self.nonchains = None
        self.atoms = {}
        self.parseFile(filename)

    def parseFile(self, filename):
        for event, element in ET.iterparse(filename):
            tag = element.tag
            if self.prefix is None:
                self.prefix = tag[:tag.find('}')+1]
            if tag == self.prefix+"atom_site":
                atom_spec = self.parseAtom(element)
                if (atom_spec['alt_id'] is None or
                    atom_spec['alt_id'] == self.alt_id) \
                       and atom_spec['model_num'] == self.model_num:
                    atom = Atom()
                    self.atoms[atom_spec['atom_id']] = atom
                element.clear()
            elif tag == self.prefix+'chem_compCategory':
                self.chem_comp = element
            elif tag == self.prefix+'pdbx_entity_nameCategory':
                self.entity_names = element
            elif tag == self.prefix+'entityCategory':
                self.entities = element
            elif tag == self.prefix+'entity_poly_seqCategory':
                self.chain_entities = None
            elif tag == self.prefix+'pdbx_poly_seq_schemeCategory':
                self.chains = element
            elif tag == self.prefix+'pdbx_nonpoly_schemeCategory':
                self.nonchains = element
            elif tag == self.prefix+'atom_siteCategory':
                # This event happens when all atoms have been treated
                element.clear()

    def parseAtom(self, element):
        atom_spec = {
            'atom_id': element.get('id', None),
            'element': element.find(self.prefix+'type_symbol').text,
            'name': element.find(self.prefix+'label_atom_id').text,
            'alt_id': element.find(self.prefix+'label_alt_id').text,
            'model_num': int(element.find(self.prefix+'pdbx_PDB_model_num').text),
            'comp_id': element.find(self.prefix+'label_comp_id').text,
            'asym_id': element.find(self.prefix+'label_asym_id').text,
            'entity_id': element.find(self.prefix+'label_entity_id').text,
            'position': Vector(float(element.find(self.prefix+'Cartn_x').text),
                               float(element.find(self.prefix+'Cartn_y').text),
                               float(element.find(self.prefix+'Cartn_z').text)),
            'occupancy': float(element.find(self.prefix+'occupancy').text),
            'beta': float(element.find(self.prefix+'B_iso_or_equiv').text),
            }
        seq_id = element.find(self.prefix+'label_seq_id').text
        if seq_id is None:
            atom_spec['seq_id'] = None
        else:
            atom_spec['seq_id'] = int(seq_id)
        return atom_spec


class Atom:

    def __init__(self):
        pass


if __name__ == '__main__':
    
    s = Structure("/Users/hinsen/Desktop/1BTY.xml")
    #s = Structure("/Users/hinsen/Desktop/2BG9.xml")
