# Create molecules from scratch
#
# This example shows how a molecular system (SPCE water) can be set up using
# molecule factories instead of the molecular database.
#

from MMTK import *
from MMTK.MoleculeFactory import MoleculeFactory
from MMTK.ForceFields import SPCEForceField

# Create a new molecule factory.
factory = MoleculeFactory()

# Create the empty molecule. There is no distinction between groups
# and molecules at this level. Everything is a building block.
factory.createGroup('water')

# Add the atoms.
factory.addAtom('water', 'O', 'O')
factory.addAtom('water', 'H1', 'H')
factory.addAtom('water', 'H2', 'H')

# Add the bonds.
factory.addBond('water', 'O', 'H1')
factory.addBond('water', 'O', 'H2')

# Define the atom positions
factory.setPosition('water', 'O', Vector(0., 0., 0.00655616814675))
factory.setPosition('water', 'H1',
                    Vector(-0.0756950327264, 0., -0.0520320595151))
factory.setPosition('water', 'H2',
                    Vector(0.0756950327264, 0., -0.0520320595151))

# Define atom types and charges for the SPCE force field.
factory.setAttribute('water', 'O.spce_atom_type', 'O')
factory.setAttribute('water', 'H1.spce_atom_type', 'H')
factory.setAttribute('water', 'H2.spce_atom_type', 'H')
factory.setAttribute('water', 'O.spce_charge', -0.8476)
factory.setAttribute('water', 'H1.spce_charge', 0.4238)
factory.setAttribute('water', 'H2.spce_charge', 0.4238)

# Define the PDB atom name map
factory.setAttribute('water', 'pdbmap',
                     [('HOH', {'O': factory.getAtomReference('water', 'O'),
                               'H1': factory.getAtomReference('water', 'H1'),
                               'H2': factory.getAtomReference('water', 'H2')})])

# Create the universe and add water molecules.
universe = OrthorhombicPeriodicUniverse((5., 5., 5.), SPCEForceField())
universe.addObject(factory.retrieveMolecule('water'))
universe.addObject(factory.retrieveMolecule('water'))
universe[0].translateBy(Vector(1., 0., 0.))
universe[1].translateBy(Vector(0., 1., 0.))

# Print energy terms
print universe.energyTerms()

# Write the factory to an XML file
factory.writeXML(file('water_factory.xml', 'w'))

# Write the universe to an XML file
universe.writeXML(file('water_universe.xml', 'w'))

# Write the universe to a PDB file
universe.writeToFile('water_universe.pdb')
