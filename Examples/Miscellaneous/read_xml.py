# Read a universe from an XML file produced by the writeXML() method.
#
# The file used here is produced by running the example molecule_factory.py
#

from MMTK.XML import XMLMoleculeFactory

factory = XMLMoleculeFactory('water_universe.xml')
universe = factory.universe
