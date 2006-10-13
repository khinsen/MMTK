# An example of the use of the CCPN Data Model interface.
# Note that this interface is a first draft and thus likely to change!

from MMTK.CCPNDataModel import CCPNMolSystemAdapter
from memops.general.Io import loadXmlProjectFile

project = loadXmlProjectFile(file="/Users/hinsen/Programs/CCPN/projects/Test.xml")
mol_system = project.molSystems[0]

adapter = CCPNMolSystemAdapter(mol_system)
objects = adapter.makeAllMolecules()
for o in objects:
    print o.numberOfAtoms()
    print o.atomsWithDefinedPositions().centerOfMass()
