# This example shows how a universe can be built from a PDB file in such
# a way that all objects in the PDB file are represented as well as
# possible, using AtomCluster objects when nothing more specific can
# be constructed.
#
# This procedure is the only way to construct a universe that uses the
# same internal atom order as the PDB file, which is important for
# data exchange with other programs.

from MMTK import *
from MMTK.PDB import PDBConfiguration

configuration = PDBConfiguration('some_file.pdb')
universe = InfiniteUniverse()
universe.addObject(configuration.createAll(None, 1))
