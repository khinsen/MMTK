# Conversion of a DCD trajectory (CHARMM/XPlor) to MMTK's trajectory format.
#
# Note: In order to create the files "rotation.dcd" and "rotation.pdb" which
#       are read by this example, you must run the example dcd_export.py.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.PDB import PDBConfiguration
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.DCD import DCDReader

# Create the universe.
universe = InfiniteUniverse()

# Create all objects from the PDB file. The PDB file *must* match the
# the DCD file (same atom order), and you *must* create all objects
# defined in the PDB file, otherwise interpretation of the DCD file
# is not possible.
conf = PDBConfiguration('rotation.pdb')
chain = conf.peptide_chains[0].createPeptideChain(c_terminus=1)
universe.addObject(Protein([chain]))

# Create the trajectory object for output.
t = Trajectory(universe, "rotation_from_dcd.nc", "w", "Converted from DCD")

# Create a DCD reader that reads all configurations.
dcd_reader = DCDReader(universe, dcd_file = 'rotation.dcd',
                       actions = [TrajectoryOutput(t, "all", 0, None, 1)])
# Run the reader...
dcd_reader()
# ... and close the output trajectory.
t.close()
