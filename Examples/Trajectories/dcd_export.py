# Conversion of an MMTK trajectory to DCD format (CHARMM/XPlor).
# A PDB file with the same atom order is created at the same time,
# because otherwise it would be impossible to interpret the data
# in the DCD file.
#
# Note: In order to create the trajectory file "rotation.nc" which is
#       read by this example, you must run the example snapshot.py.
#

from MMTK import *
from MMTK.Trajectory import Trajectory
from MMTK.DCD import writeDCDPDB


trajectory = Trajectory(None, "rotation.nc")
universe = trajectory.universe

writeDCDPDB(trajectory.configuration, 'rotation.dcd', 'rotation.pdb')
