# Conversion of an MMTK trajectory into a sequence of PDB files.
#
# Note: In order to create the trajectory file "rotation.nc" which is
#       read by this example, you must run the example snapshot.py.
#

from MMTK import *
from MMTK.Trajectory import Trajectory

trajectory = Trajectory(None, "rotation.nc")
universe = trajectory.universe

for i, conf in enumerate(trajectory.configuration):
    universe.writeToFile("conf%d.pdb" % i, conf)
