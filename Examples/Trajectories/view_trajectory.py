# This program will show an animation of the trajectory file given as
# an argument. Optional parameters indicate the configurations to be
# included (first, last, step).
#

from MMTK.Visualization import viewTrajectory
import string, sys

first, last, step = 0, None, 1
argc = len(sys.argv)
if argc > 2:
    first = string.atoi(sys.argv[2])
    if argc > 3:
	last = string.atoi(sys.argv[3])
	if argc > 4:
	    step = string.atoi(sys.argv[4])

viewTrajectory(sys.argv[1], first, last, step)
