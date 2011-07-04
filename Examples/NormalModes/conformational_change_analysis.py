# Analyze a conformational change in terms of the contributions of
# each normal mode.
#
# Note: download 1SU4 and 1T5T from the PDB before running this script

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.NormalModes import EnergeticModes
import pylab
import numpy

# Make a universe for the first configuration
universe = InfiniteUniverse(CalphaForceField())
universe.protein = Protein('1SU4.pdb', model='calpha')
conf_1SU4 = universe.copyConfiguration()

# Apply the second configuration and do a rigid-body superposition fit
PDBConfiguration('1T5T.pdb').applyTo(universe.protein)
tr, rms = universe.findTransformation(conf_1SU4)
universe.applyTransformation(tr)
conf_1T5T = universe.copyConfiguration()

# Set first configuration and calculate normal modes
universe.setConfiguration(conf_1SU4)
modes = EnergeticModes(universe, 300.*Units.K)

# Calculate the normalized difference vector
diff = (conf_1T5T - conf_1SU4).scaledToNorm(1.)

# Calculate the squared overlaps for all modes
mode_numbers = numpy.arange(6, len(modes))
overlaps = [modes.rawMode(i).dotProduct(diff)**2
            for i in mode_numbers]

# First plot: show that the cumulative sum converges to 1. This means that
# the squared overlaps can be considered as a percentage that each mode
# contributes to the conformational change.
pylab.figure(1)
pylab.plot(mode_numbers+1, numpy.add.accumulate(overlaps))
pylab.xlabel('mode number')
pylab.ylabel('cumulative squared overlap')

# Second plot: show the contribution of each of the first 100 modes
# to the conformational change.
pylab.figure(2)
pylab.plot(mode_numbers[:100]+1, overlaps[:100])
pylab.xlabel('mode number')
pylab.ylabel('squared overlap')
