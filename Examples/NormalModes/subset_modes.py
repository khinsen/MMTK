# A normal mode calculation in the subspace of a part of the total system.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.NormalModes import VibrationalModes
from MMTK.Subspace import Subspace
from MMTK.Minimization import ConjugateGradientMinimizer
from MMTK.Trajectory import StandardLogOutput

from Scientific import N

# Construct system
universe = InfiniteUniverse(Amber94ForceField())
universe.protein1 = Protein('bala1')
universe.protein2 = Protein('bala1')
universe.protein2.translateBy(Vector(1., 0., 0.))

# Define normal mode subspace
class SubsetSubspace(Subspace):

    def __init__(self, universe, objects):
        vectors = []
        for a in objects.atomList():
            for d in [Vector(1.,0.,0.), Vector(0.,1.,0.), Vector(0.,0.,1.)]:
                v = ParticleProperties.ParticleVector(universe)
                v[a] = d
                vectors.append(v)
        Subspace.__init__(self, universe, vectors)

subspace = SubsetSubspace(universe, universe.protein1)

# Minimize
minimizer = ConjugateGradientMinimizer(universe,
                                       actions=[StandardLogOutput(50)])
minimizer(convergence = 1.e-3, steps = 10000)

# Calculate normal modes
modes = VibrationalModes(universe, subspace=subspace)

# Print frequencies
for mode in modes:
    print mode

# Show animation of the first mode 
modes[0].view()
