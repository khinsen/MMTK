from MMTK import *
from MMTK.Subspace import Subspace
from Scientific.Statistics.Histogram import WeightedHistogram, Histogram
from Scientific import N

modes = load('~/proteins/lysozyme/lysozyme.umodes')
universe = modes.universe

protein = universe.protein
protein.model = 'all'
protein[0].model = 'all'

helices = []
for chain in protein:
    dihedrals = chain.phiPsi()[1:-1]
    dihedrals_with_index = zip(range(1, 1+len(dihedrals)), dihedrals)
    helix_indices = [index for index, (phi, psi) in dihedrals_with_index
                           if 4.5 < phi < 5.8 and 5. < psi < 6.]
    helix_indices = N.array(helix_indices)
    breaks = N.repeat(N.arange(len(helix_indices)-1),
                      (helix_indices[1:]-helix_indices[:-1]) > 1)
    breaks = N.concatenate(([0], breaks + 1))
    backbone = chain.backbone()
    for i in range(len(breaks)-1):
        residues = N.take(backbone, helix_indices[breaks[i]:breaks[i+1]])
        helices.append(Collection(list(residues)))
helices = [h for h in helices if len(h) > 4]

#Collection(helices).view()

residue_motion_vectors = []
helix_motion_vectors = []
for helix in helices:
    end_to_end = helix[0].centerOfMass()-helix[-1].centerOfMass()
    cms, inertia = helix.centerAndMomentOfInertia()
    moments, axes = inertia.diagonalization()
    axes = [a.asVector() for a in axes]
    helix_axis = axes[N.argmax([abs(end_to_end*v) for v in axes])]
    hmv = ParticleVector(universe)
    helix_motion_vectors.append(hmv)
    for residue in helix:
        mv = ParticleVector(universe)
        mv[residue.C_alpha] = helix_axis.cross(residue.C_alpha.position()-cms)
        hmv[residue.C_alpha] = helix_axis.cross(residue.C_alpha.position()-cms)
        residue_motion_vectors.append(mv)
    
residue_subspace = Subspace(universe, residue_motion_vectors)
helix_subspace = Subspace(universe, helix_motion_vectors)

frequencies = []
projections = []
for m in modes[6:]:
    frequencies.append(m.frequency)
    ms = m.scaledToNorm(1.)
    projections.append(residue_subspace.projectionOf(ms).norm()**2
                       - helix_subspace.projectionOf(ms).norm()**2)


histo = WeightedHistogram(frequencies, projections, 100)
plot(histo)
plot(Histogram(frequencies, 100))
