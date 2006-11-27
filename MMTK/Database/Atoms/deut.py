from MMTK import *
from MMTK.ChemicalObjects import Group
from Scientific.N import sqrt

residues = ['alanine', 'arginine', 'asparagine', 'aspartic_acid',
            'cysteine', 'cystine_ss', 'glutamic_acid', 'glutamine',
            'glycine', 'histidine', 'histidine_deltah',
            'histidine_epsilonh', 'histidine_plus', 'isoleucine',
            'leucine', 'lysine', 'methionine', 'phenylalanine',
            'proline', 'serine', 'threonine', 'tryptophan',
            'tyrosine', 'valine']

for r in residues:
    binc = 0.
    binc_deut = 0.
    bcoh_sq = 0.
    bcoh_deut_sq = 0.
    mass = 0.
    for a in Group(r).atomList():
        binc = binc + a.b_incoherent/Units.fm
        binc_deut = binc_deut + a.b_incoherent_deut/Units.fm
        bcoh_sq = bcoh_sq + (a.b_coherent/Units.fm)**2
        bcoh_deut_sq = bcoh_deut_sq + (a.b_coherent_deut/Units.fm)**2
        mass = mass + a.mass()
    print r
    print "b_coherent = %f*fm" % sqrt(bcoh_sq)
    print "b_incoherent = %f*fm" % binc
    print "b_coherent_deut = %f*fm" % sqrt(bcoh_deut_sq)
    print "b_incoherent_deut = %f*fm" % binc_deut
