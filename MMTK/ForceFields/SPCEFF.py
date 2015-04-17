# SPCE force field
#
# Written by Konrad Hinsen
#

"""
SPC/E force field for water simulations
"""

__docformat__ = 'restructuredtext'

from MMTK.ForceFields import MMForceField


class SPCEParameters(object):

    atom_type_property = 'spce_atom_type'
    charge_property = 'spce_charge'
    lennard_jones_1_4 = 1.
    electrostatic_1_4 = 1.

    def ljParameters(self, type):
        if type == 'O':
            # TIP3: 0.6363936, 0.315075240657
            return (0.650169580819, 0.31655578902, 0)
        elif type == 'H':
            return (0., 0., 0)
        else:
            raise ValueError('Unknown atom type ' + type)

    # Bond and angle parameters from:
    # O. Telemann, B. Jonsson, S. Engstrom
    # Mol. Phys. 60(1), 193-203 (1987)
    def bondParameters(self, at1, at2):
        if at1 == 'O' or at2 == 'O':
            return (0.1, 463700.)
        else:
            return (0.163298086184, 0.)
    
    def bondAngleParameters(self, at1, at2, at3):
        if at2 == 'O':
            return (1.91061193216, 383.)
        else:
            return (0.615490360716, 0.)
    
    def dihedralParameters(self, at1, at2, at3, at4):
        return [(1, 0., 0.)]
    
    def improperParameters(self, at1, at2, at3, at4):
        return [(1, 0., 0.)]


class SPCEForceField(MMForceField.MMForceField):

    """
    Force field for water simulations with the SPC/E model
    """

    def __init__(self, lj_options=None, es_options=None):
        """
        :param lj_options: parameters for Lennard-Jones
                           interactions; one of:

                            * a number, specifying the cutoff
                            * None, meaning the default method
                              (no cutoff; inclusion of all
                              pairs, using the minimum-image
                              conventions for periodic universes)
                            * a dictionary with an entry "method"
                              which specifies the calculation
                              method as either "direct" (all pair
                              terms) or "cutoff", with the cutoff
                              specified by the dictionary
                              entry "cutoff".

        :param es_options: parameters for electrostatic
                           interactions; one of:

                            * a number, specifying the cutoff
                            * None, meaning the default method
                              (all pairs without cutoff for
                              non-periodic system, Ewald summation
                              for periodic systems)
                            * a dictionary with an entry "method"
                              which specifies the calculation
                              method as either "direct" (all pair
                              terms), "cutoff" (with the cutoff
                              specified by the dictionary
                              entry "cutoff"), "ewald" (Ewald
                              summation, only for periodic
                              universes), or "screened".

        """
        MMForceField.MMForceField.__init__(self, 'SPCE', SPCEParameters(),
                                           lj_options, es_options)
        self.arguments = (lj_options, es_options)
