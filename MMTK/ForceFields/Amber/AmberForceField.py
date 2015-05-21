# This file provides the Amber force field, using Amber parameter files.
#
# Written by Konrad Hinsen
#

"""
Amber force field implementation

General comments about parameters for Lennard-Jones and electrostatic
interactions:

Pair interactions in periodic systems are calculated using the
minimum-image convention; the cutoff should therefore never be
larger than half the smallest edge length of the elementary
cell.

For Lennard-Jones interactions, all terms for pairs whose distance
exceeds the cutoff are set to zero, without any form of correction.
For electrostatic interactions, a charge-neutralizing surface charge
density is added around the cutoff sphere in order to reduce
cutoff effects (see Wolf et al., J. Chem. Phys. 17, 8254 (1999)).

For Ewald summation, there are some additional parameters that can
be specified by dictionary entries:

* "beta" specifies the Ewald screening parameter
* "real_cutoff" specifies the cutoff for the real-space sum.
   It should be significantly larger than 1/beta to ensure that
   the neglected terms are small.
* "reciprocal_cutoff" specifies the cutoff for the reciprocal-space sum.
   Note that, like the real-space cutoff, this is a distance; it describes
   the smallest wavelength of plane waves to take into account.
   Consequently, a smaller value means a more precise (and more expensive)
   calculation.

MMTK provides default values for these parameter which are calculated
as a function of the system size. However, these parameters are
exaggerated in most cases of practical interest and can lead to excessive
calculation times for large systems. It is preferable to determine
suitable values empirically for the specific system to be simulated.

The method "screened" uses the real-space part of the Ewald sum
with a charge-neutralizing surface charge density around the
cutoff sphere, and no reciprocal sum (see article cited above).
It requires the specification of the dictionary entries "cutoff"
and "beta".
"""

__docformat__ = 'restructuredtext'

from MMTK.ForceFields import MMForceField
from MMTK.ForceFields.Amber import AmberData
import os

#
# Read parameter files
#
Amber12SB = None
Amber14SB = None
Amber99 = None
Amber94 = None
Amber91 = None
OPLS = None

this_directory = os.path.split(__file__)[0]

def fullModFilePath(modfile):
    if not isinstance(modfile, basestring):
        return modfile
    if os.path.exists(modfile):
        return modfile
    return os.path.join(this_directory, os.path.basename(modfile))

def readAmberFiles(main_file, std_mod_file,
                   mod_files, atom_type_property, charge_property,
                   with_gaff):
    main_file = os.path.join(this_directory, main_file)
    if mod_files is None:
        mod_files = []
    else:
        mod_files = [(fullModFilePath(mf), 'MOD4') for mf in mod_files]
    if std_mod_file is not None:
        mod_files.insert(0, (os.path.join(this_directory, std_mod_file), 'MOD4'))
    params = AmberData.AmberParameters(main_file, mod_files)
    if with_gaff:
        gaff_params = AmberData.AmberParameters(os.path.join(this_directory,
                                                             "gaff.dat"), [])
        params.merge(gaff_params)
    params.lennard_jones_1_4 = 0.5
    params.electrostatic_1_4 = 1./1.2
    params.default_ljpar_set = params.ljpar_sets['MOD4']
    params.atom_type_property = atom_type_property
    params.charge_property = charge_property
    return params

def readAmber94(mod_files = None, with_gaff=False):
    global Amber94
    standard = mod_files is None
    if standard and Amber94 is not None:
        return Amber94
    params = readAmberFiles("parm94.dat", "frcmod.heme_ff94", mod_files,
                            'amber_atom_type', 'amber_charge', with_gaff)
    if standard:
        Amber94 = params
    return params

def readAmber99(mod_files = None, with_gaff=False):
    global Amber99
    standard = mod_files is None
    if standard and Amber99 is not None:
        return Amber99
    params = readAmberFiles("parm99.dat", None, mod_files,
                            'amber_atom_type', 'amber_charge', with_gaff)
    if standard:
        Amber99 = params
    return params

def readAmber12SB(mod_files = None, with_gaff=False):
    global Amber12SB
    standard = mod_files is None
    if standard and Amber12SB is not None:
        return Amber12SB
    params = readAmberFiles("parm10.dat", "frcmod.ff12SB", mod_files,
                            'amber12_atom_type', 'amber_charge', with_gaff)
    if standard:
        Amber12SB = params
    return params

def readAmber14SB(mod_files = None, with_gaff=False):
    global Amber14SB
    standard = mod_files is None
    if standard and Amber14SB is not None:
        return Amber14SB
    params = readAmberFiles("parm10.dat", "frcmod.ff14SB", mod_files,
                            'amber12_atom_type', 'amber_charge', with_gaff)
    if standard:
        Amber14SB = params
    return params

def readAmber91():
    global Amber91
    if Amber91 is None:
        Amber91 = AmberData.AmberParameters(os.path.join(this_directory,
                                                         "parm91.dat"))
        Amber91.lennard_jones_1_4 = 0.5
        Amber91.electrostatic_1_4 = 0.5
        Amber91.default_ljpar_set = Amber91.ljpar_sets['STDA']
        Amber91.atom_type_property = 'amber91_atom_type'
        Amber91.charge_property = 'amber91_charge'

def readOPLS(mod_files = None):
    global OPLS
    if mod_files is None:
        if OPLS is None:
            paramfile = os.path.join(this_directory, "opls_parm.dat")
            OPLS = AmberData.AmberParameters(paramfile)
            OPLS.lennard_jones_1_4 = 0.125
            OPLS.electrostatic_1_4 = 0.5
            OPLS.default_ljpar_set = OPLS.ljpar_sets['OPLS']
            OPLS.atom_type_property = 'opls_atom_type'
            OPLS.charge_property = 'opls_charge'
        return OPLS
    else:
        main_file = os.path.join(this_directory, "opls_parm.dat")
        if mod_files is None:
            mod_files = []
        else:
            mod_files = map(lambda mf: (fullModFilePath(mf), 'OPLS'), mod_files)
        params = AmberData.AmberParameters(main_file, mod_files)
        params.lennard_jones_1_4 = 0.125
        params.electrostatic_1_4 = 0.5
        params.default_ljpar_set = params.ljpar_sets['OPLS']
        params.atom_type_property = 'opls_atom_type'
        params.charge_property = 'opls_charge'
        return params

#
# The total force field
#
class Amber94ForceField(MMForceField.MMForceField):
    """
    Amber 94 force field
    """

    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor = 1., **kwargs):
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

        :keyword mod_files: a list of parameter modification files. The file
                            format is the one defined by AMBER. Each item
                            in the list can be either a file object
                            or a filename, filenames are looked up
                            first relative to the current directory and then
                            relative to the directory containing MMTK's
                            AMBER parameter files.

        """
        global Amber94
        mod_files = kwargs.get('mod_files', None)
        gaff = kwargs.get('gaff', False)
        parameters = readAmber94(mod_files, gaff)
        MMForceField.MMForceField.__init__(self, 'Amber94', parameters,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)

AmberForceField = Amber94ForceField


class Amber99ForceField(MMForceField.MMForceField):

    """
    Amber 99 force field
    """
    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor = 1., **kwargs):
 
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
    
        :keyword mod_files: a list of parameter modification files. The file
                            format is the one defined by AMBER. Each item
                            in the list can be either a file object
                            or a filename, filenames are looked up
                            first relative to the current directory and then
                            relative to the directory containing MMTK's
                            AMBER parameter files.
        """

        mod_files = kwargs.get('mod_files', None)
        gaff = kwargs.get('gaff', False)
        parameters = readAmber99(mod_files, gaff)
        MMForceField.MMForceField.__init__(self, 'Amber99', parameters,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)


class Amber12SBForceField(MMForceField.MMForceField):

    """
    Amber force field ff12SB
    """
    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor = 1., **kwargs):
 
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
    
        :keyword mod_files: a list of parameter modification files. The file
                            format is the one defined by AMBER. Each item
                            in the list can be either a file object
                            or a filename, filenames are looked up
                            first relative to the current directory and then
                            relative to the directory containing MMTK's
                            AMBER parameter files.
        """

        mod_files = kwargs.get('mod_files', None)
        gaff = kwargs.get('gaff', False)
        parameters = readAmber12SB(mod_files, gaff)
        MMForceField.MMForceField.__init__(self, 'Amber12SB', parameters,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)


class Amber14SBForceField(MMForceField.MMForceField):

    """
    Amber force field ff14SB
    """
    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor = 1., **kwargs):
 
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
    
        :keyword mod_files: a list of parameter modification files. The file
                            format is the one defined by AMBER. Each item
                            in the list can be either a file object
                            or a filename, filenames are looked up
                            first relative to the current directory and then
                            relative to the directory containing MMTK's
                            AMBER parameter files.
        """

        mod_files = kwargs.get('mod_files', None)
        gaff = kwargs.get('gaff', False)
        parameters = readAmber14SB(mod_files, gaff)
        MMForceField.MMForceField.__init__(self, 'Amber14SB', parameters,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)


class Amber91ForceField(MMForceField.MMForceField):

    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor=1., **kwargs):
        readAmber91()
        MMForceField.MMForceField.__init__(self, 'Amber91', Amber91,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)

class OPLSForceField(MMForceField.MMForceField):

    def __init__(self, lj_options = None, es_options = None,
                 bonded_scale_factor = 1., **kwargs):
        mod_files = kwargs.get('mod_files', None)
        parameters = readOPLS(mod_files)
        MMForceField.MMForceField.__init__(self, 'OPLS', parameters,
                                           lj_options, es_options,
                                           bonded_scale_factor)
        self.arguments = (lj_options, es_options, bonded_scale_factor)

#
# The following classes provides access to individual terms of the
# Amber force field. They are not used normally, but can be useful
# for testing.
#

#
# Bonded interactions
#
class AmberBondedForceField(MMForceField.MMBondedForceField):

    def __init__(self):
        readAmber94()
        MMForceField.MMBondedForceField.__init__(self, 'Amber_bonded', Amber94)
        self.arguments = ()

#
# Nonbonded interactions
#
class AmberLJForceField(MMForceField.MMLJForceField):

    def __init__(self, cutoff = None):
        readAmber94()
        MMForceField.MMLJForceField.__init__(self, 'Amber_LJ', Amber94, cutoff)
        self.arguments = (cutoff,)

class AmberESForceField(MMForceField.MMESForceField):

    def __init__(self, cutoff = None):
        readAmber94()
        MMForceField.MMESForceField.__init__(self, 'Amber_ES', Amber94, cutoff)
        self.arguments = (cutoff,)

class AmberEwaldESForceField(MMForceField.MMEwaldESForceField):

    def __init__(self, options = {}):
        readAmber94()
        MMForceField.MMEwaldForceESField.__init__(self, 'Amber_ES',
                                                  Amber94, options)
        self.arguments = (options,)

class AmberNonbondedForceField(MMForceField.MMNonbondedForceField):

    def __init__(self, lj_options = None, es_options = None):
        readAmber94()
        MMForceField.MMNonbondedForceField.__init__(self, 'Amber_NB', Amber94,
                                                    lj_options, es_options)
        self.arguments = (lj_options, es_options)
