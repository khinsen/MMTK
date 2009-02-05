# Deformation energy module
#
# Written by Konrad Hinsen
# last revision: 2009-2-5
#

"""
Deformation energies in proteins

This module implements deformational energies for use in the analysis
of motions and conformational changes in macromolecules. A description
of the techniques can be found in the following articles:

 1. K. Hinsen
    Analysis of domain motions by approximate normal mode calculations
    Proteins 33 (1998): 417-429

 2. K. Hinsen, A. Thomas, M.J. Field
    Analysis of domain motions in large proteins
    Proteins 34 (1999): 369-382
"""

__docformat__ = 'epytext'

try:
    from MMTK_forcefield import NonbondedList
    from MMTK_deformation import deformation, reduceDeformation, \
                                 reduceFiniteDeformation
except ImportError:
    pass
from MMTK import ParticleProperties
from Scientific import N

#
# Deformation energy evaluations
#
class DeformationEvaluationFunction(object):

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        self.universe = universe
        self.fc_length = fc_length
        self.cutoff = cutoff
        self.factor = factor

        nothing = N.zeros((0,2), N.Int)
        self.pairs = NonbondedList(nothing, nothing, nothing,
                                   universe._spec, cutoff)
        self.pairs.update(universe.configuration().array)
        self.normalize = 0
        try:
            self.version = self.forms.index(form)
        except ValueError:
            raise ValueError("unknown functional form")

    forms = ['exponential', 'calpha']

    def newConfiguration(self):
        self.pairs.update(self.universe.configuration().array)


class DeformationFunction(DeformationEvaluationFunction):

    """
    Infinite-displacement deformation function

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A DeformationFunction object must be called with a single parameter,
    which is a ParticleVector object containing the infinitesimal displacements
    of the atoms for which the deformation is to be evaluated.
    The return value is a ParticleScalar object containing the
    deformation value for each atom.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector):
        conf = self.universe.configuration()
        r = ParticleProperties.ParticleScalar(self.universe)
        l = deformation(conf.array, vector.array, self.pairs,
                        None, r.array, self.cutoff, self.fc_length,
                        self.factor, self.normalize, 0, self.version)
        return r


class NormalizedDeformationFunction(DeformationFunction):

    """
    Normalized infinite-displacement deformation function

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.
    The normalization is defined by equation 10 of reference 1 (see above).
    
    A NormalizedDeformationFunction object must be called with a single
    parameter, which is a ParticleVector object containing the infinitesimal
    displacements of the atoms for which the deformation is to be evaluated.
    The return value is a ParticleScalar object containing the
    deformation value for each atom.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationFunction.__init__(self, universe, fc_length,
                                     cutoff, factor, form)

    def __init__(self, *args, **kwargs):
        apply(DeformationFunction.__init__, (self, ) + args, kwargs)
        self.normalize = 1


class FiniteDeformationFunction(DeformationEvaluationFunction):

    """
    Finite-displacement deformation function

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A FiniteDeformationFunction object must be called with a single parameter,
    which is a Configuration or a ParticleVector object containing the
    alternate configuration of the universe for which the deformation is to be
    evaluated. The return value is a ParticleScalar object containing the
    deformation value for each atom.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector):
        conf = self.universe.configuration()
        vector = vector-conf
        r = ParticleProperties.ParticleScalar(self.universe)
        l = deformation(conf.array, vector.array, self.pairs, None, r.array,
                        self.cutoff, self.fc_length, self.factor, 0, 1,
                        self.version)
        return r


class DeformationEnergyFunction(DeformationEvaluationFunction):

    """
    Infinite-displacement deformation energy function

    The deformation energy is the sum of the deformation values over
    all atoms of a system.

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A DeformationEnergyFunction is called with one or two parameters.
    The first parameter is a ParticleVector object containing the
    infinitesimal displacements of the atoms for which the deformation
    energy is to be evaluated. The optional second argument can be
    set to a non-zero value to request the gradients of the energy
    in addition to the energy itself. In that case there are two
    return values (energy and the gradients in a ParticleVector
    object), otherwise only the energy is returned.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector, gradients = None):
        conf = self.universe.configuration()
        g = None
        if gradients is not None:
            if ParticleProperties.isParticleProperty(gradients):
                g = gradients
            elif isinstance(gradients, N.array_type):
                g = ParticleProperties.ParticleVector(self.universe, gradients)
            elif gradients:
                g = ParticleProperties.ParticleVector(self.universe)
        if g is None:
            g_array = None
        else:
            g_array = g.array
        l = deformation(conf.array, vector.array, self.pairs,
                        g_array, None, self.cutoff, self.fc_length,
                        self.factor, self.normalize, 0, self.version)
        if g is None:
            return l
        else:
            return l, g


class NormalizedDeformationEnergyFunction(DeformationEnergyFunction):

    """
    Normalized infinite-displacement deformation energy function

    The normalized deformation energy is the sum of the normalized
    deformation values over all atoms of a system.

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.
    The normalization is defined by equation 10 of reference 1.

    A NormalizedDeformationEnergyFunction is called with one or two parameters.
    The first parameter is a ParticleVector object containing the
    infinitesimal displacements of the atoms for which the deformation
    energy is to be evaluated. The optional second argument can be
    set to a non-zero value to request the gradients of the energy
    in addition to the energy itself. In that case there are two
    return values (energy and the gradients in a ParticleVector
    object), otherwise only the energy is returned.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEnergyFunction.__init__(self, universe, fc_length,
                                           cutoff, factor, form)
        self.normalize = 1


class FiniteDeformationEnergyFunction(DeformationEvaluationFunction):

    """
    Finite-displacement deformation energy function

    The deformation energy is the sum of the
    deformation values over all atoms of a system.

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A FiniteDeformationEnergyFunction is called with one or two parameters.
    The first parameter is a ParticleVector object containing the
    alternate configuration of the universe for which the deformation
    energy is to be evaluated. The optional second argument can be
    set to a non-zero value to request the gradients of the energy
    in addition to the energy itself. In that case there are two
    return values (energy and the gradients in a ParticleVector
    object), otherwise only the energy is returned.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector, gradients = None):
        conf = self.universe.configuration()
        g = None
        if gradients is not None:
            if ParticleProperties.isParticleProperty(gradients):
                g = gradients
            elif isinstance(gradients, N.array_type):
                g = ParticleProperties.ParticleVector(self.universe, gradients)
            elif gradients:
                g = ParticleProperties.ParticleVector(self.universe)
        if g is None:
            g_array = None
        else:
            g_array = g.array
        l = deformation(conf.array, vector.array, self.pairs, g_array, None,
                        self.cutoff, self.fc_length, self.factor,
                        0, 1, self.version)
        if g is None:
            return l
        else:
            return l, g

#
# Deformation energy minimization
#
class DeformationReducer(DeformationEvaluationFunction):

    """
    Iterative reduction of the deformation energy

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A DeformationReducer is called with two arguments. The first
    is a ParticleVector containing the initial infinitesimal displacements
    for all atoms. The second is an integer indicating the number of
    iterations. The result is a modification of the displacements
    by steepest-descent minimization of the deformation energy.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector, niter):
        conf = self.universe.configuration()
        reduceDeformation(conf.array, vector.array, self.pairs,
                          self.cutoff, self.fc_length, self.factor, niter,
                          self.version)


class FiniteDeformationReducer(DeformationEvaluationFunction):

    """
    Iterative reduction of the finite-displacement deformation energy

    The default values are appropriate for a C_alpha model of a protein
    with the global scaling described in the reference cited above.

    A FiniteDeformationReducer is called with two arguments. The first
    is a ParticleVector or Configuration containing the alternate
    configuration for which the deformation energy is evaluated.
    The second is the RMS distance that defines the termination
    condition. The return value a configuration that differs from
    the input configuration by approximately the specified RMS distance,
    and which is obtained by iterative steepest-descent minimization of
    the finite-displacement deformation energy.
    """

    def __init__(self, universe, fc_length = 0.7, cutoff = 1.2,
                 factor = 46402., form = 'exponential'):
        """
        @param universe: the universe for which the deformation function
                         is defined
        @type universe: L{MMTK.Universe.Universe}
        @param fc_length: the range parameter r_0 in the pair interaction term
        @type fc_length: C{float}
        @param cutoff: the cutoff used in the deformation calculation
        @type cutoff: C{float}
        @param factor: a global scaling factor
        @type factor: C{float}
        @param form: the functional form ('exponential' or 'calpha')
        @type form: C{str}
        """
        DeformationEvaluationFunction.__init__(self, universe, fc_length,
                                               cutoff, factor, form)

    def __call__(self, vector, rms_reduction):
        conf = self.universe.configuration()
        vector = vector-conf
        reduceFiniteDeformation(conf.array, vector.array, self.pairs,
                                self.cutoff, self.fc_length, self.factor,
                                rms_reduction, self.version)
        return conf+vector
