# Vibrational normal mode calculations.
#
# Written by Konrad Hinsen
#

"""
Vibrational normal modes
"""

__docformat__ = 'restructuredtext'

from MMTK import Features, Units, ParticleProperties
from Scientific import N

from MMTK.NormalModes import Core

#
# Class for a single mode
#
class VibrationalMode(Core.Mode):

    """
    Single vibrational normal mode

    Mode objects are created by indexing a 
    :class:`~MMTK.NormalModes.VibrationalModes.VibrationalModes` object.
    They contain the atomic displacements corresponding to a
    single mode. In addition, the frequency corresponding to the mode
    is stored in the attribute "frequency".

    Note: the normal mode vectors are *not* mass weighted, and therefore
    not orthogonal to each other.
    """

    def __init__(self, universe, n, frequency, mode):
        self.frequency = frequency
        Core.Mode.__init__(self, universe, n, mode)

    def __str__(self):
        return 'Mode ' + `self.number` + ' with frequency ' + `self.frequency`

    __repr__ = __str__

#
# Class for a full set of normal modes
#
class VibrationalModes(Core.NormalModes):

    """
    Vibrational modes describe the independent vibrational motions
    of a harmonic system. They are obtained by diagonalizing the mass-weighted
    force constant matrix.

    In order to obtain physically reasonable normal modes, the configuration
    of the universe must correspond to a local minimum of the potential
    energy.

    Individual modes (see :class:`~MMTK.NormalModes.VibrationalModes.VibrationalMode`)
    can be extracted by indexing with an integer. Looping over the modes
    is possible as well.

    """

    features = []

    def __init__(self, universe=None, temperature = 300.*Units.K,
                 subspace = None, delta = None, sparse = False):
        """
        :param universe: the system for which the normal modes are calculated;
                         it must have a force field which provides the second
                         derivatives of the potential energy
        :type universe: :class:`~MMTK.Universe.Universe`
        :param temperature: the temperature for which the amplitudes of the
                            atomic displacement vectors are calculated. A
                            value of None can be specified to have no scaling
                            at all. In that case the mass-weighted norm
                            of each normal mode is one.
        :type temperature: float
        :param subspace: the basis for the subspace in which the normal modes
                         are calculated (or, more precisely, a set of vectors
                         spanning the subspace; it does not have to be
                         orthogonal). This can either be a sequence of
                         :class:`~MMTK.ParticleProperties.ParticleVector` objects
                         or a tuple of two such sequences. In the second case,
                         the subspace is defined by the space spanned by the
                         second set of vectors projected on the complement of
                         the space spanned by the first set of vectors.
                         The first set thus defines directions that are
                         excluded from the subspace.
                         The default value of None indicates a standard
                         normal mode calculation in the 3N-dimensional
                         configuration space.
        :param delta: the rms step length for numerical differentiation.
                      The default value of None indicates analytical
                      differentiation.
                      Numerical differentiation is available only when a
                      subspace basis is used as well. Instead of calculating
                      the full force constant matrix and then multiplying
                      with the subspace basis, the subspace force constant
                      matrix is obtained by numerical differentiation of the
                      energy gradients along the basis vectors of the subspace.
                      If the basis is much smaller than the full configuration
                      space, this approach needs much less memory.
        :type delta: float
        :param sparse: a flag that indicates if a sparse representation of
                       the force constant matrix is to be used. This is of
                       interest when there are no long-range interactions and
                       a subspace of smaller size then 3N is specified. In that
                       case, the calculation will use much less memory with a
                       sparse representation.
        :type sparse: bool
        """

        if universe == None:
            return
        Features.checkFeatures(self, universe)
        Core.NormalModes.__init__(self, universe, subspace, delta, sparse,
                                  ['array', 'imaginary',
                                   'frequencies', 'omega_sq'])
        self.temperature = temperature

        self.weights = N.sqrt(self.universe.masses().array)
        self.weights = self.weights[:, N.NewAxis]

        self._forceConstantMatrix()
        ev = self._diagonalize()

        self.imaginary = N.less(ev, 0.)
        self.omega_sq = ev
        self.frequencies = N.sqrt(N.fabs(ev)) / (2.*N.pi)
        self.sort_index = N.argsort(self.frequencies)
        self.array.shape = (self.nmodes, self.natoms, 3)

        self.cleanup()

    def __setstate__(self, state):
        # This fixes unpickled objects from pre-2.5 versions.
        # Their modes were stored in un-weighted coordinates.
        if not state.has_key('weights'):
            weights = N.sqrt(state['universe'].masses().array)[:, N.NewAxis]
            state['weights'] = weights
            state['array'] *= weights[N.NewAxis, :, :]
            if state['temperature'] is not None:
                factor = N.sqrt(2.*state['temperature']*Units.k_B/Units.amu) / \
                         (2.*N.pi)
                freq = state['frequencies']
                for i in range(state['nmodes']):
                    index = state['sort_index'][i]
                    if index >= 6:
                        state['array'][index] *= freq[index]/factor
        self.__dict__.update(state)

    def __getitem__(self, item):
        index = self.sort_index[item]
        f = self.frequencies[index]
        if self.imaginary[index]:
            f = f*1.j
        # !!
        if self.temperature is None or item < 6:
            amplitude = 1.
        else:
            amplitude = N.sqrt(2.*self.temperature*Units.k_B/Units.amu) / \
                        (2.*N.pi*f)
        return VibrationalMode(self.universe, item, f,
                               amplitude*self.array[index]/self.weights)

    def rawMode(self, item):
        """
        :param item: the index of a normal mode
        :type item: int
        :returns: the unscaled mode vector
        :rtype: :class:`~MMTK.NormalModes.VibrationalModes.VibrationalMode`
        """
        index = self.sort_index[item]
        f = self.frequencies[index]
        if self.imaginary[index]:
            f = f*1.j
        return VibrationalMode(self.universe, item, f, self.array[index])

    def fluctuations(self, first_mode=6):
        f = ParticleProperties.ParticleScalar(self.universe)
        for i in range(first_mode, self.nmodes):
            mode = self.rawMode(i)
            f += (mode*mode)/mode.frequency**2
        f.array /= self.weights[:, 0]**2
        s = 1./(2.*N.pi)**2
        if self.temperature is not None:
            s *= Units.k_B*self.temperature
        f.array *= s
        return f

    def anisotropicFluctuations(self, first_mode=6):
        f = ParticleProperties.ParticleTensor(self.universe)
        for i in range(first_mode, self.nmodes):
            mode = self.rawMode(i)
            array = mode.array
            f.array += (array[:, :, N.NewAxis]*array[:, N.NewAxis, :]) \
                       / mode.frequency**2
        s = (2.*N.pi)**2
        if self.temperature is not None:
            s /= Units.k_B*self.temperature
        f.array /= s*self.universe.masses().array[:, N.NewAxis, N.NewAxis]
        return f

    def meanSquareDisplacement(self, subset=None, weights=None,
                               time_range = (0., None, None), tau=None,
                               first_mode = 6):
        """
        :param subset: the subset of the universe used in the calculation
                       (default: the whole universe)
        :type subset: :class:`~MMTK.Collections.GroupOfAtoms`
        :param weights: the weight to be given to each atom in the average
                        (default: atomic masses)
        :type weights: :class:`~MMTK.ParticleProperties.ParticleScalar`
        :param time_range: the time values at which the mean-square
                           displacement is evaluated, specified as a
                           range tuple (first, last, step).
                           The defaults are first=0, last=
                           20 times the longest vibration perdiod,
                           and step defined such that 300 points are
                           used in total.
        :type time_range: tuple
        :param tau: the relaxation time of an exponential damping factor
                    (default: no damping)
        :type tau: float
        :param first_mode: the first mode to be taken into account for
                           the fluctuation calculation. The default value
                           of 6 is right for molecules in vacuum.
        :type first_mode: int
        :returns: the averaged mean-square displacement of the
                  atoms in subset as a function of time
        :rtype: Scientific.Functions.Interpolation.InterpolatingFunction
        """
        if self.temperature is None:
            raise ValueError("no temperature available")
        if subset is None:
            subset = self.universe
        if weights is None:
            weights = self.universe.masses()
        weights = weights*subset.booleanMask()
        total = weights.sumOverParticles()
        weights = weights/total
        first, last, step = (time_range + (None, None))[:3]
        if last is None:
            last = 20./self[first_mode].frequency
        if step is None:
            step = (last-first)/300.
        time = N.arange(first, last, step)
        if tau is None: damping = 1.
        else: damping = N.exp(-(time/tau)**2)
        msd = N.zeros(time.shape, N.Float)
        for i in range(first_mode, self.nmodes):
            mode = self[i]
            omega = 2.*N.pi*mode.frequency
            d = (weights*(mode*mode)).sumOverParticles()
            N.add(msd, d*(1.-damping*N.cos(omega*time)), msd)
        return InterpolatingFunction((time,), msd)

    def EISF(self, q_range = (0., 15.), subset=None, weights = None,
             random_vectors = 15, first_mode = 6):
        """
        :param q_range: the range of angular wavenumber values
        :type q_range: tuple
        :param subset: the subset of the universe used in the calculation
                       (default: the whole universe)
        :type subset: :class:`~MMTK.Collections.GroupOfAtoms`
        :param weights: the weight to be given to each atom in the average
                        (default: incoherent scattering lengths)
        :type weights: :class:`~MMTK.ParticleProperties.ParticleScalar`
        :param random_vectors: the number of random direction vectors
                               used in the orientational average
        :type random_vectors: int
        :param first_mode: the first mode to be taken into account for
                           the fluctuation calculation. The default value
                           of 6 is right for molecules in vacuum.
        :type first_mode: int
        :returns: the Elastic Incoherent Structure Factor (EISF) as a
                  function of angular wavenumber
        :rtype: Scientific.Functions.Interpolation.InterpolatingFunction
        """
        if self.temperature is None:
            raise ValueError("no temperature available")
        if subset is None:
            subset = self.universe
        if weights is None:
            weights = self.universe.getParticleScalar('b_incoherent')
            weights = weights*weights
        weights = weights*subset.booleanMask()
        total = weights.sumOverParticles()
        weights = weights/total
    
        first, last, step = (q_range+(None,))[:3]
        if step is None:
            step = (last-first)/50.
        q = N.arange(first, last, step)
    
        f = MMTK.ParticleTensor(self.universe)
        for i in range(first_mode, self.nmodes):
            mode = self[i]
            f = f + mode.dyadicProduct(mode)
    
        eisf = N.zeros(q.shape, N.Float)
        for i in range(random_vectors):
            v = MMTK.Random.randomDirection()
            for a in subset.atomList():
                exp = N.exp(-v*(f[a]*v))
                N.add(eisf, weights[a]*exp**(q*q), eisf)
        return InterpolatingFunction((q,), eisf/random_vectors)

    def incoherentScatteringFunction(self, q, time_range = (0., None, None),
                                     subset=None, weights=None, tau=None,
                                     random_vectors=15, first_mode = 6):
        """
        :param q: the angular wavenumber
        :type q: float
        :param time_range: the time values at which the mean-square
                           displacement is evaluated, specified as a
                           range tuple (first, last, step).
                           The defaults are first=0, last=
                           20 times the longest vibration perdiod,
                           and step defined such that 300 points are
                           used in total.
        :type time_range: tuple
        :param subset: the subset of the universe used in the calculation
                       (default: the whole universe)
        :type subset: :class:`~MMTK.Collections.GroupOfAtoms`
        :param weights: the weight to be given to each atom in the average
                        (default: incoherent scattering lengths)
        :type weights: :class:`~MMTK.ParticleProperties.ParticleScalar`
        :param tau: the relaxation time of an exponential damping factor
                    (default: no damping)
        :type tau: float
        :param random_vectors: the number of random direction vectors
                               used in the orientational average
        :type random_vectors: int
        :param first_mode: the first mode to be taken into account for
                           the fluctuation calculation. The default value
                           of 6 is right for molecules in vacuum.
        :type first_mode: int
        :returns: the Incoherent Scattering Function as a function of time
        :rtype: Scientific.Functions.Interpolation.InterpolatingFunction
        """
        if self.temperature is None:
            raise ValueError("no temperature available")
        if subset is None:
            subset = self.universe
        if weights is None:
            weights = self.universe.getParticleScalar('b_incoherent')
            weights = weights*weights
        weights = weights*subset.booleanMask()
        total = weights.sumOverParticles()
        weights = weights/total
    
        first, last, step = (time_range + (None, None))[:3]
        if last is None:
            last = 20./self[first_mode].frequency
        if step is None:
            step = (last-first)/400.
        time = N.arange(first, last, step)
    
        if tau is None:
            damping = 1.
        else:
            damping = N.exp(-(time/tau)**2)
        finc = N.zeros(time.shape, N.Float)
        random_vectors = MMTK.Random.randomDirections(random_vectors)
        for v in random_vectors:
            for a in subset.atomList():
                msd = 0.
                for i in range(first_mode, self.nmodes):
                    mode = self[i]
                    d = (v*mode[a])**2
                    omega = 2.*N.pi*mode.frequency
                    msd = msd+d*(1.-damping*N.cos(omega*time))
                N.add(finc, weights[a]*N.exp(-q*q*msd), finc)
        return InterpolatingFunction((time,), finc/len(random_vectors))
