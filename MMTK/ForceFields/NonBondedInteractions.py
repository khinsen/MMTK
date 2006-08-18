# This module implements classes that represent force fields
# for non-bonded interactions
#
# Written by Konrad Hinsen
# last revision: 2006-8-18
#

_undocumented = 1

from MMTK import Units, Utility
from ForceField import ForceField
from Scientific.Geometry import Vector
import Numeric

# Class definitions

#
# The base class NonBondedForceField provides the common
# functionality for all non-bonded interactions. The derived
# classes have to deal with determining functional forms
# and parameters and providing the evaluation code
#
class NonBondedForceField(ForceField):

    def __init__(self, name):
        ForceField.__init__(self, name)
        self.type = 'nonbonded'

    def ready(self, global_data):
        return 'bonded' in global_data.get('initialized')

    def excludedPairs(self, subset1, subset2, global_data):
        if 'excluded_pairs' not in global_data.get('initialized'):
            excluded_pairs = global_data.get('excluded_pairs')
            if subset1 is not None:
                for s1, s2 in [(subset1, subset2), (subset2, subset1)]:
                    set = {}
                    for a in s1.atomList():
                        set[a.index] = None
                    for a in s2.atomList():
                        try:
                            del set[a.index]
                        except KeyError: pass
                    excluded_pairs = excluded_pairs + \
                                     Utility.pairs(set.keys())
                set = {}
                for a in subset1.atomList():
                    set[a.index] = None
                for a in subset2.atomList():
                    set[a.index] = None
                atom_subset = set.keys()
                atom_subset.sort()
            else:
                atom_subset = None
            global_data.set('atom_subset', atom_subset)
            excluded_pairs = map(_normalizePair, excluded_pairs)
            excluded_pairs.sort(_cmpPair)
            _makeUnique(excluded_pairs)
            global_data.set('excluded_pairs', excluded_pairs)
            one_four_pairs = map(_normalizePair,
                                 global_data.get('1_4_pairs'))
            one_four_pairs.sort(_cmpPair)
            _makeUnique(one_four_pairs)
            _eliminateExcluded(one_four_pairs, excluded_pairs)
            global_data.set('1_4_pairs', one_four_pairs)
            global_data.add('initialized', 'excluded_pairs')
        return global_data.get('excluded_pairs'), \
               global_data.get('1_4_pairs'), \
               global_data.get('atom_subset')

    def nonbondedList(self, universe, subset1, subset2, global_data):
        try:
            from MMTK_forcefield import NonbondedList, NonbondedListTerm
        except ImportError:
            return None, None
        nbl = None
        update = None
        if 'nonbondedlist' in global_data.get('initialized'):
            nbl, update, cutoff = global_data.get('nonbondedlist')
        if nbl is None:
            excluded_pairs, one_four_pairs, atom_subset = \
                            self.excludedPairs(subset1, subset2, global_data)
            excluded_pairs = Numeric.array(excluded_pairs)
            one_four_pairs = Numeric.array(one_four_pairs)
            if atom_subset is not None:
                atom_subset = Numeric.array(atom_subset)
            else:
                atom_subset = Numeric.array([], Numeric.Int)
            nbl = NonbondedList(excluded_pairs, one_four_pairs, atom_subset,
                                universe._spec, self.cutoff)
            update = NonbondedListTerm(nbl)
            update.info = 0
            global_data.set('nonbondedlist', (nbl, update, self.cutoff))
            global_data.add('initialized', 'nonbondedlist')
        else:
            if cutoff is not None and \
                       (self.cutoff is None or self.cutoff > cutoff):
                nbl.setCutoff(self.cutoff)
        return nbl, update

    # the following methods must be overridden by derived classes
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        raise AttributeError

#
# Lennard-Jones force field
#
class LJForceField(NonBondedForceField):

    def __init__(self, name, cutoff, scale_factor=1.):
        NonBondedForceField.__init__(self, name)
        self.cutoff = cutoff
        self.scale_factor = scale_factor

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        n = universe.numberOfPoints()
        lj_type = -Numeric.ones((n,), Numeric.Int)
        atom_types = {}
        for o in universe:
            for a in o.atomList():
                atom_types[self._atomType(o, a, global_data)] = 1
        i = 0
        for t in atom_types.keys():
            atom_types[t] = i
            i = i + 1
        n_types = i
        for o in universe:
            for a in o.atomList():
                lj_type[a.index] = atom_types[self._atomType(o,a,global_data)]
        eps_sigma = Numeric.zeros((n_types, n_types, 2), Numeric.Float)
        for t, i in atom_types.items():
            eps, sigma, mix = self._ljParameters(t, global_data)
            eps = eps*self.scale_factor
            eps_sigma[i,i] = Numeric.array([eps, sigma])
        eps = eps_sigma[:,:,0]
        sigma = eps_sigma[:,:,1]
        for i in range(n_types):
            for j in range(i+1, n_types):
                eps[i,j] = Numeric.sqrt(eps[i,i]*eps[j,j])
                eps[j,i] = eps[i,j]
                if mix == 0:
                    sigma[i,j] = 0.5*(sigma[i,i]+sigma[j,j])
                elif mix == 1:
                    sigma[i,j] = Numeric.sqrt(sigma[i,i]*sigma[j,j])
                else:
                    raise ValueError("undefined Lennard-Jones mixing rule")
                sigma[j,i] = sigma[i,j]
        global_data.lj_type = lj_type
        global_data.lj_14_factor = self.lj_14_factor
        global_data.eps_sigma = eps_sigma
        nblist, update = \
                   self.nonbondedList(universe, subset1, subset2, global_data)
        if self.cutoff is None:
            cutoff = 0.
        else:
            cutoff = self.cutoff
        from MMTK_forcefield import LennardJonesTerm
        ev = LennardJonesTerm(universe._spec, nblist, eps_sigma,
                              lj_type, cutoff, self.lj_14_factor)
        update.addTerm(ev, 0);
        if update.info:
            return [ev]
        else:
            update.info = 1
            return [update, ev]

    # the following methods must be overridden by derived classes
    def _atomType(self, o, a, global_data):
        raise AttributeError
    def _ljParameters(self, t, global_data):
        raise AttributeError

def _normalizePair(pair):
    i, j = pair
    if i > j:
        return j, i
    else:
        return i, j

def _cmpPair(p1, p2):
    cmp0 = cmp(p1[0], p2[0])
    if not cmp0:
        return cmp(p1[1], p2[1])
    else:
        return cmp0

def _makeUnique(pair_list):
    i = 0
    while i < len(pair_list):
        item = pair_list[i]
        j = i + 1
        while j < len(pair_list) and pair_list[j] == item:
            del pair_list[j]
        i = i + 1

def _eliminateExcluded(one_four_list, excluded_list):
    i = 0
    j = 0
    n = len(excluded_list)
    while i < len(one_four_list):
        item = one_four_list[i]
        while j < n and item > excluded_list[j]:
            j = j + 1
        if j == n:
            break
        if item == excluded_list[j]:
            del one_four_list[i]
        else:
            i = i + 1

#
# Electrostatic force field
#
class ElectrostaticForceField(NonBondedForceField):

    def __init__(self, name, cutoff, scale_factor=1.):
        NonBondedForceField.__init__(self, name)
        self.cutoff = cutoff
        self.scale_factor = scale_factor

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        n = universe.numberOfPoints()
        charge = Numeric.zeros((n,), Numeric.Float)
        for o in universe:
            for a in o.atomList():
                charge[a.index] = self._charge(o, a, global_data)
        charge = charge*Numeric.sqrt(self.scale_factor)
        nblist, update = \
                     self.nonbondedList(universe, subset1, subset2,
                                        global_data)
        if self.cutoff is None:
            cutoff = 0.
        else:
            cutoff = self.cutoff
        from MMTK_forcefield import ElectrostaticTerm
        ev = ElectrostaticTerm(universe._spec, nblist, charge,
                               cutoff, self.es_14_factor)
        update.addTerm(ev, 1);
        if update.info:
            return [ev]
        else:
            update.info = 1
            return [update, ev]

    # the following method must be overridden by derived classes
    def _charge(self, o, a, global_data):
        raise AttributeError

#
# Ewald evaluator for electrostatic interactions
#
class ESEwaldForceField(NonBondedForceField):

    def __init__(self, name, options = {}):
        NonBondedForceField.__init__(self, name)
        self.cutoff = options.get('real_cutoff', None)
        self.scale_factor = options.get('scale_factor', 1.)
        self.options = options
        for key in options.keys():
            if key not in self.known_options:
                raise ValueError(key + " is not a recognized option")

    known_options = ['beta', 'real_cutoff', 'cutoff', 'reciprocal_cutoff',
                     'ewald_precision', 'no_reciprocal_sum', 'method',
                     'scale_factor']

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        rsum = not self.options.get('no_reciprocal_sum', 0)
        if not universe.is_periodic and rsum:
            raise ValueError("Ewald methods accepts only periodic universes")
        n = universe.numberOfPoints()
        charge = Numeric.zeros((n,), Numeric.Float)
        for o in universe:
            for a in o.atomList():
                charge[a.index] = self._charge(o, a, global_data)
        charge = charge*Numeric.sqrt(self.scale_factor)
        precision = self.options.get('ewald_precision', 1.e-6)
        p = Numeric.sqrt(-Numeric.log(precision))
        if rsum:
            beta_opt = Numeric.sqrt(Numeric.pi) * \
                (5.*universe.numberOfAtoms()/universe.cellVolume()**2)**(1./6.)
            max_cutoff = universe.largestDistance()
            beta_opt = max(p/max_cutoff, beta_opt)
        else:
            beta_opt = 0.01
        options = {}
        options['beta'] = beta_opt
        options['real_cutoff'] = p/beta_opt
        options['reciprocal_cutoff'] = Numeric.pi/(beta_opt*p)
        options['no_reciprocal_sum'] = 0
        for key, value in self.options.items():
            options[key] = value
        lx = universe.boxToRealCoordinates(Vector(1., 0., 0.)).length()
        ly = universe.boxToRealCoordinates(Vector(0., 1., 0.)).length()
        lz = universe.boxToRealCoordinates(Vector(0., 0., 1.)).length()
        kmax = Numeric.array([lx,ly,lz])/options['reciprocal_cutoff']
        kmax = Numeric.ceil(kmax).astype(Numeric.Int)
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        if atom_subset is not None:
            raise ValueError("Ewald summation not available for subsets")
        nblist, update = \
                self.nonbondedList(universe, subset1, subset2,global_data)
        if options['no_reciprocal_sum']:
            kcutoff = 0.
            shape = Numeric.zeros((3, 3), Numeric.Float)
        else:
            kcutoff = (2.*Numeric.pi/options['reciprocal_cutoff'])**2
            shape = universe.basisVectors()
            if shape is None:
                raise ValueError("Ewald evaluator needs periodic universe")
            shape = Numeric.array(shape, Numeric.Float)
        from MMTK_forcefield import EsEwaldTerm
        ev = EsEwaldTerm(universe._spec, shape, nblist, charge,
                         options['real_cutoff'], kcutoff, kmax,
                         self.es_14_factor, options['beta'])
        update.addTerm(ev, 2);
        if update.info:
            return [ev]
        else:
            update.info = 1
            return [update, ev]

    # the following methods must be overridden by derived classes
    def _charge(self, o, a, global_data):
        raise AttributeError

#
# Multipole evaluator for electrostatic interactions
#
class ESMPForceField(NonBondedForceField):

    def __init__(self, name, options = {}):
        NonBondedForceField.__init__(self, name)
        self.cutoff = None
        self.scale_factor = options.get('scale_factor', 1.)
        self.options = options
        for key in options.keys():
            if key not in self.known_options:
                raise ValueError(key + " is not a recognized option")

    known_options = ['spatial_decomposition_levels',
                     'multipole_expansion_terms',
                     'use_fft',
                     'fft_blocking_factor',
                     'macroscopic_expansion_terms',
                     'multipole_acceptance',
                     'method',
                     'scale_factor']

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if universe.is_periodic:
            try:
                shape = universe.basisVectors()
            except AttributeError:
                raise ValueError("Multipole method implemented only " +
                                  "for orthorhombic universes.")
        else:
            shape = None
        n = universe.numberOfPoints()
        charge = Numeric.zeros((n,), Numeric.Float)
        atom_types = {}
        for o in universe:
            for a in o.atomList():
                charge[a.index] = self._charge(o, a, global_data)
        charge = Numeric.zeros((n,), Numeric.Float)
        options = {}
        if n < 10000:
            options['spatial_decomposition_levels'] = 4
        elif n < 100000:
            options['spatial_decomposition_levels'] = 5
        else:
            options['spatial_decomposition_levels'] = 6
        options['multipole_expansion_terms'] = 8
        options['use_fft'] = 0
        options['fft_blocking_factor'] = 4
        options['macroscopic_expansion_terms'] = 6
        options['multipole_acceptance'] = 0.5
        for key, value in self.options.items():
            options[key] = value
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        excluded_pairs = Numeric.array(excluded_pairs)
        one_four_pairs = Numeric.array(one_four_pairs)
        if atom_subset is None:
            atom_subset = Numeric.array([], Numeric.Int)
        else:
            atom_subset = Numeric.array(atom_subset)
        nbinfo = [excluded_pairs, one_four_pairs, atom_subset]
        if shape is None:
            shape = Numeric.zeros((), Numeric.Float)
        else:
            shape = Numeric.array(shape, Numeric.Float)
        from MMTK_forcefield import EsMPTerm
        ev = EsMPTerm(universe._spec, shape, nbinfo,
                      charge, self.es_14_factor,
                      options['spatial_decomposition_levels'],
                      options['multipole_expansion_terms'],
                      options['use_fft'],
                      options['fft_blocking_factor'],
                      options['macroscopic_expansion_terms'],
                      options['multipole_acceptance'])
        return [ev]

    # the following methods must be overridden by derived classes
    def _charge(self, o, a, global_data):
        raise AttributeError
