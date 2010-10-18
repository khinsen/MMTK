# This module implements classes that represent force fields
# for non-bonded interactions
#
# Written by Konrad Hinsen
#

from MMTK import Units, Utility
from MMTK.ForceFields.ForceField import ForceField
from Scientific.Geometry import Vector
from Scientific import N

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

    def supportsPathIntegrals(self):
        return True

    def excludedPairs(self, subset1, subset2, global_data):
        if 'excluded_pairs' not in global_data.get('initialized'):
            excluded_pairs = set(global_data.get('excluded_pairs'))
            if subset1 is not None:
                set1 = set(a.index for a in subset1.atomList())
                set2 = set(a.index for a in subset2.atomList())
                added_excluded_pairs = set(Utility.orderedPairs(list(set1-set2))) \
                                       | set(Utility.orderedPairs(list(set2-set1)))
                added_excluded_pairs = sum(([(i+oi, j+oj)
                                             for oi, oj in self.beadOffsetsAndFactor((i, j),
                                                                                     global_data)[1]]
                                            for i, j in added_excluded_pairs), [])
                excluded_pairs |= set(added_excluded_pairs)
                bead_data = global_data.get('bead_data')
                atom_subset = sum((range(i, i+bead_data[i, 1]) for i in list(set1 | set2)), [])
                atom_subset.sort()
            else:
                atom_subset = None
            global_data.set('atom_subset', atom_subset)
            global_data.set('excluded_pairs', list(excluded_pairs))
            one_four_pairs = set(global_data.get('1_4_pairs')) \
                               - excluded_pairs
            global_data.set('1_4_pairs', list(one_four_pairs))
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
            excluded_pairs = N.array(excluded_pairs)
            one_four_pairs = N.array(one_four_pairs)
            if atom_subset is not None:
                atom_subset = N.array(atom_subset)
            else:
                atom_subset = N.array([], N.Int)
            nbl = NonbondedList(excluded_pairs, one_four_pairs, atom_subset,
                                global_data.get('nbeads'), global_data.get('bead_data'),
                                universe._spec, self.cutoff)
            update = NonbondedListTerm(nbl)
            update.info = False
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

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        n = universe.numberOfPoints()
        lj_type = -N.ones((n,), N.Int)
        atom_types = {}
        atom_type_names = []
        for o in universe:
            for a in o.atomIterator():
                t = self._atomType(o, a, global_data)
                if not atom_types.has_key(t):
                    atom_types[t] = len(atom_types)
                    atom_type_names.append(t)
                for b in a.beads():
                    lj_type[b.index] = atom_types[t]
        n_types = len(atom_types)
        eps_sigma = N.zeros((n_types, n_types, 2), N.Float)
        for t, i in atom_types.items():
            eps, sigma, mix = self._ljParameters(t, global_data)
            eps *= self.scale_factor
            eps_sigma[i, i] = N.array([eps, sigma])
        eps = eps_sigma[:, :, 0]
        sigma = eps_sigma[:, :, 1]
        for i in range(n_types):
            for j in range(i+1, n_types):
                eps[i,j] = N.sqrt(eps[i,i]*eps[j,j])
                eps[j,i] = eps[i,j]
                if mix == 0:
                    sigma[i,j] = 0.5*(sigma[i,i]+sigma[j,j])
                elif mix == 1:
                    sigma[i,j] = N.sqrt(sigma[i,i]*sigma[j,j])
                else:
                    raise ValueError("unknown Lennard-Jones mixing rule")
                sigma[j,i] = sigma[i,j]
        global_data.lj_type = lj_type
        global_data.lj_14_factor = self.lj_14_factor
        global_data.eps_sigma = eps_sigma
        if self.cutoff is None:
            cutoff = 0.
        else:
            cutoff = self.cutoff
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        return {'lennard_jones': {'type': lj_type,
                                  'type_names': atom_type_names,
                                  'epsilon_sigma': eps_sigma,
                                  'one_four_factor': self.lj_14_factor,
                                  'cutoff': cutoff},
                'nonbonded': {'excluded_pairs': excluded_pairs,
                              'one_four_pairs': one_four_pairs,
                              'atom_subset': atom_subset}
               }

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['lennard_jones']
        nblist, update = \
                   self.nonbondedList(universe, subset1, subset2, global_data)
        from MMTK_forcefield import LennardJonesTerm
        ev = LennardJonesTerm(universe._spec, nblist, params['epsilon_sigma'],
                              params['type'], params['cutoff'],
                              params['one_four_factor'])
        update.addTerm(ev, 0);
        if update.info:
            return [ev]
        else:
            update.info = True
            return [update, ev]

    # the following methods must be overridden by derived classes
    def _atomType(self, o, a, global_data):
        raise AttributeError
    def _ljParameters(self, t, global_data):
        raise AttributeError

#
# Electrostatic force field
#
class ElectrostaticForceField(NonBondedForceField):

    def __init__(self, name, cutoff, scale_factor=1.):
        NonBondedForceField.__init__(self, name)
        self.cutoff = cutoff
        self.scale_factor = scale_factor

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        n = universe.numberOfPoints()
        charge = N.zeros((n,), N.Float)
        for o in universe:
            for a in o.atomIterator():
                c = self._charge(o, a, global_data)
                for b in a.beads():
                    charge[b.index] = c
        charge = charge*N.sqrt(self.scale_factor)
        if self.cutoff is None:
            cutoff = 0.
        else:
            cutoff = self.cutoff
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        return {'electrostatic': {'algorithm':'direct',
                                  'charge': charge,
                                  'cutoff': cutoff,
                                  'one_four_factor': self.es_14_factor},
                'nonbonded': {'excluded_pairs': excluded_pairs,
                              'one_four_pairs': one_four_pairs,
                              'atom_subset': atom_subset}
               }

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['electrostatic']
        assert params['algorithm'] == 'direct'
        nblist, update = \
                    self.nonbondedList(universe, subset1, subset2, global_data)
        from MMTK_forcefield import ElectrostaticTerm
        ev = ElectrostaticTerm(universe._spec, nblist, params['charge'],
                               params['cutoff'], params['one_four_factor'])
        update.addTerm(ev, 1)
        if update.info:
            return [ev]
        else:
            update.info = True
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

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        rsum = not self.options.get('no_reciprocal_sum', False)
        if not universe.is_periodic and rsum:
            raise ValueError("Ewald method accepts only periodic universes")
        n = universe.numberOfPoints()
        charge = N.zeros((n,), N.Float)
        for o in universe:
            for a in o.atomIterator():
                c = self._charge(o, a, global_data)
                for b in a.beads():
                    charge[b.index] = c
        charge = charge*N.sqrt(self.scale_factor)
        precision = self.options.get('ewald_precision', 1.e-6)
        p = N.sqrt(-N.log(precision))
        if rsum:
            beta_opt = N.sqrt(N.pi) * \
                (5.*universe.numberOfAtoms()/universe.cellVolume()**2)**(1./6.)
            max_cutoff = universe.largestDistance()
            beta_opt = max(p/max_cutoff, beta_opt)
        else:
            beta_opt = 0.01
        options = {}
        options['beta'] = beta_opt
        options['real_cutoff'] = p/beta_opt
        options['reciprocal_cutoff'] = N.pi/(beta_opt*p)
        options['no_reciprocal_sum'] = False
        for key, value in self.options.items():
            options[key] = value
        kx, ky, kz = [1./(rbv.length()*options['reciprocal_cutoff'])
                      for rbv in universe.reciprocalBasisVectors()]
        kmax = N.ceil([kx, ky, kz]).astype(N.Int)
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        if atom_subset is not None:
            raise ValueError("Ewald summation not available for subsets")
        if options['no_reciprocal_sum']:
            kcutoff_sq = 0.
        else:
            kcutoff_sq = (2.*N.pi/options['reciprocal_cutoff'])**2
        return {'electrostatic': {'algorithm': 'ewald',
                                  'charge': charge,
                                  'real_cutoff': options['real_cutoff'],
                                  'k_cutoff_sq': kcutoff_sq,
                                  'beta': options['beta'],
                                  'k_max': kmax,
                                  'one_four_factor': self.es_14_factor},
                'nonbonded': {'excluded_pairs': excluded_pairs,
                              'one_four_pairs': one_four_pairs,
                              'atom_subset': atom_subset}
               }

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['electrostatic']
        assert params['algorithm'] == 'ewald'
        nblist, update = \
                self.nonbondedList(universe, subset1, subset2, global_data)
        if params['k_cutoff_sq'] == 0.:
            shape = N.zeros((3, 3), N.Float)
        else:
            shape = universe.basisVectors()
            if shape is None:
                raise ValueError("Ewald evaluator needs periodic universe")
            shape = N.array(shape, N.Float)
        from MMTK_forcefield import EsEwaldTerm
        ev = EsEwaldTerm(universe._spec, shape, nblist, params['charge'],
                         params['real_cutoff'], params['k_cutoff_sq'],
                         params['k_max'], params['one_four_factor'],
                         params['beta'],
                         global_data.get('bead_data'),
                         global_data.get('nbeads'))
        update.addTerm(ev, 2)
        if update.info:
            return [ev]
        else:
            update.info = True
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

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        n = universe.numberOfPoints()
        charge = N.zeros((n,), N.Float)
        atom_types = {}
        for o in universe:
            for a in o.atomIterator():
                c = self._charge(o, a, global_data)
                for b in a.beads():
                    charge[b.index] = c
        charge = N.zeros((n,), N.Float)
        params = {}
        if n < 10000:
            params['spatial_decomposition_levels'] = 4
        elif n < 100000:
            params['spatial_decomposition_levels'] = 5
        else:
            params['spatial_decomposition_levels'] = 6
        params['multipole_expansion_terms'] = 8
        params['use_fft'] = 0
        params['fft_blocking_factor'] = 4
        params['macroscopic_expansion_terms'] = 6
        params['multipole_acceptance'] = 0.5
        for key, value in self.options.items():
            params[key] = value
        params['algorithm'] = 'dpmta'
        params['charge'] = charge
        params['one_four_factor'] = self.es_14_factor
        excluded_pairs, one_four_pairs, atom_subset = \
                             self.excludedPairs(subset1, subset2, global_data)
        return {'electrostatic': params,
                'nonbonded': {'excluded_pairs': excluded_pairs,
                              'one_four_pairs': one_four_pairs,
                              'atom_subset': atom_subset}
               }

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)
        pe = params['electrostatic']
        assert pe['algorithm'] == 'dpmta'
        pn = params['nonbonded']
        excluded_pairs = N.array(pn['excluded_pairs'])
        one_four_pairs = N.array(pn['one_four_pairs'])
        if atom_subset is None:
            atom_subset = N.array([], N.Int)
        else:
            atom_subset = N.array(pn['atom_subset'])
        nbinfo = [excluded_pairs, one_four_pairs, atom_subset]

        if universe.is_periodic:
            try:
                shape = universe.basisVectors()
            except AttributeError:
                raise ValueError("Multipole method implemented only " +
                                 "for orthorhombic universes.")
            shape = N.array(shape, N.Float)
        else:
            shape = N.zeros((), N.Float)

        from MMTK_forcefield import EsMPTerm
        ev = EsMPTerm(universe._spec, shape, nbinfo,
                      pe['charge'], pe['one_four_factor'],
                      pe['spatial_decomposition_levels'],
                      pe['multipole_expansion_terms'],
                      pe['use_fft'],
                      pe['fft_blocking_factor'],
                      pe['macroscopic_expansion_terms'],
                      pe['multipole_acceptance'])
        return [ev]

    # the following methods must be overridden by derived classes
    def _charge(self, o, a, global_data):
        raise AttributeError
