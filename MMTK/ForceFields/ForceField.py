# This module implements classes that represent general force fields
# and force field evaluators.
#
# Written by Konrad Hinsen
#

from MMTK import Environment, ParticleProperties, Units, Universe, Utility
from Scientific import N
import copy, itertools, operator
from MMTK_energy_term import PyEnergyTerm as EnergyTerm
from MMTK_forcefield import HarmonicDistanceTerm

# Class definitions

#
# The base class ForceField contains common operations for all force fields
#
class ForceField(object):

    def __init__(self, name):
        self.name = name
        self.type = None

    def __reduce__(self):
        return type(self), self.arguments

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        # must be defined by derived classes
        raise NotImplementedError

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # must be defined by derived classes
        raise NotImplementedError

    def __add__(self, other):
        return CompoundForceField(self, other)

    def getOptions(self, object, options):
        attr = self.type + '_options'
        if hasattr(object, attr):
            for key, value in getattr(object, attr).items():
                options[key] = value
        attr = self.name + '_options'
        if hasattr(object, attr):
            for key, value in getattr(object, attr).items():
                options[key] = value

    def ready(self, global_data):
        return True

    def declareDependencies(self, global_data):
        pass

    def bondedForceFields(self):
        return []

    def bondLengthDatabase(self, universe):
        return None

    def description(self):
        return self.__class__.__module__ + '.' + \
               self.__class__.__name__ + `self.arguments`

    def getAtomParameterIndices(self, atoms):
        universe = None
        indices = []
        for a in atoms:
            if isinstance(a, int):
                indices.append(a)
            else:
                if universe is None:
                    universe = a.universe()
                    universe.configuration()
                else:
                    if universe is not a.universe():
                        raise ValueError("Atom parameters %s are not all "
                                         "in the same universe" % repr(atoms))
                indices.append(a.index)
        return tuple(indices)

    def supportsPathIntegrals(self):
        return False

    def beadOffsetsAndFactor(self, atom_indices, global_data):
        # Maximum number of beads in any atom.
        nbeads = global_data.get('nbeads')
        bead_data = global_data.get('bead_data')
        # Bead counts for the requested atoms.
        nb = [bead_data[i, 1] for i in atom_indices]

        return Utility.beadOffsetsAndFactor(nb, nbmax=nbeads)

#
# A CompoundForceField represents the sum of its component force fields
#
class CompoundForceField(ForceField):

    def __init__(self, *args):
        self.fflist = list(args)
        self.name = '*'
        self.type = 'compound'
        self.arguments = args

    def declareDependencies(self, global_data):
        for ff in self.fflist:
            ff.declareDependencies(global_data)

    def flatFFList(self):
        flat_list = []
        for ff in self.fflist:
            if isinstance(ff, CompoundForceField):
                flat_list.extend(ff.flatFFList())
            else:
                flat_list.append(ff)
        return flat_list

    # def evaluatorParameters(self, system, subset1, subset2, global_data):
    #     self.declareDependencies()
    #     parameters = {}
    #     remaining = self.flatFFList()
    #     while remaining:
    #         done = []
    #         for ff in remaining:
    #             if ff.ready(global_data):
    #                 params = ff.evaluatorParameters(system, subset1, subset2,
    #                                                 global_data)
    #                 if not isinstance(params, dict):
    #                     raise ValueError("evaluator parameters are not a dict")
    #                 for key, value in params.items():
    #                     parameters[key] = _combine(value,
    #                                                parameters.get(key, None))
    #                 done.append(ff)
    #         if not done:
    #             raise TypeError("Cyclic force field dependence")
    #         for ff in done:
    #             remaining.remove(ff)
    #     return parameters

    # def evaluatorTerms(self, system, subset1, subset2, global_data):
    #     self.declareDependencies()
    #     eval_objects = []
    #     remaining = self.flatFFList()
    #     while remaining:
    #         done = []
    #         for ff in remaining:
    #             if ff.ready(global_data):
    #                 terms = ff.evaluatorTerms(system, subset1, subset2,
    #                                           global_data)
    #                 if not isinstance(terms, list):
    #                     raise ValueError("evaluator term list not a list")
    #                 eval_objects.extend(terms)
    #                 done.append(ff)
    #         if not done:
    #             raise TypeError("Cyclic force field dependence")
    #         for ff in done:
    #             remaining.remove(ff)
    #     return eval_objects

    def evaluatorParameters(self, system, subset1, subset2, global_data):
        parameters = {}
        def add_item(params):
            if not isinstance(params, dict):
                raise ValueError("evaluator parameters are not a dict")
            for key, value in params.items():
                parameters[key] = _combine(value,
                                           parameters.get(key, None))
        self._compoundEvaluator(system, subset1, subset2, global_data,
                                'evaluatorParameters', add_item)
        return parameters

    def evaluatorTerms(self, system, subset1, subset2, global_data):
        eval_objects = []
        def add_item(terms):
            if not isinstance(terms, list):
                raise ValueError("evaluator term list not a list")
            eval_objects.extend(terms)
        self._compoundEvaluator(system, subset1, subset2, global_data,
                                'evaluatorTerms', add_item)
        return eval_objects

    def _compoundEvaluator(self, system, subset1, subset2, global_data,
                           method, add_item):
        self.declareDependencies(global_data)
        remaining = self.flatFFList()
        while remaining:
            done = []
            for ff in remaining:
                if ff.ready(global_data):
                    item = getattr(ff, method)(system, subset1, subset2,
                                               global_data)
                    add_item(item)
                    done.append(ff)
            if not done:
                raise TypeError("Cyclic force field dependence")
            for ff in done:
                remaining.remove(ff)

    def bondedForceFields(self):
        return reduce(operator.add,
                      [f.bondedForceFields() for f in self.fflist])

    def description(self):
        return '+'.join([f.description() for f in self.fflist])

    def supportsPathIntegrals(self):
        for ff in self.fflist:
            if not ff.supportsPathIntegrals():
                return False
        return True

def _combine(params1, params2):
    if params1 is None:
        return params2
    if params2 is None:
        return params1
    if isinstance(params1, list) and isinstance(params2, list):
        return params1 + params2
    if params1 == params2:
        return params1
    raise TypeError, "_combine not implemented yet"

#
# This class serves to define data containers used in
# force field initialization.
#
class ForceFieldData(object):

    def __init__(self):
        self.dict = {}

    def add(self, tag, value):
        try:
            self.dict[tag].append(value)
        except KeyError:
            self.dict[tag] = [value]

    def get(self, tag):
        return self.dict.get(tag, [])

    def set(self, tag, value):
        self.dict[tag] = value


# Type check functions

def isForceField(x):
    return isinstance(x, ForceField)

def isCompoundForceField(x):
    return isinstance(x, CompoundForceField)

#
# Force field evaluator support functions.
# These functions are meant to be used with force field evaluators
# implemented in Python using automatic derivatives. They put the
# first and second derivative information into the right places.
#
def addToGradients(coordinates, indices, vectors):
    for index, vector in zip(indices, vectors):
        coordinates[index] += vector.array

def addToForceConstants(total_fc, indices, small_fc):
    indices = zip(indices, range(len(indices)))
    for i1, i2 in itertools.chain(itertools.izip(indices, indices),
                                  Utility.pairs(indices)):
        ii1, jj1 = i1
        ii2, jj2 = i2
        if ii1 > ii2:
            ii1, ii2 = ii2, ii1
            jj1, jj2 = jj2, jj1
        jj1 = 3*jj1
        jj2 = 3*jj2
        total_fc[ii1,:,ii2,:] += small_fc[jj1:jj1+3, jj2:jj2+3]

#
# High-level energy evaluator (i.e. the Python interface)
#
class _EnergyEvaluator(object):

    def __init__(self, universe, force_field, subset1=None, subset2=None):
        if not Universe.isUniverse(universe):
            raise TypeError("energy evaluator defined only for universes")
        self.universe = universe
        self.universe_version = self.universe._version
        self.ff = force_field
        self.configuration = self.universe.configuration()
        self.global_data = ForceFieldData()
        self.global_data.set('universe', universe)
        self.subset1 = subset1
        self.subset2 = subset2 if subset2 is not None else subset1

        self.spring_parameters = []

        if self.subset1 is None:
            pi_atoms = [a for a in self.universe.atomIterator() if a.numberOfBeads() > 1]
        else:
            pi_atoms = list(set(a for a in self.subset1.atomIterator() if a.numberOfBeads() > 1)
                            | set(a for a in self.subset2.atomIterator() if a.numberOfBeads() > 1))

        if pi_atoms:
            if not self.ff.supportsPathIntegrals():
                raise ValueError("Some force field term doesn't handle path integrals")
            nbeads = self.universe.maxNumberOfBeads()
            pi_environments = self.universe.environmentObjectList(
                                                Environment.PathIntegrals)
            if len(pi_environments) == 1:
                beta = pi_environments[0].beta
            else:
                raise ValueError('exactly one path integral environment required')
            if pi_environments[0].include_spring_terms:
                for a in pi_atoms:
                    nb = a.numberOfBeads()
                    k = float(nbeads*nb)*a.mass() / (beta*beta*Units.hbar*Units.hbar*2.)
                    for b in range(nb):
                        self.spring_parameters.append((a.index+b, a.index+(b+1)%nb, 0., k))
        else:
            nbeads = 1

        self.global_data.set('nbeads', nbeads)
        bead_data = N.zeros((self.universe.numberOfPoints(), 2), N.Int32)
        for b in universe.beadIterator():
            bead_data[b.index, 0] = b.bead_number
            bead_data[b.index, 1] = b.atom.numberOfBeads()
        self.global_data.set('bead_data', bead_data)


class EnergyEvaluatorParameters(_EnergyEvaluator):

    def __init__(self, universe, force_field, subset1=None, subset2=None):
        _EnergyEvaluator.__init__(self, universe, force_field, subset1, subset2)

        self.params = self.ff.evaluatorParameters(self.universe,
                                                  self.subset1, self.subset2,
                                                  self.global_data)
        key = 'harmonic_distance_term'
        self.params[key] = _combine(self.spring_parameters,
                                    self.params.get(key, None))


class EnergyEvaluator(_EnergyEvaluator):

    def __init__(self, universe, force_field, subset1=None, subset2=None,
                 threads=None, mpi_communicator=None):
        _EnergyEvaluator.__init__(self, universe, force_field, subset1, subset2)

        terms = []
        if self.spring_parameters:
            indices = [x[:2] for x in self.spring_parameters]
            parameters = [x[2:] for x in self.spring_parameters]
            terms.append(HarmonicDistanceTerm(universe._spec,
                                              N.array(indices),
                                              N.array(parameters),
                                              "path integral spring"))
        terms.extend(self.ff.evaluatorTerms(self.universe,
                                            self.subset1, self.subset2,
                                            self.global_data))

        from MMTK_forcefield import Evaluator
        if threads is None:
            import MMTK.ForceFields
            threads = MMTK.ForceFields.default_energy_threads
        self.evaluator = Evaluator(N.array(terms),
                                   1./self.global_data.get('nbeads'),
                                   threads, mpi_communicator)

    def checkUniverseVersion(self):
        if self.universe_version != self.universe._version:
            raise ValueError('the universe has been modified')

    def CEvaluator(self):
        return self.evaluator

    def __call__(self, gradients = None, force_constants = None,
                 small_change=False):
        self.checkUniverseVersion()
        args = [self.configuration.array]

        if ParticleProperties.isParticleProperty(gradients):
            args.append(gradients.array)
        elif isinstance(gradients, N.array_type):
            gradients = \
                 ParticleProperties.ParticleVector(self.universe, gradients)
            args.append(gradients.array)
        elif gradients:
            gradients = ParticleProperties.ParticleVector(self.universe, None)
            args.append(gradients.array)
        else:
            args.append(None)

        if ParticleProperties.isParticleProperty(force_constants):
            args.append(force_constants.array)
        elif isinstance(force_constants, N.array_type):
            force_constants = \
                ParticleProperties.SymmetricPairTensor(self.universe,
                                                       force_constants)
            args.append(force_constants.array)
        else:
            sparse_type = None
            try:
                from MMTK_forcefield import SparseForceConstants
                sparse_type = type(SparseForceConstants(2, 2))
            except ImportError: pass
            if isinstance(force_constants, sparse_type):
                args.append(force_constants)
            elif force_constants:
                force_constants = \
                    ParticleProperties.SymmetricPairTensor(self.universe, None)
                args.append(force_constants.array)
            else:
                args.append(None)

        args.append(small_change)
        self.universe.acquireReadStateLock()
        try:
            energy = apply(self.evaluator, tuple(args))
        finally:
            self.universe.releaseReadStateLock()
        if force_constants:
            return energy, gradients, force_constants
        elif gradients:
            return energy, gradients
        else:
            return energy

    def lastEnergyTerms(self):
        dict = {}
        values = self.evaluator.last_energy_values
        i = 0
        for term_object in self.evaluator:
            for name in term_object.term_names:
                dict[name] = dict.get(name, 0.) + values[i]
                i = i + 1
        for name, value in dict.items():
            index = name.find('/')
            if index >= 0:
                category = name[:index]
                dict[category] = dict.get(category, 0.) + value
        for name, value in dict.items():
            if '/' in name and value == 0.:
                del dict[name]
        return dict

    def lastVirial(self):
        return self.evaluator.last_energy_values[-1]
