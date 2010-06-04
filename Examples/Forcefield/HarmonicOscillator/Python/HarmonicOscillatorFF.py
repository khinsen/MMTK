# Harmonic potential with respect to a fixed point in space

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK import ParticleVector, SymmetricPairTensor
from Scientific.Geometry import delta

#
# Each force field term requires two classes. One represents the
# abstract force field (HarmonicOscillator in this example), it knows
# nothing about any molecular system to which it might be applied.
# The other one (HarmonicOscillatorTerm) stores all the parameters required
# for energy evaluation, and provides the actual code.
#
# When a force field is evaluated for the first time for a given universe,
# the method evaluatorTerms() of the force field is called. It creates
# the EnergyTerm objects. The method evaluate() of the energy term object
# is called for every energy evaluation. Consequently, as much of the
# computation as possible should be done outside of this method.
#

class HarmonicOscillatorTerm(EnergyTerm):

    # The __init__ method only remembers parameters. Note that
    # EnergyTerm.__init__ takes care of storing the name and the
    # universe object.
    def __init__(self, universe, atom_index, center, force_constant):
        EnergyTerm.__init__(self, 'harmonic oscillator', universe)
        self.atom_index = atom_index
        self.center = center
        self.force_constant = force_constant

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    def evaluate(self, configuration, do_gradients, do_force_constants):
        results = {}
        d = configuration[self.atom_index]-self.center
        results['energy'] = 0.5*self.force_constant*(d*d)
        if do_gradients:
            gradients = ParticleVector(self.universe)
            gradients[self.atom_index] = self.force_constant*d
            results['gradients'] = gradients
        if do_force_constants:
            force_constants = SymmetricPairTensor(self.universe)
            force_constants[self.atom_index, self.atom_index] = self.force_constant*delta
            results['force_constants'] = force_constants
        return results

class HarmonicOscillatorForceField(ForceField):

    """
    Harmonic potential with respect to a fixed point in space
    """

    def __init__(self, atom, center, force_constant):
        """
        @param atom: the atom on which the force field acts
        @type atom: L{MMTK.ChemicalObjects.Atom}
        @param center: the point to which the atom is attached by
                       the harmonic potential
        @type center: L{Scientific.Geometry.Vector}
        @param force_constant: the force constant of the harmonic potential
        @type force_constant: C{float}
        """
        # Get the internal index of the atom if the argument is an
        # atom object. It is the index that is stored internally,
        # and when the force field is recreated from a specification
        # in a trajectory file, it is the index that is passed instead
        # of the atom object itself. Calling this method takes care
        # of all the necessary checks and conversions.
        self.atom_index = self.getAtomParameterIndices([atom])[0]
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (self.atom_index, center, force_constant)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'harmonic_oscillator')
        # Store the parameters for later use.
        self.center = center
        self.force_constant = force_constant

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return True

    # The following method is called by the energy evaluation engine
    # to obtain a list of the low-level evaluator objects (the C routines)
    # that handle the calculations.
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # The subset evaluation mechanism does not make much sense for
        # this force field, so we just signal an error if someone
        # uses it by accident.
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        # Here we pass all the parameters to the code
        # that handles energy calculations.
        return [HarmonicOscillatorTerm(universe,
                                       self.atom_index,
                                       self.center,
                                       self.force_constant)]
