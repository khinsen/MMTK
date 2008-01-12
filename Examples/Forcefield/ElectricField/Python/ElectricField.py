# Electric field term representing a uniform external electric field.
# It can be added to any other force field.

from MMTK import ParticleScalar, ParticleVector, SymmetricPairTensor
from MMTK.ForceFields.ForceField import ForceField, EnergyTerm

#
# Each force field term requires two classes. One represents the
# abstract force field (ElectricField in this example), it knows
# nothing about any molecular system to which it might be applied.
# The other one (ElectricFieldTerm) stores all the parameters required
# for energy evaluation, and provides the actual code.
#
# When a force field is evaluated for the first time for a given universe,
# the method evaluatorTerms() of the force field is called. It creates
# the EnergyTerm objects. The method evaluate() of the energy term object
# is called for every energy evaluation. Consequently, as much of the
# computation as possible should be done outside of this method.
#

class ElectricFieldTerm(EnergyTerm):

    # The __init__ method only remembers parameters. Note that
    # EnergyTerm.__init__ takes care of storing the name and the
    # universe object.
    def __init__(self, universe, strength, charges):
        EnergyTerm.__init__(self, 'electric field', universe)
        self.strength = strength
        self.charges = charges

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    def evaluate(self, configuration, do_gradients, do_force_constants):
        results = {}
        energy = (self.charges*(configuration*self.strength)) \
                 .sumOverParticles()
        results['energy'] = energy
        results['virial'] = -energy
        if do_gradients:
            gradients = ParticleVector(self.universe)
            for atom in self.universe.atomList():
                gradients[atom] = self.charges[atom]*self.strength
            results['gradients'] = gradients
        if do_force_constants:
            # The force constants are zero -> nothing to calculate.
            results['force_constants'] = SymmetricPairTensor(self.universe)
        else:
            force_constants = None            
        return results


class ElectricField(ForceField):

    """Uniform electric field

    Constructor: ElectricField(|strength|, |charge_property|='amber_charge')

    Arguments:

    |strength| -- a vector object representing the electric field vector

    |charge_property| -- a string indicating the name of the atomic
                         property in the database that is used to retrieve
                         the atomic charges. The default is 'amber_charge',
                         the charge property for the Amber94 force field.
    """

    def __init__(self, strength, charge_property='amber_charge'):
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (strength, charge_property)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'electric_field')
        # Store the parameters for later use.
        self.strength = strength
        self.charge_property = charge_property

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return 1

    # The following method is called by the energy evaluation engine
    # to obtain a list of the evaluator objects
    # that handle the calculations.
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # The energy for subsets is defined as consisting only
        # of interactions within that subset, so the contribution
        # of an external field is zero. Therefore we just return
        # an empty list of energy terms.
        if subset1 is not None or subset2 is not None:
            return []
        # Collect the charges into an array
        charges = ParticleScalar(universe)
        for o in universe:
            for a in o.atomList():
                charges[a] = o.getAtomProperty(a, self.charge_property)
        # Here we pass all the parameters to
        # the energy term code that handles energy calculations.
        return [ElectricFieldTerm(universe, self.strength, charges)]
