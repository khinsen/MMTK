# Electric field term representing a uniform external electric field.
# can be added to any other force field

from MMTK import ParticleScalar
from MMTK.ForceFields.ForceField import ForceField
from MMTK_electric_field import ElectricFieldTerm

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
    # to obtain a list of the low-level evaluator objects (the C routines)
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
        # Here we pass all the parameters as "simple" data types to
        # the C code that handles energy calculations.
        return [ElectricFieldTerm(universe._spec, charges.array,
                                  self.strength[0], self.strength[1],
                                  self.strength[2])]

    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return self.__class__.__name__ + `self.arguments`
