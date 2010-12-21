# Electric field term representing a uniform external electric field.
# can be added to any other force field

from MMTK import ParticleScalar
from MMTK.ForceFields.ForceField import ForceField
from MMTK_electric_field import ElectricFieldTerm

class ElectricField(ForceField):

    """
    Uniform electric field
    """

    def __init__(self, strength, charge_property='amber_charge'):
        """
        @param strength: the electric field vector
        @type strength: L{Scientific.Geometry.Vector}
        @charge_property: the name of the atomic property in the database
                          that is used to retrieve the atomic charges.
                          The default is 'amber_charge', the charge property
                          for Amber94 and Amber99.
        @type charge_property: C{str}
        """
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
        return True

    # For a force field that supports path integrals (i.e. calculates
    # the right terms for all beads etc.), the method supportsPathIntegrals,
    # whose default implementation returns False, must be overridden.
    def supportsPathIntegrals(self):
        return True

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
                c = o.getAtomProperty(a, self.charge_property)
                for b in a.beads():
                    charges[b] = c
        # Here we pass all the parameters to
        # the Cython code that handles energy calculations.
        return [ElectricFieldTerm(universe,
                                  charges.array,
                                  self.strength)]
