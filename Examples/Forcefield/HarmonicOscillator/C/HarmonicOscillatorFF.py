# Harmonic potential with respect to a fixed point in space

from MMTK.ForceFields.ForceField import ForceField
from MMTK_harmonic_oscillator import HarmonicOscillatorTerm

class HarmonicOscillatorForceField(ForceField):

    """Harmonic potential with respect to a fixed point in space

    Constructor: HarmonicOscillatorForceField(|atom|, |center|,
                                              |force_constant|)

    Arguments:

    |atom| -- an atom object, specifying the
              atom on which the force field acts

    |center| -- a vector defining the point to which the atom is
                attached by the harmonic potential

    |force_constant| -- the force constant of the harmonic potential
                        (a real number)
    """

    def __init__(self, atom, center, force_constant):
        # We call this for its side effect: it makes sure that all
        # atoms have a correct index assigned to them.
        atom.universe().configuration()
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (atom, center, force_constant)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'harmonic_oscillator')
        # Store the parameters for later use.
        self.atom = atom
        self.center = center
        self.force_constant = force_constant

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
        # The subset evaluation mechanism does not make much sense for
        # this force field, so we just signal an error if someone
        # uses it by accident.
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        # Here we pass all the parameters as "simple" data types to
        # the C code that handles energy calculations.
        return [HarmonicOscillatorTerm(universe._spec,
                                       self.atom.index,
                                       self.center[0], self.center[1],
                                       self.center[2], self.force_constant)]

    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return self.__class__.__name__ + `self.arguments`
