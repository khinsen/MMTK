# Example for a forcefield implementation in Cython.


#
# Get all the required declarations
#
include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

#
# The force field term implementation.
# The rules:
#
# - The class must inherit from EnergyTerm.
#
# - EnergyTerm.__init__() must be called with the arguments
#   shown here. The third argument is the name of the EnergyTerm
#   object, the fourth a tuple of the names of all the terms it
#   implements (one object can implement several terms).
#   The assignment to self.eval_func is essential, without it
#   any energy evaluation will crash.
#
# - The function "evaluate" must have exactly the parameter
#   list from this example.
#
cdef class ElectricFieldTerm(EnergyTerm):

    cdef ArrayType charges
    cdef vector3 strength

    def __init__(self, universe, charges, strength):
        EnergyTerm.__init__(self, universe,
                            "electric_field", ("electric_field",))
        self.eval_func = <void *>ElectricFieldTerm.evaluate
        self.charges = charges
        self.strength[0] = strength[0]
        self.strength[1] = strength[1]
        self.strength[2] = strength[2]

    # The function evaluate is called for every single energy
    # evaluation and should therefore be optimized for speed.
    # Its first argument is the global energy evaluator object,
    # which is needed only for parallelized energy terms.
    # The second argument is a C structure that contains all the
    # input data, in particular the particle configuration.
    # The third argument is a C structure that contains the
    # energy term fields and gradient arrays for storing the results.
    # For details, see MMTK_forcefield.pxi.
    cdef void evaluate(self, EnergyEvaluator eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef double *q
        cdef double e
        cdef int natoms, i

        coordinates = <vector3 *>input.coordinates.data
        natoms = input.coordinates.dimensions[0]
        q = <double *>self.charges.data

        energy.energy_terms[self.index] = 0.
        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

        for i from 0 <= i < natoms:
            e = q[i]*(self.strength[0]*coordinates[i][0]
                      + self.strength[1]*coordinates[i][1]
                      + self.strength[2]*coordinates[i][2])
            energy.energy_terms[self.index] = \
                                energy.energy_terms[self.index] + e
            energy.energy_terms[self.virial_index] = \
                                energy.energy_terms[self.virial_index] - e
            if energy.gradients != NULL:
                gradients[i][0] += q[i]*self.strength[0]
                gradients[i][1] += q[i]*self.strength[1]
                gradients[i][2] += q[i]*self.strength[2]
