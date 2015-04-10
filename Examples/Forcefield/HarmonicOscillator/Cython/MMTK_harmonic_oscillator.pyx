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
#   list given in this example.
#
cdef class HarmonicOscillatorTerm(EnergyTerm):

    cdef int atom_index
    cdef double ref_x, ref_y, ref_z
    cdef double k

    def __init__(self, universe, atom_index, reference, force_constant):
        cdef ArrayType ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_oscillator", ("harmonic_oscillator",))
        self.eval_func = <void *>HarmonicOscillatorTerm.evaluate
        self.atom_index = atom_index
        ref_array = reference.array
        self.ref_x = (<double *>ref_array.data)[0]
        self.ref_y = (<double *>ref_array.data)[1]
        self.ref_z = (<double *>ref_array.data)[2]
        self.k = force_constant

    # The function evaluate is called for every single energy
    # evaluation and should therefore be optimized for speed.
    # Its first argument is the global energy evaluator object,
    # which is needed only for parallelized energy terms.
    # The second argument is a C structure that contains all the
    # input data, in particular the particle configuration.
    # The third argument is a C structure that contains the
    # energy term fields and gradient arrays for storing the results.
    # For details, see MMTK_forcefield.pxi.
    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates
        cdef vector3 *gradients
        cdef double *fc
        cdef double dx, dy, dz
        cdef int n, offset
        coordinates = <vector3 *>input.coordinates.data
        dx = coordinates[self.atom_index][0] - self.ref_x
        dy = coordinates[self.atom_index][1] - self.ref_y
        dz = coordinates[self.atom_index][2] - self.ref_z
        energy.energy_terms[self.index] = 0.5*self.k*(dx*dx + dy*dy + dz*dz)
        energy.energy_terms[self.virial_index] -= self.k*(dx*dx + dy*dy + dz*dz)
        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
            gradients[self.atom_index][0] += self.k*dx
            gradients[self.atom_index][1] += self.k*dy
            gradients[self.atom_index][2] += self.k*dz
        if energy.force_constants != NULL:
            fc = <double *>(<PyArrayObject *> energy.force_constants).data
            n = (<PyArrayObject *> energy.force_constants).dimensions[0]
            offset = (9*n+3)*self.atom_index
            fc[offset + 3*n*0 + 0] += self.k
            fc[offset + 3*n*1 + 1] += self.k
            fc[offset + 3*n*2 + 2] += self.k
