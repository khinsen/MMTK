# Center-of-mass restraints
#
# Written by Konrad Hinsen

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

#
# Trap potential (harmonic restraint to a fixed point in space)
#
cdef class HarmonicTrapTerm(EnergyTerm):

    cdef int atom_index
    cdef double ref_x, ref_y, ref_z
    cdef double k

    def __init__(self, universe, atom_index, reference, force_constant):
        cdef ArrayType ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_trap", ("harmonic_trap",))
        self.eval_func = <void *>HarmonicTrapTerm.evaluate
        self.atom_index = atom_index
        ref_array = reference.array
        self.ref_x = (<double *>ref_array.data)[0]
        self.ref_y = (<double *>ref_array.data)[1]
        self.ref_z = (<double *>ref_array.data)[2]
        self.k = force_constant

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
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
