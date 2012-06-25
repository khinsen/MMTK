# Center-of-mass restraints
#
# Written by Konrad Hinsen

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

import numpy as N
cimport numpy as N

#
# Trap potential (harmonic restraint to a fixed point in space)
#
cdef class HarmonicTrapTerm(EnergyTerm):

    cdef N.ndarray atom_indices, weights
    cdef double ref_x, ref_y, ref_z
    cdef double k

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices,
                 N.ndarray[N.float64_t] weights,
                 reference, force_constant):
        cdef N.ndarray[N.float64_t] ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_trap", ("harmonic_trap",))
        self.eval_func = <void *>HarmonicTrapTerm.evaluate
        self.atom_indices = atom_indices
        self.weights = weights
        ref_array = reference.array
        self.ref_x = (<double *>ref_array.data)[0]
        self.ref_y = (<double *>ref_array.data)[1]
        self.ref_z = (<double *>ref_array.data)[2]
        self.k = force_constant

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef N.ndarray[N.int_t] atom_indices
        cdef N.ndarray[N.float_t] weights
        cdef N.ndarray[N.float_t, ndim=4] fc
        cdef double cx, cy, cz, dx, dy, dz, tweight, kw
        cdef int i, j, n, offset

        coordinates = <vector3 *>input.coordinates.data
        atom_indices = self.atom_indices
        weights = self.weights

        tweight = 0.
        cx = cy = cz = 0.
        for i in atom_indices:
            tweight += weights[i]
            cx += weights[i]*coordinates[i][0]
            cy += weights[i]*coordinates[i][1]
            cz += weights[i]*coordinates[i][2]
        cx /= tweight
        cy /= tweight
        cz /= tweight

        dx = cx - self.ref_x
        dy = cy - self.ref_y
        dz = cz - self.ref_z

        energy.energy_terms[self.index] = 0.5*self.k*(dx*dx + dy*dy + dz*dz)
        energy.energy_terms[self.virial_index] -= self.k*(dx*dx + dy*dy + dz*dz)

        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
            kw = self.k/tweight
            for i in atom_indices:
                gradients[i][0] += kw*dx*weights[i]
                gradients[i][1] += kw*dy*weights[i]
                gradients[i][2] += kw*dz*weights[i]

        if energy.force_constants != NULL:
            fc = <object>energy.force_constants
            kw = self.k/tweight/tweight
            for i in atom_indices:
                for j in atom_indices:
                    fc [i, 0, j, 0] += kw*weights[i]*weights[j]
                    fc [i, 1, j, 1] += kw*weights[i]*weights[j]
                    fc [i, 2, j, 2] += kw*weights[i]*weights[j]
