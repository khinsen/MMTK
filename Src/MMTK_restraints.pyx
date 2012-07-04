# Center-of-mass restraints
#
# Written by Konrad Hinsen

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

cdef extern from "math.h":
    double sqrt(double)

import numpy as N
cimport numpy as N

#
# Trap potential (harmonic restraint to a fixed point in space)
#
cdef class HarmonicCMTrapTerm(EnergyTerm):

    cdef N.ndarray atom_indices, weights
    cdef double ref_x, ref_y, ref_z
    cdef double k

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices,
                 N.ndarray[N.float64_t] weights,
                 reference, force_constant):
        cdef N.ndarray[N.float64_t] ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_cm_trap", ("harmonic_cm_trap",))
        self.eval_func = <void *>HarmonicCMTrapTerm.evaluate
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

        energy.energy_terms[self.index] = self.k*(dx*dx + dy*dy + dz*dz)
        energy.energy_terms[self.virial_index] -= 2.*self.k * \
                                                       (dx*dx + dy*dy + dz*dz)

        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
            kw = 2.*self.k/tweight
            for i in atom_indices:
                gradients[i][0] += kw*dx*weights[i]
                gradients[i][1] += kw*dy*weights[i]
                gradients[i][2] += kw*dz*weights[i]

        if energy.force_constants != NULL:
            fc = <object>energy.force_constants
            kw = 2.*self.k/tweight/tweight
            for i in atom_indices:
                for j in atom_indices:
                    fc [i, 0, j, 0] += kw*weights[i]*weights[j]
                    fc [i, 1, j, 1] += kw*weights[i]*weights[j]
                    fc [i, 2, j, 2] += kw*weights[i]*weights[j]


#
# Harmonic potential between two centers of mass
#
cdef class HarmonicCMDistanceTerm(EnergyTerm):

    cdef N.ndarray atom_indices_1, atom_indices_2, weights
    cdef double k, l0

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices_1,
                 N.ndarray[N.int_t] atom_indices_2,
                 N.ndarray[N.float64_t] weights,
                 distance,
                 force_constant):
        EnergyTerm.__init__(self, universe,
                            "harmonic_cm_distance", ("harmonic_cm_distance",))
        self.eval_func = <void *>HarmonicCMDistanceTerm.evaluate
        self.atom_indices_1 = atom_indices_1
        self.atom_indices_2 = atom_indices_2
        self.weights = weights
        self.l0 = distance
        self.k = force_constant

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef N.ndarray[N.int_t] atom_indices_1
        cdef N.ndarray[N.int_t] atom_indices_2
        cdef N.ndarray[N.float_t] weights
        cdef N.ndarray[N.float_t, ndim=4] fc
        cdef tensor3 fcij
        cdef vector3 cm1, cm2, d
        cdef double m1, m2, lsq, l, dl, deriv, w
        cdef int i, j, n, offset

        coordinates = <vector3 *>input.coordinates.data
        atom_indices_1 = self.atom_indices_1
        atom_indices_2 = self.atom_indices_2
        weights = self.weights

        # This code will never be run for a periodic universe, so
        # it's safe to ignore periodic boundary conditions.
        m1 = 0.
        cm1[0] = cm1[1] = cm1[2] = 0.
        for i in atom_indices_1:
            m1 += weights[i]
            for j in range(3):
                cm1[j] += weights[i]*coordinates[i][j]
        for j in range(3):
            cm1[j] /= m1

        m2 = 0.
        cm2[0] = cm2[1] = cm2[2] = 0.
        for i in atom_indices_2:
            m2 += weights[i]
            for j in range(3):
                cm2[j] += weights[i]*coordinates[i][j]
        for j in range(3):
            cm2[j] /= m2

        for j in range(3):
            d[j] = cm2[j]-cm1[j]
        lsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2]
        l = sqrt(lsq)
        dl = l - self.l0
        energy.energy_terms[self.index] = self.k*dl*dl
        energy.energy_terms[self.virial_index] -= 2.*self.k*dl*dl

        if energy.gradients != NULL or energy.force_constants != NULL:
            deriv = 0. if l == 0. else 2.*self.k*dl/l

        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
            deriv = 0. if l == 0. else 2.*self.k*dl/l
            f1 = self.k*deriv/m1
            f2 = self.k*deriv/m2
            for i in atom_indices_1:
                for j in range(3):
                    gradients[i][j] -= f1*d[j]*weights[i]
            for i in atom_indices_2:
                for j in range(3):
                    gradients[i][j] += f2*d[j]*weights[i]

        if energy.force_constants != NULL:
            fc = <object>energy.force_constants
            for m in range(3):
                for n in range(3):
                    fcij[m][n] = (2.*self.k-deriv)*d[m]*d[n]/lsq
                fcij[m][m] += deriv
            for i in atom_indices_1:
                for j in atom_indices_1:
                    w = weights[i]*weights[j]/(m1*m1)
                    for m in range(3):
                        for n in range(3):
                            fc [i, m, j, n] += w*fcij[m][n]
            for i in atom_indices_2:
                for j in atom_indices_2:
                    w = weights[i]*weights[j]/(m2*m2)
                    for m in range(3):
                        for n in range(3):
                            fc [i, m, j, n] += w*fcij[m][n]
            for i in atom_indices_1:
                for j in atom_indices_2:
                    w = weights[i]*weights[j]/(m1*m2)
                    for m in range(3):
                        for n in range(3):
                            if i < j:
                                fc [i, m, j, n] -= w*fcij[m][n]
                            else:
                                fc [j, m, i, n] -= w*fcij[m][n]
