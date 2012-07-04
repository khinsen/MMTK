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

    cdef N.ndarray atom_indices, masses
    cdef vector3 ref
    cdef double k

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices,
                 N.ndarray[N.float64_t] masses,
                 reference, force_constant):
        cdef N.ndarray[N.float64_t] ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_cm_trap", ("harmonic_cm_trap",))
        self.eval_func = <void *>HarmonicCMTrapTerm.evaluate
        self.atom_indices = atom_indices
        self.masses = masses
        ref_array = reference.array
        for i in range(3):
            self.ref[i] = (<double *>ref_array.data)[i]
        self.k = force_constant

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef N.ndarray[N.int_t] atom_indices
        cdef N.ndarray[N.float_t] masses
        cdef N.ndarray[N.float_t, ndim=4] fc
        cdef vector3 cm, d
        cdef double m, lsq, kw
        cdef int i, j, n, offset

        coordinates = <vector3 *>input.coordinates.data
        atom_indices = self.atom_indices
        masses = self.masses

        # This code will never be run for a periodic universe, so
        # it's safe to ignore periodic boundary conditions.
        m = 0.
        cm[0] = cm[1] = cm[2] = 0.
        for i in atom_indices:
            m += masses[i]
            for j in range(3):
                cm[j] += masses[i]*coordinates[i][j]
        for j in range(3):
            cm[j] /= m

        for j in range(3):
            d[j] = cm[j]-self.ref[j]
        lsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2]

        energy.energy_terms[self.index] = self.k*lsq
        energy.energy_terms[self.virial_index] -= 2.*self.k*lsq

        if energy.gradients != NULL:
            gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
            kw = 2.*self.k/m
            for i in atom_indices:
                for j in range(3):
                    gradients[i][j] += kw*d[j]*masses[i]

        if energy.force_constants != NULL:
            fc = <object>energy.force_constants
            kw = 2.*self.k/(m*m)
            for i in atom_indices:
                for j in atom_indices:
                    for n in range(3):
                        fc [i, n, j, n] += kw*masses[i]*masses[j]


#
# Harmonic potential between two centers of mass
#
cdef class HarmonicCMDistanceTerm(EnergyTerm):

    cdef N.ndarray atom_indices_1, atom_indices_2, masses
    cdef double k, l0

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices_1,
                 N.ndarray[N.int_t] atom_indices_2,
                 N.ndarray[N.float64_t] masses,
                 distance,
                 force_constant):
        EnergyTerm.__init__(self, universe,
                            "harmonic_cm_distance", ("harmonic_cm_distance",))
        self.eval_func = <void *>HarmonicCMDistanceTerm.evaluate
        self.atom_indices_1 = atom_indices_1
        self.atom_indices_2 = atom_indices_2
        self.masses = masses
        self.l0 = distance
        self.k = force_constant

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates, *gradients
        cdef N.ndarray[N.int_t] atom_indices_1
        cdef N.ndarray[N.int_t] atom_indices_2
        cdef N.ndarray[N.float_t] masses
        cdef N.ndarray[N.float_t, ndim=4] fc
        cdef tensor3 fcij
        cdef vector3 cm1, cm2, d
        cdef double m1, m2, lsq, l, dl, deriv, w
        cdef int i, j, n, offset

        coordinates = <vector3 *>input.coordinates.data
        atom_indices_1 = self.atom_indices_1
        atom_indices_2 = self.atom_indices_2
        masses = self.masses

        # This code will never be run for a periodic universe, so
        # it's safe to ignore periodic boundary conditions.
        m1 = 0.
        cm1[0] = cm1[1] = cm1[2] = 0.
        for i in atom_indices_1:
            m1 += masses[i]
            for j in range(3):
                cm1[j] += masses[i]*coordinates[i][j]
        for j in range(3):
            cm1[j] /= m1

        m2 = 0.
        cm2[0] = cm2[1] = cm2[2] = 0.
        for i in atom_indices_2:
            m2 += masses[i]
            for j in range(3):
                cm2[j] += masses[i]*coordinates[i][j]
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
                    gradients[i][j] -= f1*d[j]*masses[i]
            for i in atom_indices_2:
                for j in range(3):
                    gradients[i][j] += f2*d[j]*masses[i]

        if energy.force_constants != NULL:
            fc = <object>energy.force_constants
            for m in range(3):
                for n in range(3):
                    fcij[m][n] = (2.*self.k-deriv)*d[m]*d[n]/lsq
                fcij[m][m] += deriv
            for i in atom_indices_1:
                for j in atom_indices_1:
                    w = masses[i]*masses[j]/(m1*m1)
                    for m in range(3):
                        for n in range(3):
                            fc [i, m, j, n] += w*fcij[m][n]
            for i in atom_indices_2:
                for j in atom_indices_2:
                    w = masses[i]*masses[j]/(m2*m2)
                    for m in range(3):
                        for n in range(3):
                            fc [i, m, j, n] += w*fcij[m][n]
            for i in atom_indices_1:
                for j in atom_indices_2:
                    w = masses[i]*masses[j]/(m1*m2)
                    for m in range(3):
                        for n in range(3):
                            if i < j:
                                fc [i, m, j, n] -= w*fcij[m][n]
                            else:
                                fc [j, m, i, n] -= w*fcij[m][n]
