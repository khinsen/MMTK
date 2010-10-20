cdef extern from "MMTK/forcefield.h":

    void import_MMTK_forcefield()

    # This structure combines all the input data required by an
    # energy calculation.
    ctypedef struct energy_spec:
        # particle configuration
        PyArrayObject *coordinates
        # number of particles
        int natoms
        # thread number (only for multi-threaded energy terms)
        int thread_id
        # processor number (only for parallelized energy terms)
        int proc_id
        # slice number (only for parallelized energy terms)
        # A slice is the part of the total that one thread on
        # one processor works on. The number of slices is
        # the number of processors times the number of threads
        # per processor, where one processor is really one node
        # in a parallel system, i.e. a multiprocessor node with
        # shared memory is counted as *one* processor.
        int slice_id
        # number of threads (only for multi-threaded energy terms)
        int nthreads
        # number of processors (or rather nodes)
        int nprocs
        # number of slices
        int nslices
        # Set to 1 if the input configuration is close to the one of
        # the previous call. This is only a hint for choosing the best
        # algorithm.
        int small_change

    ctypedef struct energy_data
    ctypedef int gradient_function(energy_data *energy,
                                   int i, vector3 gradient)
    ctypedef int fc_function(energy_data *energy,
                             int i, int j, tensor3 fc,
                             double r_sq)
    ctypedef struct energy_data:
        # An array or a sparse gradient storage.
        PyObject *gradients
        # A pointer to a function that gets/sets the gradient components,
        # used when gradients is NOT an array.
        gradient_function *gradient_fn
        # An array or a sparse force constant matrix.
        PyObject *force_constants
        # A pointer to a function that gets/sets the force constant components,
        # used when force_constannt is NOT an array.
        fc_function *fc_fn
        # An array of energy term values. Each energy term should
        # access only its alloted slots in there, plus the virial
        # term indicated by self->virial_index
        double *energy_terms
        # The total energy, calculated after all energy terms
        # have been evaluated. Don't touch inside energy term routines
        double energy
        # The virial. The same comment applies. Energy terms should add
        # their contribution to energy_terms[self->virial_index].
        double virial
        # If a term cannot calculate a virial, it should set this
        # variable to 0, in which case the total virial is considered
        # invalid and pressure-related calculations are not possible.
        int virial_available
        # Set this to 1 if there is an error during processing.
        int error

    ctypedef class MMTK_forcefield.EnergyTerm [object PyFFEnergyTermObject]:
        cdef int index
        cdef int virial_index
        cdef void *eval_func
        cdef int nterms
        cdef char *evaluator_name
        cdef char **term_names

    ctypedef struct PyFFEvaluatorObject
    
    ctypedef void ff_eval_function(PyFFEvaluatorObject *evaluator,
                                   energy_data *ed,
                                   N.ndarray[double, ndim=2] configuration,
                                   int small_change) nogil

    ctypedef struct PyFFEvaluatorObject:
        ff_eval_function eval_func
        double *energy_terms
        PyThreadState *tstate_save

import_MMTK_forcefield()

