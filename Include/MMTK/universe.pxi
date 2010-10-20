cdef extern from "MMTK/universe.h": 

    void import_MMTK_universe()

    ctypedef void distance_fn(vector3 d, vector3 r1, vector3 r2,
                              double *data) nogil
    ctypedef void correction_fn(vector3 *x, int natoms, double *data) nogil
    ctypedef double volume_fn(double scale_factor, double *data) nogil
    ctypedef void box_fn(vector3 *x, vector3 *b, int n, double *data,
                         int to_box) nogil
    ctypedef void bounding_box_fn(vector3 *box1, vector3 *box2, vector3 *x,
                                  int n, double *data) nogil

    ctypedef struct PyUniverseSpecObject:
        PyArrayObject *geometry
        double *geometry_data
        distance_fn *distance_function
        correction_fn *correction_function
        volume_fn *volume_function
        box_fn *box_function
        box_fn *trajectory_function
        bounding_box_fn *bounding_box_function
        int is_periodic
        int is_orthogonal
        int geometry_data_length

    cdef int PyUniverseSpec_StateLock(PyUniverseSpecObject *universe,
                                      int action) nogil

    cdef void **PyUniverse_API
    
import_MMTK_universe()
