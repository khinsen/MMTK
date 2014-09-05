cdef extern from "MMTK/arrayobject.h":

    void import_array()

    cdef enum Pyarray_TYPES:
        PyArray_CHAR, PyArray_UBYTE, PyArray_SBYTE,
        PyArray_SHORT, PyArray_USHORT, 
        PyArray_INT, PyArray_UINT, 
        PyArray_LONG,
        PyArray_FLOAT, PyArray_DOUBLE, 
        PyArray_CFLOAT, PyArray_CDOUBLE,
        PyArray_OBJECT,
        PyArray_NTYPES, PyArray_NOTYPE

    struct PyArray_Descr: 
        int type_num, elsize 
        char type 

    ctypedef struct PyArrayObject:
        char *data 
        int nd 
        int *dimensions 
        int *strides 
        PyObject *base 
        PyArray_Descr *descr 
        int flags

    ctypedef class Scientific.N.ArrayType [object PyArrayObject]: 
        cdef char *data 
        cdef int nd 
        cdef int *dimensions 
        cdef int *strides 
        cdef object base 
        cdef PyArray_Descr *descr 
        cdef int flags

    object PyArray_FromDims(int n_dimensions, int dimensions[], int item_type)

import_array()
