cdef extern from "Python.h":

    ctypedef void PyObject
    ctypedef PyObject PyListObject

    ctypedef struct PyThreadState

    cdef PyThreadState *PyEval_SaveThread()
    cdef void PyEval_RestoreThread(PyThreadState *thread)

    cdef void *PyMem_Malloc(ssize_t length)
    cdef void PyMem_Free(void *)

