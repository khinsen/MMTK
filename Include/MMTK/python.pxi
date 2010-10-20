cdef extern from "Python.h":

    ctypedef void PyObject

    ctypedef struct PyThreadState

    cdef PyThreadState *PyEval_SaveThread()
    cdef void PyEval_RestoreThread(PyThreadState *thread)
