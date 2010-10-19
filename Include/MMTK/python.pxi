cdef extern from "Python.h":

    ctypedef void PyObject

    ctypedef void PyThreadState

    cdef PyThreadState *PyEval_SaveThread()
    cdef void PyEval_RestoreThread(PyThreadState *thread)
