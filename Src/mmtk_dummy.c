/* This pseudo-module is not used at all. It only exists to make
* distutils compile lapack_dlamc.c with special compiler options. */

#include "Python.h"

extern double dlamch_(char *);

static struct PyMethodDef mmtk_dummy_module_methods[] = {
{ NULL,NULL,0}
};

void initmmtk_dummy()
{ PyObject *m,*d;
m = Py_InitModule("mmtk_dummy", mmtk_dummy_module_methods);
dlamch_("P");
}
