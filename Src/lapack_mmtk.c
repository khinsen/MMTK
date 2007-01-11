#include "Python.h"
#include "MMTK/arrayobject.h"

typedef struct { float r, i; } f2c_complex;
typedef struct { double r, i; } f2c_doublecomplex;
typedef long int (*L_fp)();

static PyObject *ErrorObject;

static PyObject *LapackError()
{  if (! ErrorObject)
      ErrorObject = PyString_FromString("LapackError");
   Py_INCREF(ErrorObject);
   return ErrorObject;
}

static PyObject *ErrorReturn(char *mes)
{  if (!ErrorObject)
      ErrorObject = PyString_FromString("LapackError");
   PyErr_SetString(ErrorObject,mes);
   return NULL;
}

#define TRY(E) if( ! (E)) return NULL

static int lapack_mmtk_CheckObject(PyObject *ob, int t, char *obname,
	char *tname, char *funname)
{	char buf[255];
	if (! PyArray_Check(ob))
	{	sprintf(buf,"Expected an array for parameter %s in lapack_mmtk.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->flags & CONTIGUOUS)) 
    {    sprintf(buf,"Parameter %s is not contiguous in lapack_mmtk.%s",obname, funname);
        ErrorReturn(buf);
        return 0;
    }
    if (!(((PyArrayObject *)ob)->descr->type_num == t))
    {    sprintf(buf,"Parameter %s is not of type %s in lapack_mmtk.%s",obname,tname, funname);
        ErrorReturn(buf);
        return 0;
    }
    return 1;
}

#define LDATA(p) ((long int *) (((PyArrayObject *)p)->data))
#define CHDATA(p) ((char *) (((PyArrayObject *)p)->data))
#define SHDATA(p) ((short int *) (((PyArrayObject *)p)->data))
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define FDATA(p) ((float *) (((PyArrayObject *)p)->data))
#define CDATA(p) ((f2c_complex *) (((PyArrayObject *)p)->data))
#define ZDATA(p) ((f2c_doublecomplex *) (((PyArrayObject *)p)->data))




static PyObject *lapack_mmtk_dsyev(PyObject *self, PyObject *args)
{
int  lapack_mmtk_status__;
char jobz;
char uplo;
long int n;
PyObject *a;
long int lda;
PyObject *w;
PyObject *work;
long int lwork;
long int info;
TRY(PyArg_ParseTuple(args,"cclOlOOll",&jobz,&uplo,&n,&a,&lda,&w,&work,&lwork,&info));

TRY(lapack_mmtk_CheckObject(a,PyArray_DOUBLE,"a","PyArray_DOUBLE","dsyev")); 
TRY(lapack_mmtk_CheckObject(w,PyArray_DOUBLE,"w","PyArray_DOUBLE","dsyev")); 
TRY(lapack_mmtk_CheckObject(work,PyArray_DOUBLE,"work","PyArray_DOUBLE","dsyev")); 

Py_BEGIN_ALLOW_THREADS;
#if defined(NO_APPEND_FORTRAN)
lapack_mmtk_status__ = dsyev(&jobz,&uplo,&n,DDATA(a),&lda,DDATA(w),DDATA(work),&lwork,&info);
#else
lapack_mmtk_status__ = dsyev_(&jobz,&uplo,&n,DDATA(a),&lda,DDATA(w),DDATA(work),&lwork,&info);
#endif
Py_END_ALLOW_THREADS;

return Py_BuildValue("{s:i,s:c,s:c,s:i,s:i,s:i,s:i}","dsyev_",lapack_mmtk_status__,"jobz",jobz,"uplo",uplo,"n",n,"lda",lda,"lwork",lwork,"info",info);
}

static PyObject *lapack_mmtk_dgesvd(PyObject *self, PyObject *args)
{
int  lapack_mmtk_status__;
char jobu;
char jobvt;
long int m;
long int n;
PyObject *a;
long int lda;
PyObject *s;
PyObject *u;
long int ldu;
PyObject *vt;
long int ldvt;
PyObject *work;
long int lwork;
long int info;
TRY(PyArg_ParseTuple(args,"ccllOlOOlOlOll",&jobu,&jobvt,&m,&n,&a,&lda,&s,&u,&ldu,&vt,&ldvt,&work,&lwork,&info));

TRY(lapack_mmtk_CheckObject(a,PyArray_DOUBLE,"a","PyArray_DOUBLE","dgesvd")); 
TRY(lapack_mmtk_CheckObject(s,PyArray_DOUBLE,"s","PyArray_DOUBLE","dgesvd")); 
TRY(lapack_mmtk_CheckObject(u,PyArray_DOUBLE,"u","PyArray_DOUBLE","dgesvd")); 
TRY(lapack_mmtk_CheckObject(vt,PyArray_DOUBLE,"vt","PyArray_DOUBLE","dgesvd")); 
TRY(lapack_mmtk_CheckObject(work,PyArray_DOUBLE,"work","PyArray_DOUBLE","dgesvd")); 


Py_BEGIN_ALLOW_THREADS;
#if defined(NO_APPEND_FORTRAN)
lapack_mmtk_status__ = dgesvd(&jobu,&jobvt,&m,&n,DDATA(a),&lda,DDATA(s),DDATA(u),&ldu,DDATA(vt),&ldvt,DDATA(work),&lwork,&info);
#else
lapack_mmtk_status__ = dgesvd_(&jobu,&jobvt,&m,&n,DDATA(a),&lda,DDATA(s),DDATA(u),&ldu,DDATA(vt),&ldvt,DDATA(work),&lwork,&info);
#endif
Py_END_ALLOW_THREADS;

return Py_BuildValue("{s:i,s:c,s:c,s:i,s:i,s:i,s:i,s:i,s:i,s:i}","dgesvd_",lapack_mmtk_status__,"jobu",jobu,"jobvt",jobvt,"m",m,"n",n,"lda",lda,"ldu",ldu,"ldvt",ldvt,"lwork",lwork,"info",info);
}


static struct PyMethodDef lapack_mmtk_module_methods[] = {
{"dsyev",lapack_mmtk_dsyev,1},
{"dgesvd",lapack_mmtk_dgesvd,1},
{ NULL,NULL,0}
};

static PyObject *lapack_mmtkError;
void initlapack_mmtk()
{ PyObject *m,*d;
m = Py_InitModule("lapack_mmtk",lapack_mmtk_module_methods);
#ifdef import_array
import_array();
#endif
d = PyModule_GetDict(m);
ErrorObject = LapackError();
PyDict_SetItemString(d,"LapackError",ErrorObject);
}
