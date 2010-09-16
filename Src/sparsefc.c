/* Sparse force-constant matrix objects.
 *
 * Written by Konrad Hinsen
 */

#define NO_IMPORT
#define _FORCEFIELD_MODULE
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_MMTKFF_API

#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"


/* Type "SparseForceConstants"
 *
 * The type declaration is in forcefield.h.
 */

/* Set array elements to zero */

void
PySparseFC_Zero(PySparseFCObject *fc)
{
  int i, j, k;

  for (i = 0; i < fc->nalloc; i++) {
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
	fc->data[i].fc[j][k] = 0.;
  }
}

/* Allocation and deallocation */

PySparseFCObject *
PySparseFC_New(int natoms, int nalloc)
{
  PySparseFCObject *self;
  int i;
  self = PyObject_NEW(PySparseFCObject, &PySparseFC_Type);
  if (self == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  if (nalloc < natoms)
    nalloc = natoms;
  self->data = (struct pair_fc *)malloc(nalloc*sizeof(struct pair_fc));
  self->index = (struct pair_descr_list *) \
                malloc(2*natoms*sizeof(struct pair_descr_list));
  if (self->data == NULL || self->index == NULL) {
    if (self->data != NULL)
      free(self->data);
    if (self->index != NULL)
      free(self->index);
    PyObject_Del(self);
    PyErr_NoMemory();
    return NULL;
  }
  self->natoms = natoms;
  self->nalloc = nalloc;
  self->nused = natoms;
  for (i = 0; i < 2*natoms; i++) {
    self->index[i].nalloc = self->index[i].nused = 0;
    self->index[i].list = NULL;
  }
  for (i = 0; i < natoms; i++) {
    self->data[i].i = i;
    self->data[i].j = i;
  }
  PySparseFC_Zero(self);
  self->cutoff_sq = 0.;
  self->fc_fn = sparse_fc_function;
  return self;
}

static void
sparsefc_dealloc(PySparseFCObject *self)
{
  int i;
  for (i = 0; i < 2*self->natoms; i++)
    if (self->index[i].nalloc > 0)
      free(self->index[i].list);
  free(self->index);
  free(self->data);
  PyObject_Del(self);
}

/* Find an entry */

struct pair_descr *
sparsefc_find(PySparseFCObject *fc, int i, int j)
{
  int sum = i+j;
  int diff = j-i;
  struct pair_descr_list *pairs = fc->index + sum;
  struct pair_descr *entry = pairs->list;
  int k;
  for (k = 0; k < pairs->nused; k++) {
    if (entry->diffij == diff) {
      if (fc->data[entry->index].i != i || fc->data[entry->index].j != j) {
	printf("Index error\n");
      }
      return entry;
    }
    entry++;
  }
  if (k < pairs->nalloc)
    return entry;
  else
    return NULL;
}

double *
PySparseFC_Find(PySparseFCObject *fc, int i, int j)
{
  struct pair_descr *pair;
  if (i == j)
    return (double *)fc->data[i].fc;
  pair = sparsefc_find(fc, i, j);
  if (pair == NULL)
    return NULL;
  if (pair->diffij < 0)
    return NULL;
  return (double *)fc->data[pair->index].fc;
}

/* Add a term */

int
PySparseFC_AddTerm(PySparseFCObject *fc, int i, int j, double *term)
{
  double *entry;
  int l;
  if (j < i)
    return 0;
  if (i == j)
    entry = (double *)fc->data[i].fc;
  else {
    struct pair_descr *pair = sparsefc_find(fc, i, j);
    if (pair == NULL) {
      struct pair_descr_list *list = fc->index + (i+j);
      int incr = (fc->nalloc/(2*fc->natoms));
      void *new;
      if (incr < 1)
	incr = 1;
#if 0
      printf("Extending list %d from %d to %d\n", i+j,
	     list->nalloc, list->nalloc+incr);
#endif
      new = realloc(list->list, (list->nalloc+incr)*sizeof(struct pair_descr));
      if (new == NULL)
	return 0;
      list->list = (struct pair_descr *)new;
      list->nalloc += incr;
      for (l = list->nused; l < list->nalloc; l++)
	list->list[l].diffij = -1;
      pair = list->list + list->nused;
    }
    if (pair->diffij < 0) {
      if (fc->nused == fc->nalloc) {
	int incr = fc->nalloc/10;
	void *new;
	if (incr < 1)
	  incr = 1;
#if 0
	printf("Extending data space from %d to %d\n",
	       fc->nalloc, fc->nalloc+incr);
#endif
	new = realloc(fc->data, (fc->nalloc+incr)*sizeof(struct pair_fc));
	if (new == NULL)
	  return 0;
	fc->data = (struct pair_fc *)new;
	fc->nalloc += incr;
	for (l = fc->nused; l < fc->nalloc; l++) {
	  double *data = (double *)fc->data[l].fc;
	  int k;
	  for (k = 0; k < 9; k++)
	    *data++ = 0.;
	}
      }
      pair->index = fc->nused;
      fc->nused++;
      pair->diffij = j-i;
      fc->index[i+j].nused++;
      fc->data[pair->index].i = i;
      fc->data[pair->index].j = j;
    }
#if 0
    if (fc->data[pair->index].i != i || fc->data[pair->index].j != j)
      printf("Incorrect pair entry: %d/%d, %d/%d\n",
	     fc->data[pair->index].i, i, fc->data[pair->index].j, j);
#endif
    entry = (double *)fc->data[pair->index].fc;
  }
  for (l = 0; l < 9; l++)
    *entry++ += *term++;
  return 1;
}

/* Extract a subarray */

void
PySparseFC_CopyToArray(PySparseFCObject *fc, double *data, int lastdim,
		       int from1, int to1, int from2, int to2)
{
  int i, j;
  double *temp;
  int pairs;

  temp = data;
  for (i = 0; i < 3*(to2-from2); i++) {
    for (j = 0; j < 3*(to1-from1); j++)
      temp[j] = 0.;
    temp += lastdim;
  }

  pairs = (to1-from1)*(to2-from2);
  for (i = 0; i < fc->nused; i++) {
    if (fc->data[i].i >= from1 && fc->data[i].i < to1 &&
	fc->data[i].j >= from2 && fc->data[i].j < to2) {
      int offset = 3*lastdim*(fc->data[i].i-from1)+3*(fc->data[i].j-from2);
      int l, m;
      for (l = 0; l < 3; l++) {
	for (m = 0; m < 3; m++)
	  data[offset+m] = fc->data[i].fc[l][m];
	offset += lastdim;
      }
      pairs--;
    }
    if (fc->data[i].i != fc->data[i].j &&
	fc->data[i].j >= from1 && fc->data[i].j < to1 &&
	fc->data[i].i >= from2 && fc->data[i].i < to2) {
      int offset = 3*lastdim*(fc->data[i].j-from1)+3*(fc->data[i].i-from2);
      int l, m;
      for (l = 0; l < 3; l++) {
	for (m = 0; m < 3; m++)
	  data[offset+m] = fc->data[i].fc[m][l];
	offset += lastdim;
      }
      pairs--;
    }
    if (pairs == 0)
      break;
  }
}

PyObject *
PySparseFC_AsArray(PySparseFCObject *fc, int from1, int to1, int from2, int to2)
{
  PyArrayObject *array;
#if defined(NUMPY)
  npy_intp dims[4];
#else
  int dims[4];
#endif
  dims[0] = to1-from1;
  if (dims[0] < 0) dims[0] = 0;
  dims[2] = to2-from2;
  if (dims[2] < 0) dims[2] = 0;
  dims[1] = dims[3] = 3;
#if defined(NUMPY)
  array = (PyArrayObject *)PyArray_SimpleNew(4, dims, PyArray_DOUBLE);
#else
  array = (PyArrayObject *)PyArray_FromDims(4, dims, PyArray_DOUBLE);
#endif
  if (array == NULL)
    return NULL;
  PySparseFC_CopyToArray(fc, (double *)array->data, 3*dims[2],
			 from1, to1, from2, to2);
  return (PyObject *)array;
}

/* Multiply with a vector */

void
PySparseFC_VectorMultiply(PySparseFCObject *fc, double *result, double *vector,
			 int from_i, int to_i, int from_j, int to_j)
{
  struct pair_fc *entry = fc->data;
  int i;

  for (i = 3*(to_i-from_i)-1; i >= 0; i--)
    result[i] = 0.;
  for (i = 0; i < fc->nused; i++) {
    int l, m;
    if (entry->i >= from_i && entry->i < to_i
	&& entry->j >= from_j && entry->j < to_j) {
      for (l = 0; l < 3; l++)
	for (m = 0; m < 3; m++)
	  result[3*(entry->i-from_i)+l] +=
	    entry->fc[l][m]*vector[3*(entry->j-from_j)+m];
    }
    if (entry->i != entry->j &&
	entry->j >= from_i && entry->j < to_i
	&& entry->i >= from_j && entry->i < to_j) {
      for (l = 0; l < 3; l++)
	for (m = 0; m < 3; m++)
	  result[3*(entry->j-from_i)+m] +=
	    entry->fc[l][m]*vector[3*(entry->i-from_j)+l];
    }
    entry++;
  }  
}

/* Multiply with a vector */

static void
solve_3x3(double *fc, double *vector, double *result)
{
  double a = fc[0*3+0];
  double b = fc[1*3+1];
  double c = fc[2*3+2];
  double d = fc[0*3+1];
  double e = fc[0*3+2];
  double f = fc[1*3+2];
  double o = vector[0];
  double p = vector[1];
  double q = vector[2];

  double af_de = a*f-d*e;
  double aq_eo = a*q-e*o;
  double ab_dd = a*b-d*d;
  double ac_ee = a*c-e*e;

  double z = (af_de*(a*p-d*o)-ab_dd*aq_eo) / (af_de*af_de-ab_dd*ac_ee);
  double y = (aq_eo - z*ac_ee)/af_de;
  double x = (o - d*y - e*z)/a;

  result[0] = x;
  result[1] = y;
  result[2] = z;
}

int
PySparseFC_VectorSolve(PySparseFCObject *fc, double *result, double *vector,
		       double tolerance, int max_iter)
{
  int natoms = fc->natoms;
  int nc = 3*natoms;
  double *work, *r, *z, *p, *q;
  double rho, rho_prev, alpha, beta, residual;
  int i, iter;

  work = (double *)malloc(4*nc*sizeof(double));
  if (work == NULL)
    return -1;
  r = work; z = r + nc; p = z + nc; q = p + nc;
  for (i = 0; i < nc; i++) {
    r[i] = vector[i];
    result[i] = 0.;
  }
  iter = 0;
  while (1) {
    for (i = 0; i < natoms; i++) {
#if 1
      /* Use 3x3 diagonal as preconditioner */
      double *diag = PySparseFC_Find(fc, i, i);
      solve_3x3(diag, &r[3*i], &z[3*i]);
#else
      /* No preconditioning */
      z[3*i] = r[3*i];
      z[3*i+1] = r[3*i+1];
      z[3*i+2] = r[3*i+2];      
#endif
    }
    rho_prev = rho;
    rho = 0.;
    for (i = 0; i < nc; i++)
      rho += r[i]*z[i];
    if (iter == 0) {
      for (i = 0; i < nc; i++)
	p[i] = z[i];
    }
    else {
      beta = rho/rho_prev;
      for (i = 0; i < nc; i++)
	p[i] = z[i] + beta*p[i];
    }
    PySparseFC_VectorMultiply(fc, q, p, 0, natoms, 0, natoms);
    alpha = 0.;
    for (i = 0; i < nc; i++)
      alpha += p[i]*q[i];
    alpha = rho/alpha;
    residual = 0.;
    for (i = 0; i < nc; i++) {
      result[i] += alpha*p[i];
      r[i] -= alpha*q[i];
      residual += r[i]*r[i];
    }
    residual = sqrt(residual/natoms);
#ifdef DEBUG
    printf("Residual: %lg\n", residual);
#endif
    iter ++;
    if (iter > 2 && residual < tolerance)
      break;
    if (iter > max_iter) {
      free(work);
      return 0;
    }
  }
#ifdef DEBUG
  printf("CG: %d iterations\n", iter);
  fflush(stdout);
#endif
  free(work);
  return 1;
}


/* Scale by weight vector */

void
PySparseFC_Scale(PySparseFCObject *fc, PyArrayObject *factors)
{
  struct pair_fc *entry = fc->data;
  double *f = (double *)factors->data;
  int i;

  for (i = 0; i < fc->nused; i++) {
    int l, m;
    for (l = 0; l < 3; l++)
      for (m = 0; m < 3; m++)
	entry->fc[l][m] *= f[entry->i]*f[entry->j];
    entry++;
  }
}

/* Methods */

static PyObject *
setCutoff(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  double cutoff;
  if (!PyArg_ParseTuple(args, "d", &cutoff))
    return NULL;
  fc->cutoff_sq = cutoff*cutoff;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
accessFunction(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  return PyCObject_FromVoidPtr((void *)fc->fc_fn, NULL);
}

static PyObject *
multiplyVector(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  PyArrayObject *vector = NULL;
  PyObject *result = NULL;
  PyArrayObject *result_array;
  int from_i = 0, to_i = fc->natoms;
  int from_j = 0, to_j = fc->natoms;
#if defined(NUMPY)
  npy_intp dims[2];
#else
  int dims[2];
#endif
  if (!PyArg_ParseTuple(args, "O!|Oiiii", &PyArray_Type, &vector, &result,
			&from_i, &to_i, &from_j, &to_j))
    return NULL;
  if (result == Py_None)
    result = NULL;
  if (result != NULL) {
    if (!PyArray_Check(result)) {
      PyErr_SetString(PyExc_TypeError, "result must be array");
      return NULL;
    }
    result_array = (PyArrayObject *)result;
    if (result_array->nd != 2 || result_array->dimensions[0] != to_i-from_i
	|| result_array->dimensions[1] != 3) {
      PyErr_SetString(PyExc_ValueError, "illegal array shape");
      return NULL;
    }
  }
  if (vector->nd != 2 || vector->dimensions[0] != to_j-from_j
      || vector->dimensions[1] != 3) {
    PyErr_SetString(PyExc_ValueError, "illegal array shape");
    return NULL;
  }
  if (from_i < 0 || to_i > fc->natoms || to_i < from_i
      || from_j < 0 || to_j > fc->natoms || to_j < from_j) {
    PyErr_SetString(PyExc_ValueError, "illegal subset");
    return NULL;
  }
  if (result == NULL) {
    dims[0] = to_i-from_i;
    dims[1] = 3;
#if defined(NUMPY)
    result = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
#else
    result = PyArray_FromDims(2, dims, PyArray_DOUBLE);
#endif
    if (result == NULL)
      return NULL;
  }
  else
    Py_INCREF(result);
  result_array = (PyArrayObject *)result;
  PySparseFC_VectorMultiply(fc, (double *)result_array->data,
			    (double *)vector->data,
			    from_i, to_i, from_j, to_j);
  return (PyObject *)result;
}

static PyObject *
solveForVector(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  PyArrayObject *vector = NULL;
  PyObject *result = NULL;
  PyArrayObject *result_array;
  double tolerance = 1.e-8;
#if defined(NUMPY)
  npy_intp dims[2];
#else
  int dims[2];
#endif
  int max_iter = 0;
  int ret;
  if (!PyArg_ParseTuple(args, "O!|Odi", &PyArray_Type, &vector, &result,
			                &tolerance, &max_iter))
    return NULL;
  if (result == Py_None)
    result = NULL;
  if (result != NULL) {
    if (!PyArray_Check(result)) {
      PyErr_SetString(PyExc_TypeError, "result must be array");
      return NULL;
    }
    result_array = (PyArrayObject *)result;
    if (result_array->nd != 2 || result_array->dimensions[0] != fc->natoms
	|| result_array->dimensions[1] != 3) {
      PyErr_SetString(PyExc_ValueError, "illegal array shape");
      return NULL;
    }
  }
  if (vector->nd != 2 || vector->dimensions[0] != fc->natoms
      || vector->dimensions[1] != 3) {
    PyErr_SetString(PyExc_ValueError, "illegal array shape");
    return NULL;
  }
  if (result == NULL) {
    dims[0] = fc->natoms;
    dims[1] = 3;
#if defined(NUMPY)
    result = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
#else
    result = PyArray_FromDims(2, dims, PyArray_DOUBLE);
#endif
    if (result == NULL)
      return NULL;
  }
  else
    Py_INCREF(result);
  result_array = (PyArrayObject *)result;
  if (max_iter == 0)
    max_iter = 4*fc->natoms;
  ret = PySparseFC_VectorSolve(fc, (double *)result_array->data,
			       (double *)vector->data,
			       tolerance, max_iter);
  if (ret == -1) {
    PyErr_NoMemory();
    Py_DECREF(result);
    return NULL;
  }
  if (ret == 0) {
    PyErr_SetString(PyExc_ValueError, "no convergence");
    Py_DECREF(result);
    return NULL;
  }
  return (PyObject *)result;
}

static PyObject *
asArray(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  return PySparseFC_AsArray(fc, 0, fc->natoms, 0, fc->natoms);
}

static PyObject *
scale(PyObject *self, PyObject *args)
{
  PySparseFCObject *fc = (PySparseFCObject *)self;
  PyArrayObject *factors;
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &factors))
    return NULL;
  PySparseFC_Scale(fc, factors);
  Py_INCREF(Py_None);
  return Py_None;
}

/* Documentation string */

static char PySparseFC_Type__doc__[] = 
  "sparse force constant matrix";

/* Method table */

static struct PyMethodDef sparsefc_methods[] = {
  {"setCutoff", setCutoff, 1},
  {"accessFunction", accessFunction, 1},
  {"multiplyVector", multiplyVector, 1},
  {"solveForVector", solveForVector, 1},
  {"asArray", asArray, 1},
  {"scale", scale, 1},
  {NULL, NULL} /* sentinel */
};

/* Attribute access. */

static PyObject *
sparsefc_getattr(PySparseFCObject *self, char *name)
{
  return Py_FindMethod(sparsefc_methods, (PyObject *)self, name);
}

/* Sequence protocol */

static Py_ssize_t
sparsefc_length(PySparseFCObject *self)
{
  return (Py_ssize_t)self->nused;
}

static PyObject *
sparsefc_item(PySparseFCObject *self, Py_ssize_t i)
{
  if (i < 0 || i >= self->nused) {
    PyErr_SetString(PyExc_IndexError,"index out of bounds");
    return NULL;
  }
  else {
    PyArrayObject *array;
    PyObject *ret;
#if defined(NUMPY)
    npy_intp dims[2];
#else
    int dims[2];
#endif
    dims[0] = 3; dims[1] = 3;
#if defined(NUMPY)
    array = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
#else
    array = (PyArrayObject *)PyArray_FromDims(2, dims, PyArray_DOUBLE);
#endif
    if (array == NULL)
      return NULL;
    memcpy(array->data, self->data[i].fc, 9*sizeof(double));
    ret = Py_BuildValue("iiO", self->data[i].i, self->data[i].j,
			(PyObject *)array);
    Py_DECREF(array);
    return ret;
  }
}

/* Mapping protocol */

static PyObject *
sparsefc_subscript(PySparseFCObject *self, PyObject *index)
{
  PyArrayObject *array;
  Py_ssize_t from[2], to[2], rank[2], stride;
  int i;

  if (PyInt_Check(index))
    return sparsefc_item(self, PyInt_AsLong(index));
  if (!PyTuple_Check(index) || PyTuple_Size(index) != 2) {
    PyErr_SetString(PyExc_TypeError, "index must be tuple of length 2");
    return NULL;
  }
  for (i = 0; i < 2; i++) {
    PyObject *subscript = PyTuple_GetItem(index, (Py_ssize_t)i);
    if (PyInt_Check(subscript)) {
      from[i] = PyInt_AsLong(subscript);
      to[i] = from[i] + 1;
      rank[i] = 0;
      stride = 1;
    }
    else if (PySlice_Check(subscript)) {
      PySlice_GetIndices((PySliceObject *)subscript, self->natoms,
			 &from[i], &to[i], &stride);
      rank[i] = 1;
    }
    else {
      PyErr_SetString(PyExc_TypeError, "illegal subscript type");
      return NULL;
    }
    if (from[i] < 0 || to[i] > self->natoms || to[i] < from[i] || stride != 1) {
      PyErr_SetString(PyExc_IndexError, "illegal subscript");
      return NULL;
    }
  }
  if (rank[0] != rank[1]) {
    PyErr_SetString(PyExc_IndexError, "illegal subscript");
    return NULL;
  }
  array = (PyArrayObject *)PySparseFC_AsArray(self, from[0], to[0],
					      from[1], to[1]);
  if (array == NULL)
    return NULL;
  if (rank[0] == 0) {
    PyObject *shape = PyTuple_New((Py_ssize_t)2);
    if (shape == NULL)
      return NULL;
    PyTuple_SetItem(shape, (Py_ssize_t)0, PyInt_FromLong(3));
    PyTuple_SetItem(shape, (Py_ssize_t)1, PyInt_FromLong(3));
    array = (PyArrayObject *)PyArray_Reshape(array, shape);
    Py_DECREF(shape);
  }
  return (PyObject *)array;
}

/* Type object */

static PySequenceMethods sparsefc_as_sequence = {
  (lenfunc)sparsefc_length,   /*sq_length*/
  0,                          /*nb_add, concat is numeric add*/
  0,                          /*nb_multiply, repeat is numeric multiply*/
  (ssizeargfunc)sparsefc_item,/* sq_item*/
  0,                          /*sq_slice*/
  0,                          /*sq_ass_item*/
  0,                          /*sq_ass_slice*/
};

static PyMappingMethods sparsefc_as_mapping = {
  (lenfunc)sparsefc_length,	    /*mp_length*/
  (binaryfunc)sparsefc_subscript,   /*mp_subscript*/
  0,                                /*mp_ass_subscript*/
};

PyTypeObject PySparseFC_Type = {
  PyObject_HEAD_INIT(NULL)
  0,			          /*ob_size*/
  "SparseForceConstants",	  /*tp_name*/
  sizeof(PySparseFCObject),       /*tp_basicsize*/
  0,			          /*tp_itemsize*/
  /* methods */
  (destructor)sparsefc_dealloc,   /*tp_dealloc*/
  0,			          /*tp_print*/
  (getattrfunc)sparsefc_getattr,  /*tp_getattr*/
  0, 			          /*tp_setattr*/
  0,			          /*tp_compare*/
  0,                              /*tp_repr*/
  0,                              /*tp_as_number*/
  &sparsefc_as_sequence,	  /*tp_as_sequence*/
  &sparsefc_as_mapping,	          /*tp_as_mapping*/
  0,			          /*tp_hash*/
  0,	                          /*tp_call*/
  0,                              /*tp_str*/
  0,                              /*tp_getattro*/
  0,                              /*tp_setattro*/
  /* Space for future expansion */
  0L,0L,
  /* Documentation string */
  PySparseFC_Type__doc__
};

/* Interface to force field evaluators */

int
sparse_fc_function(energy_data *energy,
		   int i, int j, tensor3 term, double r_sq)
{
  PySparseFCObject *fc = (PySparseFCObject *)energy->force_constants;
  if (i < 0) {
    PySparseFC_Zero(fc);
    return 1;
  }
  else if (term == NULL) {
    return r_sq < fc->cutoff_sq || fc->cutoff_sq == 0.;
  }
  else if (r_sq < fc->cutoff_sq || fc->cutoff_sq == 0.) {
    if (!PySparseFC_AddTerm(fc, i, j, (double *)term)) {
      energy->error = 1;
      PyErr_SetString(PyExc_IndexError, "couldn't access sparse array");
    }
    return 1;
  }
  else
    return 0;
}

/* Constructor function */

PyObject *
SparseForceConstants(PyObject *dummy, PyObject *args)
{
  int natoms, nalloc;

  if (!PyArg_ParseTuple(args, "ii", &natoms, &nalloc)) {
    return NULL;
  }
  return (PyObject *)PySparseFC_New(natoms, nalloc);
}
