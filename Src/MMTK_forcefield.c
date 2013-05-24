/* Low-level force field calculations
 *
 * Written by Konrad Hinsen
 */

#define _FORCEFIELD_MODULE
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_MMTKFF_API

#ifdef WITH_MPI
#define PyMPI_API PyMPI_API_MMTKFF
#endif
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

#define DEBUG 0
#define THREAD_DEBUG 0
#define MPI_DEBUG 0

/* Windows header that defines Sleep() */

#ifdef WITH_THREAD
#ifdef MS_WINDOWS
#include "Windows.h"
#endif
#endif

/* Global variables */

double electrostatic_energy_factor;
distance_fn *distance_vector_pointer;
distance_fn *orthorhombic_distance_vector_pointer;
distance_fn *parallelepipedic_distance_vector_pointer;

#ifdef WITH_MPI
static PyObject *PyExc_MPIError;
#endif

staticforward ff_eval_function evaluator;


/*
 * Threading support
 */

#ifdef WITH_THREAD
#ifdef _POSIX_THREADS
static int
allocate_barrier(barrierinfo *binfo) {
  if (pthread_mutex_init(&binfo->lock, NULL) != 0)
    return 0;
  if (pthread_cond_init(&binfo->cond, NULL) != 0)
    return 0;
  binfo->n = 0;
  return 1;
}

static void
deallocate_barrier(barrierinfo *binfo)
{
  pthread_mutex_destroy(&binfo->lock);
  pthread_cond_destroy(&binfo->cond);
  free(binfo);
}

void
barrier(barrierinfo *binfo, int thread_id, int nthreads)
{
  if (nthreads > 1) {
#if THREAD_DEBUG
    printf("Thread %d arrived at barrier\n", thread_id);
#endif
    pthread_mutex_lock(&binfo->lock);
    if (binfo->n == nthreads)
      binfo->n = 1;
    else
      binfo->n++;
    pthread_cond_broadcast(&binfo->cond);
#if THREAD_DEBUG
    printf("Thread %d: barrier count is %d\n", thread_id, binfo->n);
#endif
    while (binfo->n != nthreads) {
      pthread_cond_wait(&binfo->cond, &binfo->lock);
    }
    pthread_mutex_unlock(&binfo->lock);
#if THREAD_DEBUG
    printf("Thread %d continuing after barrier\n", thread_id);
#endif
  }
}
#else
static int
allocate_barrier(barrierinfo *binfo)
{
  binfo->lock = PyThread_allocate_lock();
  binfo->n = 0;
  return (binfo->lock != NULL);
}

static void
deallocate_barrier(barrierinfo *binfo)
{
  PyThread_free_lock(binfo->lock);
  free(binfo);
}

void
barrier(barrierinfo *binfo, int thread_id, int nthreads)
{
  int done = 0;
  if (nthreads > 1) {
    PyThread_acquire_lock(binfo->lock, 1);
    if (binfo->n == nthreads)
      binfo->n = 1;
    else
      binfo->n++;
    PyThread_release_lock(binfo->lock);
    while (!done) {
      PyThread_acquire_lock(binfo->lock, 1);
      done = (binfo->n == nthreads);
      PyThread_release_lock(binfo->lock);
    }
  }
}
#endif
#endif

/* String copy with memory allocation */

char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* Type "energy term"
 *
 * Objects of this type represent the individual force field terms.
 * They are called by an energy evaluator object (see below).
 *
 * The structure declaration is in mmtk_forcefield.h.
 */

/* Allocation and deallocation */

#ifdef EXTENDED_TYPES
static PyObject *
PyFFEnergyTerm_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
  PyFFEnergyTermObject *self;
  int i;
  self = (PyFFEnergyTermObject *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->user_info = NULL;
    self->universe_spec = NULL;
    self->scratch = NULL;
    for (i = 0; i < MMTK_MAX_DATA; i++)
      self->data[i] = NULL;
    self->evaluator_name = NULL;
    for (i = 0; i < MMTK_MAX_TERMS; i++)
      self->term_names[i] = NULL;
    self->threaded = 0;
    self->thread_safe = 0;
    self->parallelized = 0;
    self->n = self->nterms = 0;
  }
  return (PyObject *)self;
}

static int
PyFFEnergyTerm_init(PyObject *self_po, PyObject *args, PyObject *kw)
{
  PyFFEnergyTermObject *self = (PyFFEnergyTermObject *)self_po;
  char *name;
  PyObject *universe;
  PyObject *term_names;
  int thread_safe = 0;
  int i;
  if (! PyArg_ParseTuple(args, "OsO!|i",
			 &universe, &name,
			 &PyTuple_Type, &term_names,
			 &thread_safe))
        return -1; 
  self->nterms = (int)PyTuple_Size(term_names);
  if (self->nterms == 0) {
    PyErr_SetString(PyExc_ValueError, "at least one term name required");
    return -1;
  }
  if (self->nterms > MMTK_MAX_TERMS) {
    PyErr_SetString(PyExc_ValueError, "too many terms");
    return -1;
  }
  self->universe_spec = (PyUniverseSpecObject *)
                        PyObject_GetAttrString(universe, "_spec");
  if (self->universe_spec == NULL)
    return -1;
  Py_INCREF(self->universe_spec);
  self->evaluator_name = allocstring(name);
  if (self->evaluator_name == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  for (i = 0; i < self->nterms; i++) {
    PyObject *name = PyTuple_GetItem(term_names, (Py_ssize_t)i);
    if (!PyString_Check(name)) {
      PyErr_SetString(PyExc_TypeError, "term names must be strings");
      return -1;
    }
    self->term_names[i] = allocstring(PyString_AsString(name));
    if (self->term_names[i] == NULL) {
      PyErr_NoMemory();
      return -1;
    }
  }
  self->thread_safe = thread_safe;
  return 0;
}
#endif

PyFFEnergyTermObject *
PyFFEnergyTerm_New(void)
{
  int i;
  PyFFEnergyTermObject *self;
  self = PyObject_NEW(PyFFEnergyTermObject, &PyFFEnergyTerm_Type);
  if (self == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  self->user_info = NULL;
  self->universe_spec = NULL;
  self->scratch = NULL;
  for (i = 0; i < MMTK_MAX_DATA; i++)
    self->data[i] = NULL;
  self->evaluator_name = NULL;
  for (i = 0; i < MMTK_MAX_TERMS; i++)
    self->term_names[i] = NULL;
  self->threaded = 0;
  self->parallelized = 0;
  self->n = self->nterms = 0;
  return self;
}

static void
energyterm_dealloc(PyFFEnergyTermObject *self)
{
  int i;
  for (i = 0; i < self->nterms; i++)
    free(self->term_names[i]);
  Py_XDECREF(self->user_info);
  Py_XDECREF(self->universe_spec);
  for (i = 0; i < MMTK_MAX_DATA; i++)
    Py_XDECREF(self->data[i]);
  if (self->scratch != NULL)
    free(self->scratch);
  self->ob_type->tp_free((PyObject *)self);
}

/* Add nonbonded term to NonbondedListTerm */

static PyObject *
add_term(PyObject *self, PyObject *args)
{
  PyFFEnergyTermObject *ev = (PyFFEnergyTermObject *)self;
  PyFFEnergyTermObject *ev_term;
  int ev_type;
  if (!PyArg_ParseTuple(args, "O!i", &PyFFEnergyTerm_Type, &ev_term, &ev_type))
    return NULL;
  if (strcmp(ev->evaluator_name, "nonbonded list summation") != 0) {
    PyErr_SetString(PyExc_ValueError, "not a NonbondedListTerm");
    return NULL;
  }
  Py_INCREF(ev_term);
  ev->data[ev_type+1] = (PyObject *)ev_term;
  Py_INCREF(Py_None);
  return Py_None;
}

/* Documentation string */

static char PyFFEnergyTerm_Type__doc__[] = 
  "energy term in force field";

/* Method table */

static struct PyMethodDef energyterm_methods[] = {
  {"addTerm", add_term, 1},
  {NULL, NULL} /* sentinel */
};

/* Attribute and element access. */

static PyObject *
energyterm_getattr(PyFFEnergyTermObject *self, char *name)
{
  if (strcmp(name, "name") == 0) {
    return PyString_FromString(self->evaluator_name);
  }
  else if (strcmp(name, "term_names") == 0) {
    PyObject *ret = PyTuple_New((Py_ssize_t)self->nterms);
    int i;
    for (i = 0; i < self->nterms; i++)
      PyTuple_SetItem(ret, (Py_ssize_t)i,
		      PyString_FromString(self->term_names[i]));
    return ret;
  }
  else if (strcmp(name, "info") == 0) {
    if (self->user_info == NULL) {
      PyErr_SetString(PyExc_AttributeError, "attribute not defined");
      return NULL;
    }
    Py_INCREF(self->user_info);
    return self->user_info;
  }
  else
    return Py_FindMethod(energyterm_methods, (PyObject *)self, name);
}

static int
energyterm_setattr(PyFFEnergyTermObject *self, char *name, PyObject *value)
{
  if (strcmp(name, "info") == 0) {
    Py_XDECREF(self->user_info);
    Py_INCREF(value);
    self->user_info = value;
    return 0;
  }
  else {
    PyErr_SetString(PyExc_AttributeError, "attribute not defined");
    return -1;
  }
}

/* The garbage collector routines do nothing, but they need to be
   present because Pyrex derived types call them. */
static int
energyterm_traverse(PyFFEnergyTermObject *self, visitproc visit, void *arg)
{
  return 0;
}

static Py_ssize_t
energyterm_clear(PyFFEnergyTermObject *self)
{
  return 0;
}

/* Type object */

PyTypeObject PyFFEnergyTerm_Type = {
  PyObject_HEAD_INIT(NULL)
  0,			          /*ob_size*/
  "FFEnergyTerm",	          /*tp_name*/
  sizeof(PyFFEnergyTermObject),   /*tp_basicsize*/
  0,			          /*tp_itemsize*/
  /* methods */
  (destructor)energyterm_dealloc, /*tp_dealloc*/
  0,			          /*tp_print*/
  (getattrfunc)energyterm_getattr,/*tp_getattr*/
  (setattrfunc)energyterm_setattr,/*tp_setattr*/
  0,			          /*tp_compare*/
  0,                              /*tp_repr*/
  0,                              /*tp_as_number*/
  0,                              /*tp_as_sequence*/
  0,			          /*tp_as_mapping*/
  0,			          /*tp_hash*/
  0,                              /*tp_call*/
  0,                              /*tp_str*/
  0,                              /*tp_getattro*/
  0,                              /*tp_setattro*/
  /* Functions to access object as input/output buffer */
  /* PyBufferProcs *tp_as_buffer */
  0,
  /* long tp_flags */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  /* Documentation string */
  PyFFEnergyTerm_Type__doc__

#ifdef EXTENDED_TYPES
  ,
  /* Assigned meaning in release 2.0 */

  /* call function for all accessible objects */
  /* traverseproc tp_traverse */
  (traverseproc)energyterm_traverse,

  /* delete references to contained objects */
  /* inquiry tp_clear */
  (inquiry)energyterm_clear,

  /* Assigned meaning in release 2.1 */

  /* rich comparisons */
  /* richcmpfunc tp_richcompare */
  0,

  /* weak reference enabler */
  /* long tp_weaklistoffset */
  0,

  /* Added in release 2.2 */

  /* Iterators */
  /* getiterfunc tp_iter */
  0,
  /* iternextfunc tp_iternext */
  0,

  /* Attribute descriptor and subclassing stuff */
  /* struct PyMethodDef *tp_methods */
  0,
  /* struct PyMemberDef *tp_members */
  0,
  /* struct PyGetSetDef *tp_getset */
  0,
  /* struct _typeobject *tp_base */
  0,
  /* PyObject *tp_dict */
  0,
  /* descrgetfunc tp_descr_get */
  0,
  /* descrsetfunc tp_descr_set */
  0,
  /* long tp_dictoffset */
  0,
  /* initproc tp_init */
  PyFFEnergyTerm_init,
  /* allocfunc tp_alloc */
  0,
  /* newfunc tp_new */
  PyFFEnergyTerm_new,
  /* freefunc tp_free */
  0,
  /* inquiry tp_is_gc */
  0,
  /* PyObject *tp_bases */
  0,
  /* PyObject *tp_mro */
  0,
  /* PyObject *tp_cache */
  0,
  /* PyObject *tp_subclasses */
  0,
  /* PyObject *tp_weaklist */
  0,
  /* destructor tp_del */
  0
#endif
};



/* Type "evaluator"
 *
 * Objects of this type are callable, expecting a coordinate array,
 * a gradient array (None if no gradients are to be calculated),
 * and a force constant array (None if no force constants are to be
 * calculated). They return a float object containing the energy.
 * For efficient access from C code, these objects have an attribute
 * cfunc, which is a pointer to an equivalent C function.
 *
 * The structure declaration is in mmtk_forcefield.h.
 */

/* Allocation and deallocation */

PyFFEvaluatorObject *
PyFFEvaluator_New(void)
{
  PyFFEvaluatorObject *self;
  self = PyObject_NEW(PyFFEvaluatorObject, &PyFFEvaluator_Type);
  if (self == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  self->universe_spec = NULL;
  self->terms = NULL;
  self->energy_terms = NULL;
  self->nterms =  self->ntermobjects = 0;
  self->scratch = NULL;
  self->nthreads = 0;
#ifdef WITH_THREAD
  self->global_lock = NULL;
  self->binfo = NULL;
#endif
  return self;
}

static void
evaluator_dealloc(PyFFEvaluatorObject *self)
{
  int i;
#ifdef WITH_THREAD
  if (self->eval_func == evaluator) {
    threadinfo *tinfo = (threadinfo *)self->scratch;
    if (self->global_lock != NULL)
      PyThread_free_lock(self->global_lock);
    if (self->binfo != NULL)
      deallocate_barrier(self->binfo);
    for (i = 1; i < self->nthreads; i++) {
      int j = 50;
      tinfo->exit = 1;
#if THREAD_DEBUG
      printf("Releasing thread %d\n", tinfo->input.thread_id);
#endif
      PyThread_release_lock(tinfo->lock);
      while (!tinfo->stop && j--) {
#if THREAD_DEBUG
        printf("evaluator_dealloc: Waiting for thread %d to stop (tinfo->stop is %d)\n",
               tinfo->input.thread_id, tinfo->stop);
#endif
#ifdef MS_WINDOWS
	Sleep(10);  /* 10 ms */
#else
        {
          struct timeval tv;
          tv.tv_sec = 0;
          tv.tv_usec = 10000;  /* 10 ms */
          select(0, NULL, NULL, NULL, &tv);
        }
#endif
      }
      Py_XDECREF(tinfo->energy.gradients);
      free(tinfo->energy.energy_terms);
      PyThread_free_lock(tinfo->lock);
      tinfo++;
    }
  }
#endif
#ifdef WITH_MPI
  if (self->energy_parts)
    free(self->energy_parts);
  if (self->gradient_parts)
    free(self->gradient_parts);
#endif
  Py_XDECREF(self->universe_spec);
  Py_XDECREF(self->terms);
  Py_XDECREF(self->energy_terms_array);
  if (self->scratch != NULL)
    free(self->scratch);
  PyObject_Del(self);
}

/* Call */

#if 0
int
test_gradient_function(struct ffeval *evaluator,
		       PyObject *data, int i, vector3 gradient)
{
  PyArrayObject *array = (PyArrayObject *)data;
  vector3 *f = (vector3 *)array->data;
  int natoms = array->dimensions[0];
  if (i < 0) {
    int i;
#if 0
    printf("Zeroing gradient array\n");
#endif
    for (i = 0; i < natoms; i++) {
      f[i][0] = 0.;
      f[i][1] = 0.;
      f[i][2] = 0.;
    }
  }
  else {
#if 0
    printf("Gradient for atom %d from %s\n", i, evaluator->name);
    printf("  %lf, %lf, %lf\n", gradient[0], gradient[1], gradient[2]);
#endif
    f[i][0] += gradient[0];
    f[i][1] += gradient[1];
    f[i][2] += gradient[2];
  }
  return 1;
}

int
test_fc_function(struct ffeval *evaluator, PyObject *data,
		 int i, int j, tensor3 fc, double r_sq)
{
  PyArrayObject *array = (PyArrayObject *)data;
  double *sd = (double *)array->data;
  int natoms = array->dimensions[0];
  if (i < 0) {
    long i;
#if 0
    printf("Zeroing force constant array\n");
#endif
    for (i = 0; i < 9L*(long)natoms*(long)natoms; i++)
      sd[i] = 0.;
  }
  else if (fc != NULL) {
    double *fcij = sd + 9L*(long)natoms*(long)i + 3L*(long)j;
    int l, m;
#if 0
    printf("Force constants for atom pair (%d, %d) with distance %lf from %s\n",
	   i, j, sqrt(r_sq), evaluator->name);
#endif
    if (i > j)
      printf("Illegal index pair: %d > %d\n", i, j);
    for (l = 0; l < 3; l++) {
      for (m = 0; m < 3; m++)
	fcij[m] += fc[l][m];
      fcij += 3*natoms;
    }
  }
  return 1;
}
#endif

static PyObject *
evaluator_call(PyFFEvaluatorObject *self, PyObject *args)
{
  PyArrayObject *coordinates = NULL;
  PyObject *gradients = NULL;
  PyObject *force_constants = NULL;
  int small_change = 0;
  gradient_function *gf = NULL;
  fc_function *fcf = NULL;
  energy_data energy;
  if (!PyArg_ParseTuple(args, "O!|OOi",
			&PyArray_Type, &coordinates,
			&gradients, &force_constants,
			&small_change))
    return NULL;
  if (gradients == Py_None)
    gradients = NULL;
  if (force_constants == Py_None)
    force_constants = NULL;
  if (gradients != NULL && !PyArray_Check(gradients)) {
    PyObject *fnptr = PyObject_CallMethod(gradients, "accessFunction", NULL);
    if (fnptr == NULL)
      return NULL;
    gf = (gradient_function *)PyCObject_AsVoidPtr(fnptr);
  }
  if (force_constants != NULL && !PyArray_Check(force_constants)) {
    PyObject *fnptr = PyObject_CallMethod(force_constants,
					  "accessFunction", NULL);
    if (fnptr == NULL)
      return NULL;
    fcf = (fc_function *)PyCObject_AsVoidPtr(fnptr);
  }
  energy.gradients = gradients;
  energy.gradient_fn = gf;
  energy.force_constants = force_constants;
  energy.fc_fn = fcf;
#ifdef WITH_THREAD
  self->tstate_save = PyEval_SaveThread();
#endif
  (*self->eval_func)(self, &energy, coordinates, small_change);
#ifdef WITH_THREAD
  PyEval_RestoreThread(self->tstate_save);
#endif
  if (energy.error)
    return NULL;
  else
    return PyFloat_FromDouble(energy.energy);
}

static void
evaluator(PyFFEvaluatorObject *self,
	  energy_data *energy,
	  PyArrayObject *coordinates, int small_change)
{
  int natoms = coordinates->dimensions[0];
  PyFFEnergyTermObject *term;
  energy_spec input;
  int i;

  input.coordinates = coordinates;
  input.natoms = natoms;
  input.small_change = small_change;
  input.nthreads = self->nthreads;
  input.nprocs = self->nprocs;
  input.nslices = self->nslices;
  input.proc_id = self->proc_id;
  input.thread_id = 0;
  input.slice_id = self->nthreads*self->proc_id;
  energy->energy_terms = self->energy_terms;
  for (i = 0; i < self->nterms+1; i++)
    energy->energy_terms[i] = 0.;
  energy->virial_available = 1;
  energy->error = 0;
  if (energy->force_constants != NULL) {
    input.nthreads = 1;
    input.nprocs = 1;
    input.nslices = 1;
    if (energy->fc_fn != NULL)
      (*energy->fc_fn)(energy, -1, -1, NULL, 0.);
    else {
      double *data = (double *)
	((PyArrayObject *)energy->force_constants)->data;
      long nelements = 9L*(long)natoms*(long)natoms;
      long j;
      for (j = 0; j < nelements; j++)
	data[j] = 0.;
    }
  }
  if (energy->gradients != NULL) {
    if (energy->gradient_fn != NULL) {
#ifdef GRADIENTFN
      (*energy->gradient_fn)(energy, -1, NULL);
#else
      PyErr_SetString(PyExc_EnvironmentError,
		      "gradient function support not available");
      energy->error = 1;
      return;
#endif
    }
    else
    {
      double *data = (double *)((PyArrayObject *)energy->gradients)->data;
      for (i = 0; i < 3*natoms; i++)
	data[i] = 0.;
#ifdef WITH_MPI
      if (self->communicator != NULL && input.nprocs > 1
	  && self->gradient_parts == NULL) {
	self->gradient_parts = malloc(3*natoms*self->nprocs*sizeof(double));
	if (self->gradient_parts == NULL) {
	  energy->error = 1;
	  return;
	}
      }
#endif
    }
  }
#ifdef WITH_THREAD
  {
    threadinfo *tinfo = (threadinfo *)self->scratch;
#if THREAD_DEBUG
    printf("%d threads to be released\n", input.nthreads-1);
#endif
    for (i = 1; i < input.nthreads; i++) {
      tinfo->input.coordinates = coordinates;
      tinfo->input.natoms = natoms;
      tinfo->input.small_change = small_change;
      tinfo->with_gradients = (energy->gradients != NULL);
      if (tinfo->with_gradients && tinfo->energy.gradients == NULL) {
#if defined(NUMPY)
	npy_intp dims[2];
	dims[0] = 3; dims[1] = natoms;
	tinfo->energy.gradients = PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
#else
	int dims[2];
	dims[0] = 3; dims[1] = natoms;
	tinfo->energy.gradients = PyArray_FromDims(2, dims, PyArray_DOUBLE);
#endif
	if (tinfo->energy.gradients == NULL) {
	  energy->error = 1;
	  return;
	}
      }
      if (tinfo->with_gradients) {
	double *data = (double *)
	  ((PyArrayObject *)tinfo->energy.gradients)->data;
	int j;
	for (j = 0; j < 3*natoms; j++)
	  data[j] = 0.;
      }
#if THREAD_DEBUG
      printf("Releasing thread %d\n", tinfo->input.thread_id);
#endif
      PyThread_release_lock(tinfo->lock);
      tinfo++;
    }
  }
#endif
  for (i = 0; i < self->ntermobjects; i++) {
    term = ((PyFFEnergyTermObject **)self->terms->data)[i];
#ifdef WITH_THREAD
    if (term->thread_safe)
      (*term->eval_func)(term, self, &input, energy);
    else {
      PyEval_RestoreThread(self->tstate_save);
      (*term->eval_func)(term, self, &input, energy);
      self->tstate_save = PyEval_SaveThread();
    }
#else
    (*term->eval_func)(term, self, &input, energy);
#endif
#if THREAD_DEBUG
    {
      int j;
      printf("Thread 0: %s: ", term->evaluator_name);
      for (j = term->index; j < term->index+term->nterms; j++)
	printf("%lf ", self->energy_terms[j]);
      printf("\n");
    }
#endif
  }
#ifdef WITH_THREAD
  if (input.nthreads > 1) {
    int done = 1;
#if THREAD_DEBUG
    printf("Collecting data from other threads...\n");
#endif
    while (done < self->nthreads) {
      threadinfo *tinfo = (threadinfo *)self->scratch;
      PyThread_acquire_lock(self->global_lock, 1);
      for (i = 1; i < self->nthreads; i++) {
	if (tinfo->done) {
	  int j;
	  for (j = 0; j < self->nterms+1; j++)
	    energy->energy_terms[j] += tinfo->energy.energy_terms[j];
	  energy->virial_available &= tinfo->energy.virial_available;
	  energy->error |= tinfo->energy.error;
	  if (energy->gradients) {
	    double *data = (double*)((PyArrayObject *)energy->gradients)->data;
	    double *tdata = (double*)((PyArrayObject *)
				      tinfo->energy.gradients)->data;
	    for (i = 0; i < 3*natoms; i++)
	      data[i] += tdata[i];
	  }
	  tinfo->done = 0;
	  done++;
	}
	tinfo++;
      }
      PyThread_release_lock(self->global_lock);
#if THREAD_DEBUG
      if (done > 1)
	printf("%d threads finished\n", done);
#endif
    }
  }
#endif
#if MPI_DEBUG
  {
    int j;
    printf("Proc %d before: ", self->proc_id);
    for (j = 0; j < self->nterms+1; j++)
      printf("%lf ", energy->energy_terms[j]);
    printf("\n");
    fflush(stdout);
  }
#endif
#ifdef WITH_MPI
  if (self->communicator != NULL && input.nprocs > 1) {
    int i, j;
    if (PyMPI_Barrier(self->communicator) != MPI_SUCCESS) {
      PyErr_SetString(PyExc_MPIError, "Error in MPI_Barrier");
      energy->error = 1;
    }
    if (PyMPI_Share(self->communicator, energy->energy_terms,
		    self->energy_parts,
		    PyArray_DOUBLE, self->nterms+1) != MPI_SUCCESS) {
      PyErr_SetString(PyExc_MPIError, "Error in MPI_Share (energy)");
      energy->error = 1;
    }
    for (j = 0; j < self->nterms+1; j++) {
      energy->energy_terms[j] = 0.;
      for (i = 0; i < self->nprocs; i++)
	energy->energy_terms[j] += self->energy_parts[i*(self->nterms+1)+j];
    }
    if (energy->gradients != NULL) {
      double *gdata = (double *)((PyArrayObject *)energy->gradients)->data;
      if (PyMPI_Share(self->communicator, gdata, self->gradient_parts,
		      PyArray_DOUBLE, 3*natoms) != MPI_SUCCESS) {
	PyErr_SetString(PyExc_MPIError, "Error in MPI_Share (gradients)");
	energy->error = 1;
      }
      for (j = 0; j < 3*natoms; j++)
	gdata[j] = 0.;
      for (i = 0; i < self->nprocs; i++)
	for (j = 0; j < 3*natoms; j++)
	  gdata[j] += self->gradient_parts[3*i*natoms+j];
    }
  }
#endif
#if MPI_DEBUG
  {
    int j;
    printf("Proc %d after: ", self->proc_id);
    for (j = 0; j < self->nterms+1; j++)
      printf("%lf ", energy->energy_terms[j]);
    printf("\n");
    fflush(stdout);
  }
#endif
  energy->energy = 0.;
  for (i = 0; i < self->nterms; i++)
    energy->energy += energy->energy_terms[i];
  energy->virial = energy->energy_terms[self->nterms];
}

/* Evaluator loop for parallel threads */

#ifdef WITH_THREAD
static void
evaluator_thread(void *arg)
{
  threadinfo *info = (threadinfo *)arg;
  PyFFEnergyTermObject *term;
  int i;
  while (1) {
#if THREAD_DEBUG
    printf("Thread %d waiting for lock...\n", info->input.thread_id);
#endif
    PyThread_acquire_lock(info->lock, 1);
#if THREAD_DEBUG
    printf("Thread %d running\n", info->input.thread_id);
#endif
    if (info->exit) {
      info->stop = 1;
      break;
    }
    for (i = 0; i < info->evaluator->nterms+1; i++)
      info->energy.energy_terms[i] = 0.;
    info->energy.energy = 0.;
    info->energy.virial_available = 1;
    info->energy.error = 0;
    if (info->with_gradients && info->energy.gradients != NULL) {
      double *data = (double *)((PyArrayObject *)info->energy.gradients)->data;
      for (i = 0; i < 3*info->input.natoms; i++)
	data[i] = 0.;
    }
    PyThread_acquire_lock(info->evaluator->global_lock, 1);
    info->done = 0;
    PyThread_release_lock(info->evaluator->global_lock);
    for (i = 0; i < info->evaluator->ntermobjects; i++) {
      term = ((PyFFEnergyTermObject **)info->evaluator->terms->data)[i];
      if (term->threaded) {
	(*term->eval_func)(term, info->evaluator, &info->input, &info->energy);
#if THREAD_DEBUG
	{
	  int j;
	  printf("Thread %d: %s: ", info->input.thread_id,
		 term->evaluator_name);
	  for (j = term->index; j < term->index+term->nterms; j++)
	    printf("%lf ", info->energy.energy_terms[j]);
	  printf("\n");
	}
#endif
      }
    }
    PyThread_acquire_lock(info->evaluator->global_lock, 1);
    info->done = 1;
    PyThread_release_lock(info->evaluator->global_lock);
  }
}
#endif


/* Return self  */

static PyObject *
C_evaluator(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  Py_INCREF(self);
  return self;
}

/* Documentation string */

static char PyFFEvaluator_Type__doc__[] = 
  "force field evaluator";

/* Method table */

static struct PyMethodDef evaluator_methods[] = {
  {"CEvaluator", C_evaluator, 1},
  {NULL, NULL} /* sentinel */
};

/* Attribute and element access. */

static PyObject *
evaluator_getattr(PyObject *self, char *name)
{
  if (strcmp(name, "last_energy_values") == 0) {
    PyFFEvaluatorObject *ev = (PyFFEvaluatorObject *)self;
    Py_INCREF(ev->energy_terms_array);
    return (PyObject *)ev->energy_terms_array;
  }
  return Py_FindMethod(evaluator_methods, self, name);
}

static Py_ssize_t
evaluator_length(PyFFEvaluatorObject *self)
{
  return self->ntermobjects;
}

static PyObject *
evaluator_item(PyFFEvaluatorObject *self, Py_ssize_t i)
{
  if (i < 0)
    i = self->ntermobjects + i;
  if (i < 0 || i >= self->ntermobjects) {
    PyErr_SetString(PyExc_IndexError, "index out of bounds");
    return NULL;
  }
  Py_INCREF(((PyObject **)self->terms->data)[i]);
  return ((PyObject **)self->terms->data)[i];
}

/* Type object */

static PySequenceMethods PyFFEvaluatorObject_as_sequence = {
  (lenfunc)evaluator_length,   /*sq_length*/
  (binaryfunc)NULL,            /*nb_add*/
  (ssizeargfunc)NULL,            /*nb_multiply*/
  (ssizeargfunc)evaluator_item,  /*sq_item*/
  (ssizessizeargfunc)NULL,     /*sq_slice*/
  (ssizeobjargproc)NULL,	       /*sq_ass_item*/
  (ssizessizeobjargproc)NULL,      /*sq_ass_slice*/
};

PyTypeObject PyFFEvaluator_Type = {
  PyObject_HEAD_INIT(NULL)
  0,			          /*ob_size*/
  "FFEvaluator",	          /*tp_name*/
  sizeof(PyFFEvaluatorObject),    /*tp_basicsize*/
  0,			          /*tp_itemsize*/
  /* methods */
  (destructor)evaluator_dealloc,  /*tp_dealloc*/
  0,			          /*tp_print*/
  (getattrfunc)evaluator_getattr, /*tp_getattr*/
  0,                              /*tp_setattr*/
  0,			          /*tp_compare*/
  0,                              /*tp_repr*/
  0,                              /*tp_as_number*/
  &PyFFEvaluatorObject_as_sequence, /*tp_as_sequence*/
  0,			          /*tp_as_mapping*/
  0,			          /*tp_hash*/
  (ternaryfunc)evaluator_call,	  /*tp_call*/
  0,                              /*tp_str*/
  0,                              /*tp_getattro*/
  0,                              /*tp_setattro*/
  /* Space for future expansion */
  0L,0L,
  /* Documentation string */
  PyFFEvaluator_Type__doc__
};


/* Type "non-bonded list"
 *
 * Objects of this type provide a list of atoms pairs whose distance
 * is below a certain cutoff and which do not appear in the excluded
 * pair list. Pairs that are in the 1-4 list are marked.
 *
 * The structure declaration is in forcefield.h.
 */

/* Allocation and deallocation */

static PyNonbondedListObject *
nblist_new(void)
{
  PyNonbondedListObject *self;
  self = PyObject_NEW(PyNonbondedListObject, &PyNonbondedList_Type);
  if (self == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  self->iterator.state = nblist_start;
  self->iterator.n = -1;
  self->excluded_pairs = NULL;
  self->one_four_pairs = NULL;
  self->atom_subset = NULL;
  self->universe_spec = NULL;
  self->box_number = NULL;
  self->box_atoms = NULL;
  self->boxes = NULL;
  self->nboxes = 0;
  self->allocated_boxes = 0;
  return self;
}

static void
nblist_dealloc(PyNonbondedListObject *self)
{
  Py_XDECREF(self->excluded_pairs);
  Py_XDECREF(self->one_four_pairs);
  Py_XDECREF(self->atom_subset);
  Py_XDECREF(self->universe_spec);
  free(self->box_number);
  free(self->boxes);
  PyObject_Del(self);
}

/* Sequence protocol */

static Py_ssize_t
nblist_length(PyNonbondedListObject *self)
{
  struct nblist_iterator iterator;
  iterator.state = nblist_start;
  while (nblist_iterate(self, &iterator))
    ;
  return iterator.n+1;
}

static PyObject *
nblist_item(PyNonbondedListObject *self, Py_ssize_t i)
{
  if (i < 0) {
    PyErr_SetString(PyExc_IndexError, "index must be positive");
    return NULL;
  }
  if (i < self->iterator.n) {
    self->iterator.state = nblist_start;
    self->iterator.n = -1;
  }
  while (i > self->iterator.n) {
    if (!nblist_iterate(self, &self->iterator)) {
      PyErr_SetString(PyExc_IndexError, "index too large");
      return NULL;
    }
  }
  return Py_BuildValue("ii", self->iterator.a1, self->iterator.a2);
}

/* Methods */

static PyObject *
nblist_update_py(PyObject *self, PyObject *args)
{
  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self;
  PyObject *conf, *geometry;
  PyArrayObject *array;
  double *geometry_data;
  geometry = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &conf, &geometry))
    return NULL;
  if (!PyArray_Check(conf)) {
    geometry = PyObject_GetAttrString(conf, "cell_parameters");
    if (geometry == NULL)
      return NULL;
    conf = PyObject_GetAttrString(conf, "array");
    if (conf == NULL)
      return NULL;
  }
  if (geometry != NULL && !PyArray_Check(geometry)) {
    if (geometry == Py_None)
      geometry = NULL;
    else {
      PyErr_SetString(PyExc_ValueError, "geometry data not an array");
      return NULL;
    }
  }
  if (geometry == NULL)
    geometry_data = nblist->universe_spec->geometry_data;
  else
    geometry_data = (double *)((PyArrayObject *)geometry)->data;
  array = (PyArrayObject *)conf;
  nblist_update(nblist,  array->dimensions[0], (double *)array->data,
		geometry_data);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
nblist_set_cutoff(PyObject *self, PyObject *args)
{
  PyObject *cutoff_ob;
  if (!PyArg_ParseTuple(args, "O", &cutoff_ob))
    return NULL;
  if (cutoff_ob == Py_None)
    ((PyNonbondedListObject *)self)->cutoff = 0.;
  else if (PyNumber_Check(cutoff_ob)) {
    cutoff_ob = PyNumber_Float(cutoff_ob);
    ((PyNonbondedListObject *)self)->cutoff = PyFloat_AsDouble(cutoff_ob);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "cutoff must be a number or None");
    return NULL;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
nblist_pair_distances(PyObject *self, PyObject *args)
{
  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self;
  struct nblist_iterator iterator;
  PyArrayObject *array;
  vector3 dv;
  double *d;
  int i;
#if defined(NUMPY)
  npy_intp n;
#else
  int n;
#endif
  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  n = nblist_length(nblist);
#if defined(NUMPY)
  array = (PyArrayObject *)PyArray_SimpleNew(1, &n, PyArray_DOUBLE);
#else
  array = (PyArrayObject *)PyArray_FromDims(1, &n, PyArray_DOUBLE);
#endif
  if (array == NULL)
    return NULL;
  d = (double *)array->data;
  iterator.state = nblist_start;
  i = 0;
  while (nblist_iterate(nblist, &iterator)) {
    nblist->universe_spec->distance_function(dv, nblist->lastx[iterator.a1],
	nblist->lastx[iterator.a2], nblist->universe_spec->geometry_data);
    d[i++] = vector_length(dv);
  }
  return (PyObject *)array;
}

static PyObject *
nblist_pair_indices(PyObject *self, PyObject *args)
{
  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self;
  struct nblist_iterator iterator;
  PyArrayObject *array;
  long *index;
  int i;
#if defined(NUMPY)
  npy_intp n[2];
#else
  int n[2];
#endif
  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  n[0] = nblist_length(nblist);
  n[1] = 2;
#if defined(NUMPY)
  array = (PyArrayObject *)PyArray_SimpleNew(2, n, PyArray_LONG);
#else
  array = (PyArrayObject *)PyArray_FromDims(2, n, PyArray_LONG);
#endif
  if (array == NULL)
    return NULL;
  index = (long *)array->data;
  iterator.state = nblist_start;
  i = 0;
  while (nblist_iterate(nblist, &iterator)) {
    index[i++] = iterator.a1;
    index[i++] = iterator.a2;
  }
  return (PyObject *)array;
}

/* Documentation string */

static char PyNonbondedList_Type__doc__[] = 
  "nonbonded list";

/* Method table */

static struct PyMethodDef nblist_methods[] = {
  {"update", nblist_update_py, 1},
  {"setCutoff", nblist_set_cutoff, 1},
  {"pairDistances", nblist_pair_distances, 1},
  {"pairIndices", nblist_pair_indices, 1},
  {NULL, NULL} /* sentinel */
};

/* Attribute access. */

static PyObject *
nblist_getattr(PyNonbondedListObject *self, char *name)
{
  return Py_FindMethod(nblist_methods, (PyObject *)self, name);
}

/* Type object */

static PySequenceMethods nblist_as_sequence = {
  (lenfunc)nblist_length,   /*sq_length*/
  0,                        /*nb_add, concat is numeric add*/
  0,                        /*nb_multiply, repeat is numeric multiply*/
  (ssizeargfunc)nblist_item,  /* sq_item*/
  0,                        /*sq_slice*/
  0,                        /*sq_ass_item*/
  0,                        /*sq_ass_slice*/
};

PyTypeObject PyNonbondedList_Type = {
  PyObject_HEAD_INIT(NULL)
  0,			          /*ob_size*/
  "NonbondedList",	          /*tp_name*/
  sizeof(PyNonbondedListObject),  /*tp_basicsize*/
  0,			          /*tp_itemsize*/
  /* methods */
  (destructor)nblist_dealloc,     /*tp_dealloc*/
  0,			          /*tp_print*/
  (getattrfunc)nblist_getattr,    /*tp_getattr*/
  0, 			          /*tp_setattr*/
  0,			          /*tp_compare*/
  0,                              /*tp_repr*/
  0,			          /*tp_as_number*/
  &nblist_as_sequence,            /*tp_as_sequence*/
  0,			          /*tp_as_mapping*/
  0,			          /*tp_hash*/
  0,                              /*tp_call*/
  0,                              /*tp_str*/
  0,                              /*tp_getattro*/
  0,                              /*tp_setattro*/
  /* Space for future expansion */
  0L,0L,
  /* Documentation string */
  PyNonbondedList_Type__doc__
};


/*
 * Force field evaluation function for combined force fields
 */

#if 0
/*
 * Force field evaluation function for Python force fields
 */

static void
python_evaluator(PyFFEnergyTermObject *self,
		 PyFFEvaluatorObject *eval,
		 energy_spec *input,
		 energy_data *energy)
{
  PyObject *py_eval = self->data[0];
  PyObject *args, *result;
  PyObject *gradients, *force_constants;

  gradients = energy->gradients;
  if (gradients == NULL)
    gradients = Py_None;
  force_constants = energy->force_constants;
  if (force_constants == NULL)
    force_constants = Py_None;
  args = PyTuple_New((Py_ssize_t)3);
  Py_INCREF(input->coordinates);
  Py_INCREF(gradients);
  Py_INCREF(force_constants);
  PyTuple_SetItem(args, (Py_ssize_t)0, (PyObject *)input->coordinates);
  PyTuple_SetItem(args, (Py_ssize_t)1, gradients);
  PyTuple_SetItem(args, (Py_ssize_t)2, force_constants);    
  result = PyObject_CallObject(py_eval, args);
  Py_DECREF(args);
  if (result == NULL) {
    energy->error = 1;
  }
  else {
    energy->energy_terms[self->index] = PyFloat_AsDouble(result);
    energy->virial_available = 0;
  }
}
#endif

/*
 * Support functions for pair force fields
 */

void
add_pair_fc(energy_data *energy, int i, int j, vector3 dr,
	    double r_sq, double f1, double f2)
{
  if (energy->fc_fn != NULL) {
    if ((*energy->fc_fn)(energy, i, j, NULL, r_sq)) {
      tensor3 fij;
      int k, l;
      for (k = 0; k < 3; k++) {
	for (l = 0; l < 3; l++)
	  fij[k][l] = (f2-f1)*dr[k]*dr[l]/r_sq;
	fij[k][k] += f1;
      }
      (*energy->fc_fn)(energy, i, i, fij, r_sq);
      (*energy->fc_fn)(energy, j, j, fij, r_sq);
      tensor_changesign(fij);
      if (i > j)
	(*energy->fc_fn)(energy, j, i, fij, r_sq);
      else
	(*energy->fc_fn)(energy, i, j, fij, r_sq);
    }
  }
  else {
    double *data = (double *)((PyArrayObject *)energy->force_constants)->data;
    int n = ((PyArrayObject *)energy->force_constants)->dimensions[0];
    double *fcii = data + 9*n*i+3*i;
    double *fcjj = data + 9*n*j+3*j;
    double *fcij;
    int k, l;
    if (i > j) {
      int temp = i;  i = j;  j = temp;
    }
    fcij = data + 9*n*i+3*j;
    for (k = 0; k < 3; k++) {
      for (l = 0; l < 3; l++) {
	double f = (f2-f1)*dr[k]*dr[l]/r_sq;
	fcii[3*n*k+l] += f;
	fcjj[3*n*k+l] += f;
	fcij[3*n*k+l] -= f;
      }
      fcii[k*(3*n+1)] += f1;
      fcjj[k*(3*n+1)] += f1;
      fcij[k*(3*n+1)] -= f1;
    }
  }
}


/*
 * Constructor functions for the various evaluators
 */

static PyObject *
ListOfNParticleTerms(PyObject *args, ff_eterm_function f,
		     char *eval_name, char *default_term_name)
{
  PyFFEnergyTermObject *self;
  char *name = default_term_name;
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!OO|s",
			&PyUniverseSpec_Type, &self->universe_spec,
			&self->data[0], &self->data[1],
			&name))
    return NULL;
  self->evaluator_name = eval_name;
  self->term_names[0] = allocstring(name);
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);  /* indices */
  Py_INCREF(self->data[1]);  /* parameters */
  self->n = ((PyArrayObject *)self->data[0])->dimensions[0];
  self->nterms = 1;
  self->eval_func = f;
  self->threaded = 1;
  self->thread_safe = 1;
  self->nbarriers = 0;
  self->parallelized = 1;
  return (PyObject *)self;
}

static PyObject *
HarmonicDistanceTerm(PyObject *dummy, PyObject *args)
{
  char *eval_name = "harmonic distance";
  char *default_name = "harmonic bond";
  return ListOfNParticleTerms(args, harmonic_bond_evaluator,
			      eval_name, default_name);
}

static PyObject *
HarmonicAngleTerm(PyObject *dummy, PyObject *args)
{
  char *eval_name = "harmonic angle";
  char *default_name = "harmonic bond angle";
  return ListOfNParticleTerms(args, harmonic_angle_evaluator,
			      eval_name,default_name);
}

static PyObject *
CosineDihedralTerm(PyObject *dummy, PyObject *args)
{
  char *eval_name = "cosine dihedral angle";
  char *default_name = "cosine dihedral angle";
  return ListOfNParticleTerms(args, cosine_dihedral_evaluator,
			      eval_name, default_name);
}

static PyObject *
LennardJonesTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!O!O!O!dd",
			&PyUniverseSpec_Type, &self->universe_spec,
			&PyNonbondedList_Type, &self->data[0],
			&PyArray_Type, &self->data[1],
			&PyArray_Type, &self->data[2],
			&self->param[0], &self->param[1]))
    return NULL;
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  Py_INCREF(self->data[1]);
  Py_INCREF(self->data[2]);
  self->eval_func = lennard_jones_evaluator;
  self->evaluator_name = "Lennard-Jones";
  self->nterms = 0;
  self->thread_safe = 1;
  return (PyObject *)self;
}

static PyObject *
ElectrostaticTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!O!O!dd",
			&PyUniverseSpec_Type, &self->universe_spec,
			&PyNonbondedList_Type, &self->data[0],
			&PyArray_Type, &self->data[1],
			&self->param[0], &self->param[1]))
    return NULL;
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  Py_INCREF(self->data[1]);
  self->eval_func = electrostatic_evaluator;
  self->evaluator_name = "electrostatic";
  self->term_names[0] = allocstring("electrostatic/neutralization");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  self->thread_safe = 1;
  return (PyObject *)self;
}

static PyObject *
EsEwaldTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  PyArrayObject *box_shape;
  long *kmax;
  int natoms, nkvect;
  size_t scratch_size;
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!O!O!O!ddO!dd",
			&PyUniverseSpec_Type, &self->universe_spec,
			&PyArray_Type, &box_shape,
			&PyNonbondedList_Type, &self->data[0],
			&PyArray_Type, &self->data[1],
			&self->param[0], &self->param[3],
			&PyArray_Type, &self->data[2],
			&self->param[1], &self->param[2]))
    return NULL;

  natoms = ((PyArrayObject *)self->data[1])->dimensions[0];
  kmax = (long *)((PyArrayObject *)self->data[2])->data;
  nkvect = init_kvectors(self->universe_spec->box_function,
			 self->universe_spec->geometry_data,
			 natoms,
			 (long *)((PyArrayObject *)self->data[2])->data,
			 self->param[3], NULL, 0);
  scratch_size = natoms*sizeof(vector3) +
    (2*natoms*(kmax[0]+2*kmax[1]+2*kmax[2]+4))*sizeof(double) +
    (3*nkvect+1)*sizeof(int);
  self->scratch = malloc(scratch_size);
  if (self->scratch == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  nkvect = init_kvectors(self->universe_spec->box_function,
			 self->universe_spec->geometry_data,
			 natoms,
			 (long *)((PyArrayObject *)self->data[2])->data,
			 self->param[3], self->scratch, nkvect);

  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  Py_INCREF(self->data[1]);
  Py_INCREF(self->data[2]);
  self->eval_func = es_ewald_evaluator;
  self->thread_safe = 1;
  self->threaded = 1;
  self->nbarriers = 2;
  self->parallelized = 1;
  self->evaluator_name = "electrostatic ewald";
  self->term_names[0] = allocstring("electrostatic/ewald self term");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->term_names[1] = allocstring("electrostatic/ewald reciprocal sum");
  if (self->term_names[1] == NULL)
    return PyErr_NoMemory();
  self->nterms = 2;
  return (PyObject *)self;
}

static PyObject *
NonbondedListTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!", &PyNonbondedList_Type, &self->data[0]))
    return NULL;
  Py_INCREF(self->data[0]);
  self->eval_func = nonbonded_evaluator;
  self->thread_safe = 1;
  self->threaded = 1;
  self->nbarriers = 1;
  self->parallelized = 1;
  self->n = 0;
  self->evaluator_name = "nonbonded list summation";
  self->term_names[0] = allocstring("Lennard-Jones");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->term_names[1] = allocstring("electrostatic/pair sum");
  if (self->term_names[1] == NULL)
    return PyErr_NoMemory();
  self->term_names[2] = allocstring("electrostatic/ewald direct sum");
  if (self->term_names[2] == NULL)
    return PyErr_NoMemory();
  self->nterms = 3;
  return (PyObject *)self;
}

/* Energy evaluator with thread support */

static PyObject *
Evaluator(PyObject *dummy, PyObject *args)
{
  PyFFEvaluatorObject *self = PyFFEvaluator_New();
#ifdef WITH_MPI
  PyMPICommunicatorObject *communicator;
#else
  PyObject *communicator;
#endif
  int nthreads = 1, nbarriers = 0;
  int error = 0;
  int i;
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!|iO",
			&PyArray_Type, &self->terms,
			&nthreads, &communicator))
    return NULL;
  Py_INCREF(self->terms);
  self->eval_func = evaluator;
  self->nthreads = nthreads;
#ifdef WITH_MPI
  if ((PyObject *)communicator == Py_None) {
    self->nprocs = 1;
    self->proc_id = 0;
    self->communicator = NULL;
    self->gradient_parts = NULL;
  }
  else {
    self->nprocs = communicator->size;
    self->proc_id = communicator->rank;
    self->communicator = communicator;
    self->gradient_parts = NULL;
    PyMPI_Barrier(communicator);
  }
#else
  self->nprocs = 1;
  self->proc_id = 0;
#endif
  self->nslices = self->nprocs*self->nthreads;
  self->ntermobjects = self->terms->dimensions[0];
  self->nterms = 0;
  for (i = 0; i < self->ntermobjects; i++) {
    PyFFEnergyTermObject *term = ((PyFFEnergyTermObject **)
				  self->terms->data)[i];
    term->index = self->nterms;
    self->nterms += term->nterms;
    if (term->threaded) {
      term->barrier_index = nbarriers;
      nbarriers += term->nbarriers;
    }
  }
  for (i = 0; i < self->ntermobjects; i++) {
    PyFFEnergyTermObject *term = ((PyFFEnergyTermObject **)
				  self->terms->data)[i];
    term->virial_index = self->nterms;
  }
  self->nterms++;
#if defined(NUMPY)
  {
    npy_intp dims = self->nterms;
    self->energy_terms_array = (PyArrayObject *)
               PyArray_SimpleNew(1, &dims, PyArray_DOUBLE);
  }
#else
  self->energy_terms_array = (PyArrayObject *)
             PyArray_FromDims(1, &self->nterms, PyArray_DOUBLE);
#endif
  self->nterms--;
  if (self->energy_terms_array == NULL) {
    nthreads = 1;
    error = 1;
  }
  else
    self->energy_terms = (double *)self->energy_terms_array->data;
#ifdef WITH_MPI
  self->energy_parts = malloc((self->nterms+1)*self->nprocs*sizeof(double));
  if (self->energy_parts == NULL) {
    self->nprocs = 1;
    error = 1;
  }
#endif
  if (nthreads > 1) {
#ifdef WITH_THREAD
    threadinfo *tinfo;
    int i;
    self->global_lock = PyThread_allocate_lock();
    if (self->global_lock == NULL) {
      PyErr_SetString(PyExc_OSError, "couldn't allocate lock");
      return NULL;
    }
    if (nbarriers > 0) {
      self->binfo = malloc(nbarriers*sizeof(barrierinfo));
      if (self->binfo == NULL)
	return PyErr_NoMemory();
      for (i = 0; i < nbarriers; i++) {
	if (!allocate_barrier(self->binfo+i)) {
	  PyErr_SetString(PyExc_OSError, "couldn't allocate barrier");
	  return NULL;
	}
      }
    }
    self->scratch = malloc((nthreads-1)*sizeof(threadinfo));
    if (self->scratch == NULL)
      return PyErr_NoMemory();
    tinfo = (threadinfo *)self->scratch;
    for (i = 1; i < nthreads; i++) {
      tinfo->evaluator = self;
      tinfo->done = 0;
      tinfo->exit = 0;
      tinfo->stop = 0;
      tinfo->energy.gradients = NULL;
      tinfo->energy.gradient_fn = NULL;
      tinfo->energy.force_constants = NULL;
      tinfo->energy.fc_fn = NULL;
      tinfo->energy.energy_terms =
	       (double *)malloc((self->nterms+1)*sizeof(double));
      if (tinfo->energy.energy_terms == NULL) {
	PyErr_NoMemory();
	error = 1;
      }
      tinfo->lock = NULL;
      tinfo++;
    }
    tinfo = (threadinfo *)self->scratch;
    for (i = 1; i < nthreads; i++) {
      tinfo->lock = PyThread_allocate_lock();
      if (tinfo->lock == NULL) {
	PyErr_SetString(PyExc_OSError, "couldn't allocate lock");
	error = 1;
	break;
      }
      PyThread_acquire_lock(tinfo->lock, 1);
      tinfo++;
    }
    tinfo = (threadinfo *)self->scratch;
    if (!error)
      for (i = 1; i < nthreads; i++) {
	tinfo->input.nthreads = nthreads;
	tinfo->input.thread_id = i;
	tinfo->input.nprocs = self->nprocs;
	tinfo->input.proc_id = self->proc_id;
	tinfo->input.nslices = self->nthreads*self->nprocs;
	tinfo->input.slice_id = self->nthreads*self->proc_id+i;
	if (!PyThread_start_new_thread(evaluator_thread, (void *)tinfo)) {
	  PyErr_SetString(PyExc_OSError, "couldn't start thread");
	  error = 1;
	  break;
	}
#if THREAD_DEBUG
	printf("Thread %d started\n", i);
#endif
	PyThread_acquire_lock(self->global_lock, 1);
	PyThread_release_lock(self->global_lock);
	tinfo++;
      }
#else
    PyErr_SetString(PyExc_OSError, "no thread support");
    error = 1;
#endif
  }
  if (error) {
    evaluator_dealloc(self);
    self = NULL;
  }
  return (PyObject *)self;
}

/*
 * Constructor function for non-bonded lists
 */

static PyObject *
NonbondedList(PyObject *dummy, PyObject *args)
{
  PyNonbondedListObject *self = nblist_new();
  PyObject *cutoff_ob = NULL;
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!O!O!O!O",
			&PyArray_Type, &self->excluded_pairs,
			&PyArray_Type, &self->one_four_pairs,
			&PyArray_Type, &self->atom_subset,
			&PyUniverseSpec_Type, &self->universe_spec,
			&cutoff_ob)) {
    nblist_dealloc(self);
    return NULL;
  }
  if (cutoff_ob == Py_None)
    self->cutoff = 0.;
  else if (PyNumber_Check(cutoff_ob)) {
    cutoff_ob = PyNumber_Float(cutoff_ob);
    self->cutoff = PyFloat_AsDouble(cutoff_ob);
  }
  else {
    PyErr_SetString(PyExc_TypeError, "cutoff must be a number or None");
    nblist_dealloc(self);
    return NULL;
  }
  Py_INCREF(self->excluded_pairs);
  Py_INCREF(self->one_four_pairs);
  Py_INCREF(self->atom_subset);
  Py_INCREF(self->universe_spec);
  return (PyObject *)self;
}

/* C-API access functions for nonbonded lists */

int
PyNonbondedListUpdate(PyNonbondedListObject *nblist,
		      int natoms, double *coordinates,
		      double *geometry_data)
{
  return nblist_update(nblist, natoms, coordinates, geometry_data);
}

int
PyNonbondedListIterate(PyNonbondedListObject *nblist,
		       struct nblist_iterator *iterator)
{
  return nblist_iterate(nblist, iterator);
}

/*
 * List of functions defined in the module
 */

static PyMethodDef forcefield_methods[] = {
  {"HarmonicDistanceTerm", HarmonicDistanceTerm, 1},
  {"HarmonicAngleTerm", HarmonicAngleTerm, 1},
  {"CosineDihedralTerm", CosineDihedralTerm, 1},
  {"LennardJonesTerm", LennardJonesTerm, 1},
  {"ElectrostaticTerm", ElectrostaticTerm, 1},
  {"EsEwaldTerm",  EsEwaldTerm, 1},
  {"NonbondedListTerm", NonbondedListTerm, 1},
  {"Evaluator", Evaluator, 1},
  {"NonbondedList", NonbondedList, 1},
  {"SparseForceConstants", SparseForceConstants, 1},
  {NULL, NULL}		/* sentinel */
};


/* Initialization function for the module */

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initMMTK_forcefield(void)
{
  PyObject *m, *d, *module;
#ifdef WITH_MPI
  PyObject *mpi_module;
#endif
  static void *PyFF_API[PyFF_API_pointers];

  /* Create the module and add the functions */
  m = Py_InitModule("MMTK_forcefield", forcefield_methods);
  
  /* Import the array and MPI modules */
#ifdef import_array
  import_array();
#endif
#ifdef WITH_MPI
  import_mpi();
  mpi_module = PyImport_ImportModule("Scientific.MPI");
  if (mpi_module != NULL) {
    PyObject *module_dict = PyModule_GetDict(mpi_module);
    PyExc_MPIError = PyDict_GetItemString(module_dict, "MPIError");
  }
#endif

  /* Add C API pointer array */ 
  PyFF_API[PyFFEnergyTerm_Type_NUM] = (void *)&PyFFEnergyTerm_Type;
  PyFF_API[PyFFEvaluator_Type_NUM] = (void *)&PyFFEvaluator_Type;
  PyFF_API[PyNonbondedList_Type_NUM] = (void *)&PyNonbondedList_Type;
  PyFF_API[PySparseFC_New_NUM] = (void *)&PySparseFC_New;
  PyFF_API[PySparseFC_Type_NUM] = (void *)&PySparseFC_Type;
  PyFF_API[PySparseFC_Zero_NUM] = (void *)&PySparseFC_Zero;
  PyFF_API[PySparseFC_Find_NUM] = (void *)&PySparseFC_Find;
  PyFF_API[PySparseFC_AddTerm_NUM] = (void *)&PySparseFC_AddTerm;
  PyFF_API[PySparseFC_CopyToArray_NUM] = (void *)&PySparseFC_CopyToArray;
  PyFF_API[PySparseFC_AsArray_NUM] = (void *)&PySparseFC_AsArray;
  PyFF_API[PySparseFC_VectorMultiply_NUM] = (void *)&PySparseFC_VectorMultiply;
  PyFF_API[PySparseFC_Scale_NUM] = (void *)&PySparseFC_Scale;
  PyFF_API[PyFFEnergyTerm_New_NUM] = (void *)&PyFFEnergyTerm_New;
  PyFF_API[PyFFEvaluator_New_NUM] = (void *)&PyFFEvaluator_New;
  PyFF_API[PyNonbondedListUpdate_NUM] = (void *)&PyNonbondedListUpdate;
  PyFF_API[PyNonbondedListIterate_NUM] = (void *)&PyNonbondedListIterate;

#ifdef EXTENDED_TYPES
  if (PyType_Ready(&PyFFEnergyTerm_Type) < 0)
    return;
  if (PyType_Ready(&PyFFEvaluator_Type) < 0)
    return;
  if (PyType_Ready(&PyNonbondedList_Type) < 0)
    return;
  if (PyType_Ready(&PySparseFC_Type) < 0)
    return;
#else
  PyFFEnergyTerm_Type.ob_type = &PyType_Type;
  PyFFEvaluator_Type.ob_type = &PyType_Type;
  PyNonbondedList_Type.ob_type = &PyType_Type;
  PySparseFC_Type.ob_type = &PyType_Type;
#endif

  d = PyModule_GetDict(m);
  PyDict_SetItemString(d, "_C_API",
		       PyCObject_FromVoidPtr((void *)PyFF_API, NULL));
  PyDict_SetItemString(d, "EnergyTerm", (PyObject *)&PyFFEnergyTerm_Type);
  PyDict_SetItemString(d, "EnergyEvaluator",
		       (PyObject *)&PyFFEvaluator_Type);

  /* Get the energy conversion factor from Units */
  module = PyImport_ImportModule("MMTK.Units");
  if (module != NULL) {
    PyObject *module_dict = PyModule_GetDict(module);
    PyObject *factor = PyDict_GetItemString(module_dict,
					    "electrostatic_energy");
    electrostatic_energy_factor = PyFloat_AsDouble(factor);
  }

  /* Get function pointers from _universe */
  module = PyImport_ImportModule("MMTK_universe");
  if (module != NULL) {
    PyObject *module_dict = PyModule_GetDict(module);
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API");
    PyObject *fn;
    if (PyCObject_Check(c_api_object))
      PyUniverse_API = (void **)PyCObject_AsVoidPtr(c_api_object);
    fn = PyDict_GetItemString(module_dict,
			      "infinite_universe_distance_function");
    distance_vector_pointer = (distance_fn *)PyCObject_AsVoidPtr(fn);
    fn = PyDict_GetItemString(module_dict,
			      "orthorhombic_universe_distance_function");
    orthorhombic_distance_vector_pointer =
                        (distance_fn *)PyCObject_AsVoidPtr(fn);
    fn = PyDict_GetItemString(module_dict,
			      "parallelepipedic_universe_distance_function");
    parallelepipedic_distance_vector_pointer =
                        (distance_fn *)PyCObject_AsVoidPtr(fn);
  }

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_forcefield");
}
