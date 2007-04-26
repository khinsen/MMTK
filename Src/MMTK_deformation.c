/* Deformation energy calculation.
 *
 * Written by Konrad Hinsen
 * last revision: 2007-4-26
 */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

/* String copy with memory allocation */

char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/*
 * Deformation evaluation
 */
static double
deformation(vector3 *x, vector3 *v, vector3 *g, double *la, int natoms,
            PyNonbondedListObject *nblist,
            double cutoff, double fc_length, double factor, int normalize,
            int version)
{
  int n;
  double cutoff_sq = sqr(cutoff);
  double norm, def;
  struct nblist_iterator iterator;

  if (normalize) {
    norm = 0.;
    for (n = 0; n < natoms; n++)
      norm += v[n][0]*v[n][0] + v[n][1]*v[n][1] + v[n][2]*v[n][2];
    norm = sqrt(norm/natoms);
  }
  else
    norm = 1.;

  if (g != NULL)
    for (n = 0; n < natoms; n++)
      g[n][0] = g[n][1] = g[n][2] = 0.;

  if (la != NULL)
    for (n = 0; n < natoms; n++)
      la[n]= 0.;

  def = 0.;
  iterator.state = nblist_start;
  while (PyNonbondedListIterate(nblist, &iterator)) {
    int i = iterator.a1;
    int j = iterator.a2;
    vector3 rij;
    double r_sq;
    nblist->universe_spec->distance_function(rij,
                             x[iterator.a2], x[iterator.a1],
                             nblist->universe_spec->geometry_data);
    r_sq = vector_length_sq(rij);
    if (r_sq <= cutoff_sq) {
      vector3 vij;
      double k, l, l2;
      vij[0] = v[i][0]-v[j][0];
      vij[1] = v[i][1]-v[j][1];
      vij[2] = v[i][2]-v[j][2];
      if (version == 0)
        /* exponential */
        k = factor * exp((0.01-r_sq)/sqr(fc_length));
      else if (version == 1) {
        /* derived from Amber94 force constant matrix */
        if (r_sq < 0.16) {
          double r = sqrt(r_sq);
          k = (860000.*r-239000.)*factor;
        }
        else
          k = 128.*factor/(r_sq*r_sq*r_sq);
      }
      else
        k = 0.;
      l = dot(rij, vij)/norm;
      l2 = k*l*l/r_sq;
      if (g != NULL) {
        g[i][0] += 2*k*l*rij[0]/(norm*natoms*r_sq);
        g[i][1] += 2*k*l*rij[1]/(norm*natoms*r_sq);
        g[i][2] += 2*k*l*rij[2]/(norm*natoms*r_sq);
        g[j][0] -= 2*k*l*rij[0]/(norm*natoms*r_sq);
        g[j][1] -= 2*k*l*rij[1]/(norm*natoms*r_sq);
        g[j][2] -= 2*k*l*rij[2]/(norm*natoms*r_sq);
      }
      if (la != NULL) {
        la[i] += 0.5*l2;
        la[j] += 0.5*l2;
      }
      def += l2;
    }
  }
  if (g != NULL) {
    if (normalize)
      for (n = 0; n < natoms; n++) {
        g[n][0] = g[n][0] - 2.*def*v[n][0]/(sqr(natoms)*sqr(norm));
        g[n][1] = g[n][1] - 2.*def*v[n][1]/(sqr(natoms)*sqr(norm));
        g[n][2] = g[n][2] - 2.*def*v[n][2]/(sqr(natoms)*sqr(norm));
      }
  }
  return def/natoms;
}

static double
finite_deformation(vector3 *x, vector3 *v, vector3 *g, double *la, int natoms,
                   PyNonbondedListObject *nblist,
                   double cutoff, double fc_length, double factor,
                   int version)
{
  int n;
  double cutoff_sq = sqr(cutoff);
  double def;
  struct nblist_iterator iterator;

  if (g != NULL)
    for (n = 0; n < natoms; n++)
      g[n][0] = g[n][1] = g[n][2] = 0.;

  if (la != NULL)
    for (n = 0; n < natoms; n++)
      la[n]= 0.;

  def = 0.;
  iterator.state = nblist_start;
  while (PyNonbondedListIterate(nblist, &iterator)) {
    int i = iterator.a1;
    int j = iterator.a2;
    vector3 rij;
    double rij2;
    nblist->universe_spec->distance_function(rij,
                             x[iterator.a2], x[iterator.a1],
                             nblist->universe_spec->geometry_data);
    rij2 = vector_length_sq(rij);
    if (rij2 <= cutoff_sq) {
      vector3 rvij;
      double rvij2, rrvij2, k, l, l2;
      rvij[0] = rij[0] + v[i][0]-v[j][0];
      rvij[1] = rij[1] + v[i][1]-v[j][1];
      rvij[2] = rij[2] + v[i][2]-v[j][2];
      if (version == 0)
        /* exponential */
        k = factor * exp((0.01-rij2)/sqr(fc_length));
      else if (version == 1) {
        /* derived from Amber94 force constant matrix */
        if (rij2 < 0.16) {
          double r = sqrt(rij2);
          k = (860000.*r-239000.)*factor;
        }
        else {
          k = 128.*factor/(rij2*rij2*rij2);
        }
      }
      else
        k = 0.;
      rvij2 = vector_length_sq(rvij);
      rrvij2 = sqrt(rvij2);
      l = rrvij2-sqrt(rij2);
      l2 = k*l*l;
      if (g != NULL) {
        double f = 2.*k*l/(natoms*rrvij2);
        g[i][0] += f*rvij[0];
        g[i][1] += f*rvij[1];
        g[i][2] += f*rvij[2];
        g[j][0] -= f*rvij[0];
        g[j][1] -= f*rvij[1];
        g[j][2] -= f*rvij[2];
      }
      if (la != NULL) {
        la[i] += 0.5*l2;
        la[j] += 0.5*l2;
      }
      def += l2;
    }
  }
  return def/natoms;
}

static PyObject *
deformation_py(PyObject *dummy, PyObject *args)
{
  PyArrayObject *configuration;
  PyArrayObject *displacement;
  PyArrayObject *gradient = NULL;
  PyArrayObject *l_atom = NULL;
  PyNonbondedListObject *nblist;

  int natoms;
  vector3 *x, *v, *g;
  double *la;
  double cutoff, fc_length, factor;
  int normalize = 0, finite = 0, version = 0;

  if (!PyArg_ParseTuple(args, "O!O!O!OOdddi|ii",
                        &PyArray_Type, &configuration,
                        &PyArray_Type, &displacement,
                        &PyNonbondedList_Type, &nblist,
                        &gradient, &l_atom,
                        &cutoff, &fc_length, &factor, &normalize,
                        &finite, &version))
    return NULL;
  natoms = configuration->dimensions[0];
  x = (vector3 *)configuration->data;
  v = (vector3 *)displacement->data;
  if ((PyObject *)gradient == Py_None)
    g = NULL;
  else {
    if (PyArray_Check(gradient))
      g = (vector3 *)gradient->data;
    else {
      PyErr_SetString(PyExc_TypeError, "not an array");
      return NULL;
    }
  }
  if ((PyObject *)l_atom == Py_None)
    la = NULL;
  else {
    if (PyArray_Check(l_atom))
      la = (double *)l_atom->data;
    else {
      PyErr_SetString(PyExc_TypeError, "not an array");
      return NULL;
    }
  }
  if (finite)
    return PyFloat_FromDouble(finite_deformation(x, v, g, la, natoms, nblist,
                                                 cutoff, fc_length, factor,
                                                 version));
  else
    return PyFloat_FromDouble(deformation(x, v, g, la, natoms, nblist,
                                          cutoff, fc_length, factor,
                                          normalize, version));
}

/* Deformation reduction */

static void
reduce_deformation(vector3 *x, vector3 *v, vector3 *g, int natoms, int niter,
                   PyNonbondedListObject *nblist,
                   double cutoff, double fc_length, double factor,
                   int version)
{
  struct nblist_iterator iterator;
  double f, min_dist_sq, max_k;
  int i, j;

  min_dist_sq = 1.e30;
  iterator.state = nblist_start;
  while (PyNonbondedListIterate(nblist, &iterator)) {
    vector3 rij;
    double r_sq;
    nblist->universe_spec->distance_function(rij,
                             x[iterator.a2], x[iterator.a1],
                             nblist->universe_spec->geometry_data);
    r_sq = vector_length_sq(rij);
    if (r_sq < min_dist_sq)
      min_dist_sq = r_sq;
  }
  if (version == 0)
    max_k = factor * exp((0.01-min_dist_sq)/sqr(fc_length));
  else if (version == 1) {
    if (min_dist_sq < 0.16) {
      double r = sqrt(min_dist_sq);
      max_k = (860000.*r-239000.)*factor;
    }
    else
      max_k = 128.*factor/(min_dist_sq*min_dist_sq*min_dist_sq);
  }
  f = 0.9/max_k;

  for (i = 0; i < niter; i++) {
    deformation(x, v, g, NULL, natoms, nblist, cutoff, fc_length, factor,
                0, version);
    for (j = 0; j < natoms; j++) {
      v[j][0] -= f*g[j][0];
      v[j][1] -= f*g[j][1];
      v[j][2] -= f*g[j][2];
    }
  }
}

static void
reduce_finite_deformation(vector3 *x, vector3 *v, vector3 *g,
                          int natoms, double rms_reduction,
                          PyNonbondedListObject *nblist,
                          double cutoff, double fc_length, double factor,
                          int version)
{
  double rms_sq, rms_limit, rms_sq_last, rms_grad, scale, step;
  int i;

  rms_sq = 0.;
  for (i = 0; i < natoms; i++)
    rms_sq += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
  rms_limit = sqrt(rms_sq/natoms)-rms_reduction;
  if (rms_limit < 0.)
    rms_limit = 0.;
  rms_limit = natoms*rms_limit*rms_limit;
  step = 0.01;
  while (1) {
    if (rms_sq <= rms_limit)
      break;
    finite_deformation(x, v, g, NULL, natoms, nblist,
                       cutoff, fc_length, factor, version);
    rms_grad = 0.;
    for (i = 0; i < natoms; i++)
      rms_grad += g[i][0]*g[i][0] + g[i][1]*g[i][1] + g[i][2]*g[i][2];
    rms_sq_last = rms_sq;
    while (1) {
      scale = step/sqrt(rms_grad);
      for (i = 0; i < natoms; i++) {
        v[i][0] -= scale*g[i][0];
        v[i][1] -= scale*g[i][1];
        v[i][2] -= scale*g[i][2];
      }
      rms_sq = 0.;
      for (i = 0; i < natoms; i++)
        rms_sq += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (rms_sq > rms_sq_last) {
        for (i = 0; i < natoms; i++) {
          v[i][0] += scale*g[i][0];
          v[i][1] += scale*g[i][1];
          v[i][2] += scale*g[i][2];
        }
        step /= 2;
      }
      else
        break;
    }
    if (fabs(rms_sq - rms_sq_last) < 1.e-14)
      break;
  }
}

static PyObject *
reduce_deformation_py(PyObject *dummy, PyObject *args)
{
  PyArrayObject *configuration;
  PyArrayObject *displacement;
  PyNonbondedListObject *nblist;
  double cutoff, fc_length, factor;
  int niter;
  int version = 0;

  int natoms;
  vector3 *x, *v;

  vector3 *g;

  if (!PyArg_ParseTuple(args, "O!O!O!dddi|i",
                        &PyArray_Type, &configuration,
                        &PyArray_Type, &displacement,
                        &PyNonbondedList_Type, &nblist,
                        &cutoff, &fc_length, &factor, &niter, &version))
    return NULL;
  natoms = configuration->dimensions[0];
  x = (vector3 *)configuration->data;
  v = (vector3 *)displacement->data;

  g = (vector3 *)malloc(natoms*sizeof(vector3));
  if (g == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  reduce_deformation(x, v, g, natoms, niter, nblist,
                     cutoff, fc_length, factor, version);
  free(g);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
reduce_finite_deformation_py(PyObject *dummy, PyObject *args)
{
  PyArrayObject *configuration;
  PyArrayObject *displacement;
  PyNonbondedListObject *nblist;
  double cutoff, fc_length, factor;
  double rms_reduction;
  int version = 0;

  int natoms;
  vector3 *x, *v;

  vector3 *g;

  if (!PyArg_ParseTuple(args, "O!O!O!dddd|i",
                        &PyArray_Type, &configuration,
                        &PyArray_Type, &displacement,
                        &PyNonbondedList_Type, &nblist,
                        &cutoff, &fc_length, &factor, &rms_reduction,
                        &version))
    return NULL;
  natoms = configuration->dimensions[0];
  x = (vector3 *)configuration->data;
  v = (vector3 *)displacement->data;

  g = (vector3 *)malloc(natoms*sizeof(vector3));
  if (g == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  reduce_finite_deformation(x, v, g, natoms, rms_reduction,
                            nblist, cutoff, fc_length, factor, version);
  free(g);

  Py_INCREF(Py_None);
  return Py_None;
}


/* Forcefield evaluator for deformation energy */

static void
pair_term(energy_data *energy,
          int i, int j, vector3 dr, double r_sq, double f2)
{
  if (energy->fc_fn != NULL) {
    if ((*energy->fc_fn)(energy, i, j, NULL, r_sq)) {
      tensor3 fij;
      int k, l;
      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++)
          fij[k][l] = f2*dr[k]*dr[l]/r_sq;
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
        double f = f2*dr[k]*dr[l]/r_sq;
        fcii[3*n*k+l] += f;
        fcjj[3*n*k+l] += f;
        fcij[3*n*k+l] -= f;
      }
    }
  }
}

void
deformation_evaluator(PyFFEnergyTermObject *self,
                      PyFFEvaluatorObject *eval,
                      energy_spec *input,
                      energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  distance_fn *d_fn = self->universe_spec->distance_function;
  double *distance_data = self->universe_spec->geometry_data;

  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self->data[0];
  struct nblist_iterator iterator;
  double fc_length = self->param[0];
  double cutoff_sq = sqr(self->param[1]);
  double scale_factor = self->param[2];
  double one_four_factor = self->param[3];
  int k;

  int states[3] = {nblist_start, nblist_start_excluded, nblist_start_14};
  double factors[3] = {1., -1., -0.5};

  factors[2] = one_four_factor-1.;
  if (energy->force_constants == NULL)
    return;

  for (k = 0; k < 3; k++) {
    iterator.state = states[k];
    while (PyNonbondedListIterate(nblist, &iterator)) {
      double r_sq;
      vector3 rij;
      (*d_fn)(rij, x[iterator.a2], x[iterator.a1], distance_data);
      r_sq = vector_length_sq(rij);
      if (cutoff_sq == 0. || r_sq <= cutoff_sq) {
        double deriv2 = factors[k]*scale_factor
                         * exp((0.01-r_sq)/sqr(fc_length));
        pair_term(energy, iterator.a1, iterator.a2, rij, r_sq, deriv2);
      }
    }
  }
}

static PyObject *
DeformationTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!Odddd",
                        &PyUniverseSpec_Type, &self->universe_spec,
                        &self->data[0], &self->param[0],
                        &self->param[1], &self->param[2], &self->param[3]))
    return NULL;
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  self->eval_func = deformation_evaluator;
  self->evaluator_name = "deformation";
  self->term_names[0] = allocstring("deformation");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  return (PyObject *)self;
}

void
calpha_evaluator(PyFFEnergyTermObject *self,
                 PyFFEvaluatorObject *eval,
                 energy_spec *input,
                 energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  distance_fn *d_fn = self->universe_spec->distance_function;
  double *distance_data = self->universe_spec->geometry_data;

  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self->data[0];
  struct nblist_iterator iterator;
  double cutoff_sq = sqr(self->param[0]);
  int version = (int)self->param[2];
  int k;

  int states[2] = {nblist_start, nblist_start_excluded};
  double factors[2] = {1., -1.};

  if (energy->force_constants == NULL)
    return;

  for (k = 0; k < 2; k++) {
    iterator.state = states[k];
    while (PyNonbondedListIterate(nblist, &iterator)) {
      double r_sq;
      vector3 rij;
      (*d_fn)(rij, x[iterator.a2], x[iterator.a1], distance_data);
      r_sq = vector_length_sq(rij);
      if (cutoff_sq == 0. || r_sq <= cutoff_sq) {
        double deriv2;
        switch (version) {
        case 0:
          /* fitted from CPC trajectory (CHARMM) */
          if (r_sq < 0.16) {
            double r = sqrt(r_sq);
            deriv2 = (2280600.*r-750400.)*self->param[1];
          }
          else {
            deriv2 = 651.*self->param[1]/(r_sq*r_sq*r_sq);
          }
          break;
        case 1:
          /* derived from Amber94 force constant matrix */
          if (r_sq < 0.16) {
            double r = sqrt(r_sq);
            deriv2 = (860000.*r-239000.)*self->param[1];
          }
          else {
            deriv2 = 128.*self->param[1]/(r_sq*r_sq*r_sq);
          }
          break;
        }
        pair_term(energy, iterator.a1, iterator.a2, rij, r_sq,
                  factors[k]*deriv2);
      }
    }
  }
}

static PyObject *
CalphaTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!Oddd",
                        &PyUniverseSpec_Type, &self->universe_spec,
                        &self->data[0], &self->param[0], &self->param[1],
                        &self->param[2]))
    return NULL;
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  self->eval_func = calpha_evaluator;
  self->evaluator_name = "calpha_deformation";
  self->term_names[0] = allocstring("calpha_deformation");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  return (PyObject *)self;
}

void
an_evaluator(PyFFEnergyTermObject *self,
	     PyFFEvaluatorObject *eval,
	     energy_spec *input,
	     energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  distance_fn *d_fn = self->universe_spec->distance_function;
  double *distance_data = self->universe_spec->geometry_data;

  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self->data[0];
  struct nblist_iterator iterator;
  double cutoff_sq = sqr(self->param[0]);
  int k;

  int states[2] = {nblist_start, nblist_start_excluded};
  double factors[2] = {1., -1.};

  if (energy->force_constants == NULL)
    return;

  for (k = 0; k < 2; k++) {
    iterator.state = states[k];
    while (PyNonbondedListIterate(nblist, &iterator)) {
      double r_sq;
      vector3 rij;
      (*d_fn)(rij, x[iterator.a2], x[iterator.a1], distance_data);
      r_sq = vector_length_sq(rij);
      if (cutoff_sq == 0. || r_sq <= cutoff_sq) {
        pair_term(energy, iterator.a1, iterator.a2, rij, r_sq,
                  factors[k]*self->param[1]);
      }
    }
  }
}

static PyObject *
ANTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  if (!PyArg_ParseTuple(args, "O!Odd",
                        &PyUniverseSpec_Type, &self->universe_spec,
                        &self->data[0], &self->param[0], &self->param[1]))
    return NULL;
  Py_INCREF(self->universe_spec);
  Py_INCREF(self->data[0]);
  self->eval_func = an_evaluator;
  self->evaluator_name = "anistropic_network";
  self->term_names[0] = allocstring("anisotropic_network");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  return (PyObject *)self;
}

/*
 * List of functions defined in the module
 */

static PyMethodDef deformation_methods[] = {
  {"deformation", deformation_py, 1},
  {"reduceDeformation", reduce_deformation_py, 1},
  {"reduceFiniteDeformation", reduce_finite_deformation_py, 1},
  {"DeformationTerm", DeformationTerm, 1},
  {"CalphaTerm", CalphaTerm, 1},
  {"ANTerm", ANTerm, 1},
  {NULL, NULL}          /* sentinel */
};

/* Initialization function for the module */

DL_EXPORT(void)
initMMTK_deformation(void)
{
  PyObject *m;

  /* Create the module and add the functions */
  m = Py_InitModule("MMTK_deformation", deformation_methods);

  /* Import the array module */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_deformation");
}
