/* C routines for HarmonicOscillatorFF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
harmonic_evaluator(PyFFEnergyTermObject *self,
                   PyFFEvaluatorObject *eval,
                   energy_spec *input,
                   energy_data *energy)
     /* The four parameters are pointers to structures that are
        defined in MMTK/forcefield.h.
        PyFFEnergyTermObject: All data relevant to this particular
                              energy term.
        PyFFEvaluatorObject:  Data referring to the global energy
                              evaluation process, e.g. parallelization
                              options. Not used here.
        energy_spec:          Input parameters for this routine, i.e.
                              atom positions and parallelization parameters.
        energy_data:          Storage for the results (energy terms,
                              gradients, second derivatives).
     */
{
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  double x = self->param[0];  /* reference point x */
  double y = self->param[1];  /* reference point y */
  double z = self->param[2];  /* reference point z */
  double k = self->param[3];  /* force constant */
  int atom_index = (int)self->param[4];  /* atom index */

  double dx = coordinates[atom_index][0] - x;
  double dy = coordinates[atom_index][1] - y;
  double dz = coordinates[atom_index][2] - z;

  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */
  energy->energy_terms[self->index] = 0.5*k*(dx*dx + dy*dy + dz*dz);
  energy->energy_terms[self->virial_index] -= k*(dx*dx + dy*dy + dz*dz);

  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL && energy->force_constants == NULL)
    return;

  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased. */
  if (energy->gradients != NULL) {
    vector3 *g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
    g[atom_index][0] += k*dx;
    g[atom_index][1] += k*dy;
    g[atom_index][2] += k*dz;
  }

  /* Add the force constant contribution to the global force constant array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased. */
  if (energy->force_constants != NULL) {
    double *fc = (double *)((PyArrayObject*)energy->force_constants)->data;
    int n = ((PyArrayObject *)energy->force_constants)->dimensions[0];
    fc += (9*n+3)*atom_index;
    fc[3*n*0+0] += k;
    fc[3*n*1+1] += k;
    fc[3*n*2+2] += k;
  }
}

/* A utility function that allocates memory for a copy of a string */
static char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* The next function is meant to be called from Python. It creates the
   energy term object at the C level and stores all the parameters in
   there in a form that is convient to access for the C routine above.
   This is the routine that is imported into and called by the Python
   module, HarmonicOscillatorFF.py. */
static PyObject *
HarmonicOscillatorTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int atom_index;
  double x, y, z;
  double force_constant;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!idddd",
                        &PyUniverseSpec_Type, &self->universe_spec,
                        &atom_index, &x, &y, &z, &force_constant))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = harmonic_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "harmonic_oscillator";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("harmonic_oscillator");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = x;
  self->param[1] = y;
  self->param[2] = z;
  self->param[3] = force_constant;
  self->param[4] = (double) atom_index;
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"HarmonicOscillatorTerm", HarmonicOscillatorTerm, 1},
  {NULL, NULL} /* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_harmonic_oscillator(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_harmonic_oscillator", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_harmonic_oscillator");
}
