/* C routines for ElectricField.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
ef_evaluator(PyFFEnergyTermObject *self,
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
  int natoms = input->coordinates->dimensions[0];
  vector3 *g;
  double ex = self->param[0];  /* electric field x */
  double ey = self->param[1];  /* electric field y */
  double ez = self->param[2];  /* electric field z */
  PyArrayObject *charge_array = (PyArrayObject *)self->data[0];
  double *charges = (double *)charge_array->data;  /* atomic charges */
  int atom_index;

  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term,
     which is initialized to zero.
     Note that the virial is also stored in the array energy_terms,
     at the index self->virial_index. However, there is only one virial
     term for the whole system, it is not added up term by term. Therefore
     we don't set it to zero, we just add to it in the loop. */
  energy->energy_terms[self->index] = 0.;

  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased.
     If energy_gradients is NULL, then the calling routine does not
     want gradients, and didn't provide storage for them.
     Second derivatives are not calculated because they are zero. */
  if (energy->gradients != NULL)
    g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
  for (atom_index = 0; atom_index < natoms; atom_index++) {
    double q = charges[atom_index];
    double e = q*(ex*coordinates[atom_index][0]
		  + ey*coordinates[atom_index][1]
		  + ez*coordinates[atom_index][2]);
    energy->energy_terms[self->index] += e;
    energy->energy_terms[self->virial_index] -= e;
    if (energy->gradients != NULL) {
      g[atom_index][0] += q*ex;
      g[atom_index][1] += q*ey;
      g[atom_index][2] += q*ez;
    }
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
   module, ElectricField.py. */
static PyObject *
ElectricFieldTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  PyArrayObject *charges;
  double x, y, z;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!O!ddd",
			&PyUniverseSpec_Type, &self->universe_spec,
			&PyArray_Type, &charges,
			&x, &y, &z))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = ef_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "electric_field";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("electric_field");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there. */
  self->param[0] = x;
  self->param[1] = y;
  self->param[2] = z;
  /* self->data is the other storage area for parameters. There are
     40 Python object slots there */
  self->data[0] = (PyObject *)charges;
  Py_INCREF(charges);
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"ElectricFieldTerm", ElectricFieldTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_electric_field(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_electric_field", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_electric_field");
}
