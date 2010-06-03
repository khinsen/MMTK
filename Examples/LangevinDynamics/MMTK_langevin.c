/* Low-level Langevin dynamics integrators
 *
 * Written by Konrad Hinsen
 * last revision: 2010-6-3
 */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/trajectory.h"

#include "ranf.h"

/* Global variables */

double kB;   /* Boltzman constant */

/* Allocate and initialize Output variable descriptors */

static PyTrajectoryVariable *
get_data_descriptors(int n, PyUniverseSpecObject *universe_spec,
		     PyArrayObject *configuration, PyArrayObject *velocities,
		     PyArrayObject *gradients, PyArrayObject *masses,
		     int *ndf, double *time,
		     double *p_energy, double *k_energy,
		     double *temperature, double *pressure,
		     double *box_size)
{
  PyTrajectoryVariable *vars = (PyTrajectoryVariable *)
                               malloc((n+1)*sizeof(PyTrajectoryVariable));
  int i = 0;
  if (vars != NULL) {
    if (time != NULL && i < n) {
      vars[i].name = "time";
      vars[i].text = "Time: %lf\n";
      vars[i].unit = time_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Time;
      vars[i].value.dp = time;
      i++;
    }
    if (p_energy != NULL && i < n) {
      vars[i].name = "potential_energy";
      vars[i].text = "Potential energy: %lf, ";
      vars[i].unit = energy_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Energy;
      vars[i].value.dp = p_energy;
      i++;
    }
    if (k_energy != NULL && i < n) {
      vars[i].name = "kinetic_energy";
      vars[i].text = "Kinetic energy: %lf\n";
      vars[i].unit = energy_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Energy;
      vars[i].value.dp = k_energy;
      i++;
    }
    if (temperature != NULL && i < n) {
      vars[i].name = "temperature";
      vars[i].text = "Temperature: %lf\n";
      vars[i].unit = temperature_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Thermodynamic;
      vars[i].value.dp = temperature;
      i++;
    }
    if (pressure != NULL && i < n) {
      vars[i].name = "pressure";
      vars[i].text = "Pressure: %lf\n";
      vars[i].unit = pressure_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Thermodynamic;
      vars[i].value.dp = pressure;
      i++;
    }
    if (configuration != NULL && i < n) {
      vars[i].name = "configuration";
      vars[i].text = "Configuration:\n";
      vars[i].unit = length_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Configuration;
      vars[i].value.array = configuration;
      i++;
    }
    if (box_size != NULL && i < n) {
      vars[i].name = "box_size";
      vars[i].text = "Box size:";
      vars[i].unit = length_unit_name;
      vars[i].type = PyTrajectory_BoxSize;
      vars[i].class = PyTrajectory_Configuration;
      vars[i].value.dp = box_size;
      vars[i].length = universe_spec->geometry_data_length;
      i++;
  }
    if (velocities != NULL && i < n) {
      vars[i].name = "velocities";
      vars[i].text = "Velocities:\n";
      vars[i].unit = velocity_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Velocities;
      vars[i].value.array = velocities;
      i++;
    }
    if (gradients != NULL && i < n) {
      vars[i].name = "gradients";
      vars[i].text = "Energy gradients:\n";
      vars[i].unit = energy_gradient_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Gradients;
      vars[i].value.array = gradients;
      i++;
    }
    if (masses != NULL && i < n) {
      vars[i].name = "masses";
      vars[i].text = "Masses:\n";
      vars[i].unit = mass_unit_name;
      vars[i].type = PyTrajectory_ParticleScalar;
      vars[i].class = PyTrajectory_Internal;
      vars[i].value.array = masses;
      i++;
    }
    if (ndf != NULL && i < n) {
      vars[i].name = "degrees_of_freedom";
      vars[i].text = "Degrees of freedom: %d\n";
      vars[i].unit = "";
      vars[i].type = PyTrajectory_IntScalar;
      vars[i].class = PyTrajectory_Internal;
      vars[i].value.ip = ndf;
      i++;
    }
    vars[i].name = NULL;
  }
  return vars;
}

/* Generate random numbers with bivariate gaussian distribution */
static void
gaussian_bivariate(double *r1, double *r2,
		   double s1, double s2, double c12, double c12x)
{
  double v1, v2, s;
  do {
    v1 = 2.*Ranf()-1.;
    v2 = 2.*Ranf()-1.;
    s = v1*v1 + v2*v2;
  } while (s >= 1. || s == 0.);
  s = sqrt(-2.*log(s)/s);
  v1 *= s;
  v2 *= s;
  *r1 = s1*v1;
  *r2 = s2*(c12*v1+c12x*v2);
}

/* Generate random forces */
static void
random_forces(vector3 *rx, vector3 *rv, double *friction, double *mass,
	      int atoms, double temperature, double delta_t)
{
  int i, j;
  for (i = 0; i < atoms; i++) {
    if (friction[i] < 1.e-8) {
      rx[i][0] = rx[i][1] = rx[i][2] = 0.;
      rv[i][0] = rv[i][1] = rv[i][2] = 0.;
    }
    else {
      double ft = delta_t*friction[i]/mass[i];
      double expft = exp(-ft);
      double ktm = kB*temperature/mass[i];
      double s1 = delta_t*sqrt(ktm*(2.-(3.+(expft-4.)*expft)/ft)/ft);
      double s2 = sqrt(ktm*(1.-sqr(expft)));
      double c12 = delta_t*ktm*(1.+expft*(expft-2.))/(ft*s1*s2);
      double c12x = sqrt(1.-c12*c12);
      for (j = 0; j < 3; j++)
	gaussian_bivariate(&rx[i][j], &rv[i][j], s1, s2, c12, c12x);
    }
  }
}

/* Langevin dynamics integrator */

static PyObject *
integrateLD(PyObject *dummy, PyObject *args)
{
  /* The parameters passed from the Python code */
  PyObject *universe;                /* universe */
  PyArrayObject *configuration;      /* array of positions */
  PyArrayObject *velocities;         /* array of velocities */
  PyArrayObject *masses;             /* array of masses */
  PyArrayObject *friction;           /* array of friction coefficients */
  PyListObject *spec_list;           /* list of periodic actions */
  PyFFEvaluatorObject *evaluator;    /* energy evaluator */
  double ext_temp;                   /* temperature of the heat bath */
  double delta_t;                    /* time step */
  int steps;                         /* number of steps */

  /* Other variables, see below for explanations */
  PyThreadState *this_thread;
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *gradients, *random1, *random2;
  PyTrajectoryOutputSpec *output;
  PyTrajectoryVariable *data_descriptors = NULL;
  vector3 *x, *v, *g, *rx, *rv;
  double *m, *f;
  energy_data p_energy;
  double time, k_energy, temperature, volume, pressure;
  int atoms, df;
  int pressure_available;
  int i, j;

  /* Get arguments passed from Python code */
  if (!PyArg_ParseTuple(args, "OO!O!O!O!O!ddiO!", &universe,
			&PyArray_Type, &configuration,
			&PyArray_Type, &velocities,
			&PyArray_Type, &masses,
			&PyArray_Type, &friction,
			&PyFFEvaluator_Type, &evaluator,
			&ext_temp, &delta_t, &steps,
			&PyList_Type, &spec_list))
    return NULL;

  /* Obtain the universe specification */
  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;

  /* Create the array for energy gradients */
#if defined(NUMPY)
  gradients = (PyArrayObject *)PyArray_Copy(configuration);
#else
  gradients = (PyArrayObject *)PyArray_FromDims(configuration->nd,
						configuration->dimensions,
						PyArray_DOUBLE);
#endif
  if (gradients == NULL)
    return NULL;

  /* Create the arrays for random forces */
#if defined(NUMPY)
  random1 = (PyArrayObject *)PyArray_Copy(configuration);
#else
  random1 = (PyArrayObject *)PyArray_FromDims(configuration->nd,
					      configuration->dimensions,
					      PyArray_DOUBLE);
#endif
  if (random1 == NULL) {
    Py_DECREF(gradients);
    return NULL;
  }
#if defined(NUMPY)
  random2 = (PyArrayObject *)PyArray_Copy(configuration);
#else
  random2 = (PyArrayObject *)PyArray_FromDims(configuration->nd,
					      configuration->dimensions,
					      PyArray_DOUBLE);
#endif
  if (random2 == NULL) {
    Py_DECREF(gradients);
    Py_DECREF(random1);
    return NULL;
  }

  /* Set some convenient variables */
  atoms = configuration->dimensions[0];    /* number of atoms */
  x = (vector3 *)configuration->data;      /* pointer to positions */
  v = (vector3 *)velocities->data;         /* pointer to velocities */
  g = (vector3 *)gradients->data;          /* pointer to energy gradients */
  rx = (vector3 *)random1->data;           /* pointer to random moves */
  rv = (vector3 *)random2->data;           /* pointer to random velocities */
  f = (double *)friction->data;            /* pointer to friction constants */
  m = (double *)masses->data;              /* pointer to masses */
  df = 3*atoms;                            /* number of degrees of freedom */

  /* Initial coordinate correction (for periodic universes etc.) */
  universe_spec->correction_function(x, atoms, universe_spec->geometry_data);

  /* Initial force calculation */
  p_energy.gradients = (PyObject *)gradients;
  p_energy.gradient_fn = NULL;
  p_energy.force_constants = NULL;
  p_energy.fc_fn = NULL;
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, 1);
  (*evaluator->eval_func)(evaluator, &p_energy, configuration, 0);
  PyUniverseSpec_StateLock(universe_spec, 2);
  Py_END_ALLOW_THREADS;
  if (p_energy.error)
    goto error2;

  /* Check if pressure can be calculated (i.e. the force field implementation
     provides the virial and the universe has a finite volume) */
  Py_BEGIN_ALLOW_THREADS;
  PyUniverseSpec_StateLock(universe_spec, 1);
  volume = universe_spec->volume_function(1., universe_spec->geometry_data);
  PyUniverseSpec_StateLock(universe_spec, 2);
  Py_END_ALLOW_THREADS;
  pressure_available = volume > 0. && p_energy.virial_available;

  /* Initialize trajectory output and periodic actions */
  data_descriptors =
    get_data_descriptors(10 + pressure_available,
			 universe_spec,
			 configuration, velocities,
			 gradients, masses, &df,
			 &time, &p_energy.energy, &k_energy,
			 &temperature,
			 pressure_available ? &pressure:NULL,
			 (universe_spec->geometry_data_length > 0) ?
			           universe_spec->geometry_data : NULL);
  if (data_descriptors == NULL)
    goto error2;
  output = PyTrajectory_OutputSpecification(universe, spec_list,
					    "Langevin Dynamics",
					    data_descriptors);
  if (output == NULL)
    goto error2;

  /* Allow parallel threads and get the universe write state lock */
  this_thread = PyEval_SaveThread();
  PyUniverseSpec_StateLock(universe_spec, -1);

  /*
   * Main integration loop
   */
  time  = 0.;
  for (i = 0; i < steps; i++) {

    /* Calculation of thermodynamic properties */
    k_energy = 0.;
    for (j = 0; j < atoms; j++)
      k_energy += m[j]*dot(v[j], v[j]);
    k_energy *= 0.5;
    temperature = 2.*k_energy/(df*kB);
    pressure = 2.*k_energy+p_energy.virial/(3.*volume);

    /* Trajectory and log output */
    if (PyTrajectory_Output(output, i, data_descriptors, &this_thread) == -1) {
      PyUniverseSpec_StateLock(universe_spec, -2);
      PyEval_RestoreThread(this_thread);
      goto error;
    }

    /* Calculation of random forces */
    random_forces(rx, rv, f, m, atoms, ext_temp, delta_t);

    /* First part of integration step */
    for (j = 0; j < atoms; j++) {
      double ft = delta_t*f[j]/m[j];
      double c0 = 1.+ft*(-1.+ft*(0.5+ft*(-1./6.+ft/24.)));
      double c1 = 1.+ft*(-0.5+ft*(1./6.-ft/24.));
      double c2 = 0.5+ft*(-1./6.+ft/24.);
      x[j][0] += delta_t*(c1*v[j][0]-c2*delta_t*g[j][0]/m[j])+rx[j][0];
      x[j][1] += delta_t*(c1*v[j][1]-c2*delta_t*g[j][1]/m[j])+rx[j][1];
      x[j][2] += delta_t*(c1*v[j][2]-c2*delta_t*g[j][2]/m[j])+rx[j][2];
      v[j][0] = c0*v[j][0]+(c2-c1)*delta_t*g[j][0]/m[j]+rv[j][0];
      v[j][1] = c0*v[j][1]+(c2-c1)*delta_t*g[j][1]/m[j]+rv[j][1];
      v[j][2] = c0*v[j][2]+(c2-c1)*delta_t*g[j][2]/m[j]+rv[j][2];
    }

    /* Coordinate correction (for periodic universes etc.) */
    universe_spec->correction_function(x, atoms, universe_spec->geometry_data);

    /* Mid-step force evaluation */
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyUniverseSpec_StateLock(universe_spec, 1);
    (*evaluator->eval_func)(evaluator, &p_energy, configuration, 1);
    PyUniverseSpec_StateLock(universe_spec, 2);
    if (p_energy.error) {
      PyEval_RestoreThread(this_thread);
      goto error;
    }
    PyUniverseSpec_StateLock(universe_spec, -1);

    /* Second part of integration step */
    for (j = 0; j < atoms; j++) {
      double ft = delta_t*f[j]/m[j];
      double c2 = 0.5+ft*(-1./6.+ft/24.);
      v[j][0] -= c2*delta_t*g[j][0]/m[j];
      v[j][1] -= c2*delta_t*g[j][1]/m[j];
      v[j][2] -= c2*delta_t*g[j][2]/m[j];
    }

    /* The End - next time step! */
    time += delta_t;
  }
  /** End of main integration loop **/

  /* Final thermodynamic property evaluation */
  k_energy = 0.;
  for (j = 0; j < atoms; j++)
    k_energy += m[j]*dot(v[j], v[j]);
  k_energy *= 0.5;
  temperature = 2.*k_energy/(df*kB);
  pressure = 2.*k_energy+p_energy.virial/(3.*volume);

  /* Final trajectory and log output */
  if (PyTrajectory_Output(output, i, data_descriptors, &this_thread) == -1) {
    PyUniverseSpec_StateLock(universe_spec, -2);
    PyEval_RestoreThread(this_thread);
    goto error;
  }

  /* Cleanup */
  PyUniverseSpec_StateLock(universe_spec, -2);
  PyEval_RestoreThread(this_thread);
  PyTrajectory_OutputFinish(output, i, 0, 1, data_descriptors);
  free(data_descriptors);
  Py_DECREF(gradients);
  Py_DECREF(random1);
  Py_DECREF(random2);
  Py_INCREF(Py_None);
  return Py_None;

  /* Error return */
error:
  PyTrajectory_OutputFinish(output, i, 1, 1, data_descriptors);
error2:
  free(data_descriptors);
  Py_DECREF(gradients);
  Py_DECREF(random1);
  Py_DECREF(random2);
  return NULL;
}


/*
 * List of functions defined in the module
 */

static PyMethodDef langevin_methods[] = {
  {"integrateLD", integrateLD, 1},
  {NULL, NULL}		/* sentinel */
};

/* Initialization function for the module */

void
initMMTK_langevin()
{
  PyObject *m, *dict;
  PyObject *universe, *trajectory, *forcefield, *units;

  /* Create the module and add the functions */
  m = Py_InitModule("MMTK_langevin", langevin_methods);
  dict = PyModule_GetDict(m);

  /* Import the array module */
  import_array();

  /* Import MMTK modules */
  import_MMTK_universe();
  import_MMTK_forcefield();
  import_MMTK_trajectory();

  /* Get the Boltzman constant factor from MMTK.Units */
  units = PyImport_ImportModule("MMTK.Units");
  if (units != NULL) {
    PyObject *module_dict = PyModule_GetDict(units);
    PyObject *factor = PyDict_GetItemString(module_dict, "k_B");
    kB = PyFloat_AsDouble(factor);
  }

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_langevin");
}
