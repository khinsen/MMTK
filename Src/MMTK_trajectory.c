/*
 * Trajectory objects using netCDF files.
 *
 * Written by Konrad Hinsen
 * last revision: 2006-5-30
 */

#define _TRAJECTORY_MODULE
#include "MMTK/trajectory.h"
#include "MMTK/universe.h"
#include <time.h>
#include <limits.h>

/* Names of standard dimensions */

char *step_number = "step_number";
char *minor_step_number = "minor_step_number";
char *atom_number = "atom_number";
char *xyz = "xyz";
char *box_size_length = "box_size_length";

/* Names of data classes */

typedef struct {
  char *name;
  int number;
} PyTrajectory_DataClassName;

static PyTrajectory_DataClassName class_names[] = {
  {"configuration", PyTrajectory_Configuration},
  {"position", PyTrajectory_Configuration},
  {"positions", PyTrajectory_Configuration},
  {"velocity", PyTrajectory_Velocities},
  {"velocities", PyTrajectory_Velocities},
  {"gradient", PyTrajectory_Gradients},
  {"gradients", PyTrajectory_Gradients},
  {"energy", PyTrajectory_Energy},
  {"energies", PyTrajectory_Energy},
  {"thermodynamic", PyTrajectory_Thermodynamic},
  {"temperature", PyTrajectory_Thermodynamic},
  {"time", PyTrajectory_Time},
  {"auxiliary", PyTrajectory_Auxiliary},
  {NULL, 0}
};


/* Destroy trajectory object */

static void
trajectory_dealloc(PyTrajectoryObject *self)
{
  if (self->file != NULL)
    PyNetCDFFile_Close(self->file);
  Py_XDECREF(self->universe);
  Py_XDECREF(self->index_map);
  Py_XDECREF(self->file);
  Py_XDECREF(self->var_step);
  Py_XDECREF(self->sbuffer);
  Py_XDECREF(self->vbuffer);
  Py_XDECREF(self->box_buffer);
  PyObject_Del(self);
}

/* Close trajectory file */

static int
PyTrajectory_Close(PyTrajectoryObject *trajectory)
{
  int ret = PyNetCDFFile_Close(trajectory->file);
  Py_DECREF(trajectory->file);
  trajectory->file = NULL;
  return ret;
}

static PyObject *
trajectory_close(PyTrajectoryObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ""))
    return NULL;
  if (PyTrajectory_Close(self) == 0) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  else
    return NULL;
}

/* Write remaining data to file */

static int
PyTrajectory_Flush(PyTrajectoryObject *trajectory)
{
  return PyNetCDFFile_Sync(trajectory->file);
}

/* Retrieve variable, or create it if it doesn't exist */

static PyObject *
PyTrajectory_GetVariable(PyTrajectoryObject *trajectory,
			 char *name, int rank, int integer_flag,
			 char *units, int trajectory_flag)
{
  PyNetCDFVariableObject *var;
  char *dimensions[4];
  int nd = 0;
  if (trajectory_flag)
    dimensions[nd++] = step_number;
  if (rank == PyTrajectory_BoxSize)
    dimensions[nd++] = box_size_length;
  else {
    if (rank != PyTrajectory_Scalar)
      dimensions[nd++] = atom_number;
    if (rank == PyTrajectory_ParticleVector)
      dimensions[nd++] = xyz;
  }
  if (trajectory_flag && trajectory->block_size > 1)
    dimensions[nd++] = minor_step_number;
  var = PyNetCDFFile_GetVariable(trajectory->file, name);
  if (var == NULL) {
    char type;
    if (integer_flag)
      type = 'l';
    else if (trajectory->floattype == PyArray_FLOAT)
      type = 'f';
    else
      type = 'd';
    var = PyNetCDFFile_CreateVariable(trajectory->file, name, type,
				      dimensions, nd);
    if (var != NULL && units != NULL)
      PyNetCDFVariable_SetAttribute(var, "units", PyString_FromString(units));
  }
  return (PyObject *)var;
}

/* Read data */

static PyArrayObject *
PyTrajectory_ReadParticleVector(PyTrajectoryObject *trajectory,
				PyObject *variable, int step)
{
  PyNetCDFVariableObject *v = (PyNetCDFVariableObject *)variable;
  PyNetCDFIndex *indices;
  PyArrayObject *data, *ret;
  int dim[2];
  int i;

  indices = PyNetCDFVariable_Indices(v);
  if (indices == NULL)
    return NULL;
  if (trajectory->block_size > 1) {
    int minor = step % trajectory->block_size;
    int major = step / trajectory->block_size;
    indices[0].start = major;
    indices[0].stop = major + 1;
    indices[0].item = 1;
    indices[v->nd-1].start = minor;
    indices[v->nd-1].stop = minor + 1;
    indices[v->nd-1].item = 1;
  }
  else {
    indices[0].start = step;
    indices[0].stop = step + 1;
    indices[0].item = 1;
  }
  data = PyNetCDFVariable_ReadAsArray(v, indices);
  if (data == NULL)
    return NULL;
  if (trajectory->natoms == trajectory->trajectory_atoms
      && data->descr->type_num == PyArray_DOUBLE)
    return data;
  dim[0] = trajectory->natoms;
  dim[1] = 3;
  ret = (PyArrayObject *)PyArray_FromDims(2, dim, PyArray_DOUBLE);
  if (ret == NULL) {
    Py_DECREF(data);
    return NULL;
  }
  if (data->descr->type_num == PyArray_DOUBLE) {
    double *s = (double *)data->data;
    double *d = (double *)ret->data;
    for (i = 0; i < 3*trajectory->trajectory_atoms; i++)
      *d++ = *s++;
    for (; i < 3*trajectory->natoms; i++)
      *d++ = undefined;
  }
  else {
    float *s = (float *)data->data;
    double *d = (double *)ret->data;
    for (i = 0; i < 3*trajectory->trajectory_atoms; i++)
      *d++ = (double)*s++;
    for (; i < 3*trajectory->natoms; i++)
      *d++ = undefined;
  }
  Py_DECREF(data);
  return ret;
}

static PyArrayObject *
trajectory_read_particle_vector(PyTrajectoryObject *self, PyObject *args)
{
  PyObject *variable;
  char *var_name;
  int step;
  if (!PyArg_ParseTuple(args, "si", &var_name, &step))
    return NULL;
  variable = PyDict_GetItemString(self->file->variables, var_name);
  if (variable == NULL)
    return NULL;
  return PyTrajectory_ReadParticleVector(self, variable, step);
}

static PyArrayObject *
PyTrajectory_ReadParticleTrajectories(PyTrajectoryObject *trajectory, int atom,
				      int natoms, char *variable,
				      int first, int last, int skip,
				      int correct, int box_coordinates)
{
  PyNetCDFVariableObject *var;
  PyNetCDFVariableObject *box_var;
  PyNetCDFIndex *indices, *box_indices;
  PyArrayObject *data, *box;
  PyUniverseSpecObject *universe_spec;
  int bs = trajectory->block_size;
  int steps = (last-first+skip-1)/skip;
  int i, j, k;

  universe_spec = (PyUniverseSpecObject *)
                  PyObject_GetAttrString(trajectory->universe, "_spec");
  if (universe_spec == NULL)
    return NULL;
  var = (PyNetCDFVariableObject *)
         PyDict_GetItemString(trajectory->file->variables, variable);
  if (var == NULL) {
    PyErr_SetString(PyExc_ValueError, "variable not in trajectory");
    return NULL;
  }
  box_var = NULL;
  box = NULL;
  if (universe_spec->is_periodic && (correct || box_coordinates)) {
    box_var = (PyNetCDFVariableObject *)
         PyDict_GetItemString(trajectory->file->variables, "box_size");
    if (box_var != NULL) {
      if (trajectory->box_buffer != NULL) {
	if (trajectory->box_buffer_first == first &&
	    trajectory->box_buffer_last == last &&
	    trajectory->box_buffer_skip == skip) {
	  box = trajectory->box_buffer;
	  Py_INCREF(box);
	}
	else {
	  Py_DECREF(trajectory->box_buffer);
	  trajectory->box_buffer = NULL;
	}
      } 
    }
    else
      PyErr_Clear();
  }
  if (bs > 1) {
    PyArrayObject *read;
    double *data_data, *box_data;
    int bfirst, bsteps, rsteps, istep;
    int dim[3];
    dim[0] = natoms;
    dim[1] = steps;
    dim[2] = 3;
    data = (PyArrayObject *)PyArray_FromDims(3, dim, PyArray_DOUBLE);
    if (data == NULL)
      return NULL;
    data_data = (double *)data->data;
    if (box_var != NULL && box == NULL) {
      dim[0] = steps;
      dim[1] = universe_spec->geometry_data_length;
      box = (PyArrayObject *)PyArray_FromDims(2, dim, PyArray_DOUBLE);
      if (box == NULL) {
	Py_DECREF(data);
	return NULL;
      }
      box_data = (double *)box->data;
    }
    else
      box_data = NULL;
    istep = first;
    rsteps = steps;
    while (istep < last) {
      bfirst = istep % bs;
      bsteps = (bs-bfirst+skip-1)/skip;
      if (bsteps > rsteps)
	bsteps = rsteps;
      indices = PyNetCDFVariable_Indices(var);
      if (indices == NULL)
	return NULL;
      indices[0].start = istep/bs;
      indices[0].stop = istep/bs+1;
      indices[0].item = 1;
      indices[1].start = atom;
      indices[1].stop = atom+natoms;
      indices[3].start = bfirst;
      indices[3].stop = bfirst+skip*bsteps;
      indices[3].stride = skip;
      read = PyNetCDFVariable_ReadAsArray(var, indices);
      if (read == NULL) {
	Py_DECREF(data);
	return NULL;
      }
      if (read->descr->type_num != PyArray_DOUBLE) {
	float *read_data = (float *)read->data;
	for (k = 0; k < natoms; k++)
	  for (j = 0; j < 3; j++)
	    for (i = 0; i < bsteps; i++)
	      data_data[(k*steps+i)*3+j] = (double)read_data[(k*3+j)*bsteps+i];
      }
      else {
	double *read_data = (double *)read->data;
	for (k = 0; k < natoms; k++)
	  for (i = 0; i < bsteps; i++)
	    for (j = 0; j < 3; j++)
	      data_data[(k*steps+i)*3+j] = read_data[(k*3+j)*bsteps+i];
      }
      Py_DECREF(read);
      data_data += 3*bsteps;
      if (box_data != NULL) {
	int dl = universe_spec->geometry_data_length;
	box_indices = PyNetCDFVariable_Indices(box_var);
	if (box_indices == NULL)
	  return NULL;
	box_indices[0].start = istep/bs;
	box_indices[0].stop = istep/bs+1;
	box_indices[0].item = 1;
	box_indices[2].start = bfirst;
	box_indices[2].stop = bfirst+skip*bsteps;
	box_indices[2].stride = skip;
	read = PyNetCDFVariable_ReadAsArray(box_var, box_indices);
	if (read == NULL) {
	  Py_DECREF(data);
	  return NULL;
	}
	if (read->descr->type_num != PyArray_DOUBLE) {
	  float *read_data = (float *)read->data;
	  for (i = 0; i < bsteps; i++)
	    for (j = 0; j < dl; j++)
	      box_data[dl*i+j] = (double)read_data[bsteps*j+i];
	}
	else {
	  double *read_data = (double *)read->data;
	  for (i = 0; i < bsteps; i++)
	    for (j = 0; j < dl; j++)
	      box_data[dl*i+j] = read_data[bsteps*j+i];
	}
	Py_DECREF(read);
	box_data += dl*bsteps;
      }
      rsteps -= bsteps;
      istep += bsteps*skip;
    }
  }
  else {
    indices = PyNetCDFVariable_Indices(var);
    if (indices == NULL)
      return NULL;
    indices[0].start = first;
    indices[0].stop = last;
    indices[0].stride = skip;
    indices[1].start = atom;
    indices[1].stop = atom+natoms;
    data = PyNetCDFVariable_ReadAsArray(var, indices);
    if (data == NULL)
      return NULL;
    if (natoms == 1) {
      data->dimensions[0] = 1;
      data->dimensions[1] = steps;
      data->dimensions[2] = 3;
      if (data->descr->type_num != PyArray_DOUBLE) {
	PyArrayObject *discard = data;
	data = (PyArrayObject *)
	  PyArray_ContiguousFromObject((PyObject *)data,
				       PyArray_DOUBLE, 3, 3);
	Py_DECREF(discard);
	if (data == NULL)
	  return NULL;
      }
    }
    else {
      int dim[3];
      double *d1;
      PyArrayObject *discard = data;
      dim[0] = natoms;
      dim[1] = steps;
      dim[2] = 3;
      data = (PyArrayObject *)PyArray_FromDims(3, dim, PyArray_DOUBLE);
      if (data == NULL) {
	Py_DECREF(discard);
	return NULL;
      }
      d1 = (double *)data->data;
      if (discard->descr->type_num == PyArray_DOUBLE) {
	double *d2 = (double *)discard->data;
	for (k = 0; k < natoms; k++)
	  for (j = 0; j < steps; j++)
	    for (i = 0; i < 3; i++)
	      d1[3*(steps*k+j)+i] = d2[3*(natoms*j+k)+i];
      }
      else {
	float *d2 = (float *)discard->data;
	for (k = 0; k < natoms; k++)
	  for (j = 0; j < steps; j++)
	    for (i = 0; i < 3; i++)
	      d1[3*(steps*k+j)+i] = (double)d2[3*(natoms*j+k)+i];
      }
      Py_DECREF(discard);
    }
    if (box_var != NULL && box == NULL) {
      box_indices = PyNetCDFVariable_Indices(box_var);
      if (box_indices == NULL)
	return NULL;
      box_indices[0].start = first;
      box_indices[0].stop = last;
      box_indices[0].stride = skip;
      box = PyNetCDFVariable_ReadAsArray(box_var, box_indices);
      if (box == NULL) {
	Py_DECREF(data);
	return NULL;
      }
      if (box->descr->type_num != PyArray_DOUBLE) {
	PyArrayObject *discard = box;
	box = (PyArrayObject *)
	  PyArray_ContiguousFromObject((PyObject *)box,
				       PyArray_DOUBLE, 2, 2);
	Py_DECREF(discard);
	if (box == NULL) {
	  Py_DECREF(data);
	  return NULL;
	}
      }
    }
  }
  if (box_var != NULL && trajectory->box_buffer == NULL) {
    trajectory->box_buffer = box;
    Py_INCREF(box);
    trajectory->box_buffer_first = first;
    trajectory->box_buffer_last = last;
    trajectory->box_buffer_skip = skip;
  }
  if (universe_spec->is_periodic && (correct || box_coordinates)) {
    vector3 *t = (vector3 *)data->data;
    if (box_var == NULL)
      universe_spec->box_function(t, t, steps*natoms,
				  universe_spec->geometry_data, 1);
    else {
      for (i = 0; i < natoms; i++)
	universe_spec->trajectory_function(t+i*steps, t+i*steps,
					   steps, (double *)box->data, 1);
    }
    if (correct)
      for (i = 0; i < natoms; i++) {
	for (j = 1; j < steps; j++) {
	  int ai = i*steps+j;
	  double dx = t[ai][0]-t[ai-1][0];
	  double dy = t[ai][1]-t[ai-1][1];
	  double dz = t[ai][2]-t[ai-1][2];
	  while (dx > 0.5) {
	    t[ai][0] -= 1.;
	    dx -= 1.;
	  }
	  while (dx < -0.5) {
	    t[ai][0] += 1.;
	    dx += 1.;
	  }
	  while (dy > 0.5) {
	    t[ai][1] -= 1.;
	    dy -= 1.;
	  }
	  while (dy < -0.5) {
	    t[ai][1] += 1.;
	    dy += 1.;
	  }
	  while (dz > 0.5) {
	    t[ai][2] -= 1.;
	    dz -= 1.;
	  }
	  while (dz < -0.5) {
	    t[ai][2] += 1.;
	    dz += 1.;
	  }
	}
      }
    if (!box_coordinates) {
      if (box_var == NULL) {
	universe_spec->box_function(t, t, steps*natoms,
				    universe_spec->geometry_data, 0);
      }
      else {
	for (i = 0; i < natoms; i++)
	  universe_spec->trajectory_function(t+i*steps, t+i*steps,
					     steps, (double *)box->data, 0);
      }
    }
    if (box_var != NULL) {
      Py_DECREF(box);
    }
  }
  return data;
}

static PyArrayObject *
trajectory_read_particle_trajectories(PyTrajectoryObject *self, PyObject *args)
{
  int atom, natoms, first, last, skip, correct=0, box=0;
  char *variable;
  if (!PyArg_ParseTuple(args, "iisiii|ii", &atom, &natoms, &variable,
			&first, &last, &skip, &correct, &box))
    return NULL;
  return PyTrajectory_ReadParticleTrajectories(self, atom, natoms, variable,
					       first, last, skip,
					       correct, box);
}

static PyArrayObject *
trajectory_read_particle_scalar(PyTrajectoryObject *self, PyObject *args)
{
  PyObject *variable;
  char *var_name;
  int step;
  if (!PyArg_ParseTuple(args, "si", &var_name, &step))
    return NULL;
  variable = PyDict_GetItemString(self->file->variables, var_name);
  if (variable == NULL)
    return NULL;
  PyErr_SetString(PyExc_SystemError, "not yet implemented");
  return NULL;
}

/* Write data */

static int
PyTrajectory_WriteArray(PyTrajectoryObject *trajectory, PyObject *variable,
			PyArrayObject *value)
{
  PyNetCDFVariableObject *v = (PyNetCDFVariableObject *)variable;
  PyNetCDFIndex *indices;
  if (trajectory->write) {
    indices = PyNetCDFVariable_Indices(v);
    if (indices == NULL)
      return 0;
    if (trajectory->block_size > 1) {
      int step = trajectory->steps-1;
      int minor = step % trajectory->block_size;
      int major = step / trajectory->block_size;
      indices[0].start = major;
      indices[0].stop = major + 1;
      indices[0].item = 1;
      indices[v->nd-1].start = minor;
      indices[v->nd-1].stop = minor + 1;
      indices[v->nd-1].item = 1;
    }
    else {
      indices[0].start = trajectory->steps-1;
      indices[0].stop = indices[0].start + 1;
      indices[0].item = 1;
    }
    return PyNetCDFVariable_WriteArray(v, indices, (PyObject *)value);    ;
  }
  else
    return 0;
}

static int
PyTrajectory_WriteFloats(PyTrajectoryObject *trajectory, PyObject *variable,
			 double *values, int n)
{
  static PyArrayObject *a[2] = {NULL, NULL};
  static int last_n[2] = {0, 0};
  int type = (trajectory->floattype == PyArray_DOUBLE) ? 1 : 0;

  if (last_n[type] != n) {
    Py_XDECREF(a[type]);
    a[type] = NULL;
  }
  if (a[type] == NULL) {
    a[type] = (PyArrayObject *)PyArray_FromDims((n == 1) ? 0 : 1, &n,
						trajectory->floattype);
    if (a[type] == NULL)
      return -1;
    last_n[type] = n;
  }

  if (trajectory->floattype == PyArray_DOUBLE) {
    double *data = (double *)(a[type]->data);
    int i;
    for (i = 0; i < n; i++)
      data[i] = values[i];
  }
  else {
    float *data = (float *)(a[type]->data);
    int i;
    for (i = 0; i < n; i++)
      data[i] = (float)values[i];
  }

  return PyTrajectory_WriteArray(trajectory, variable, a[type]);
}

static int
PyTrajectory_WriteInteger(PyTrajectoryObject *trajectory, PyObject *variable,
			  long value)
{
  static PyArrayObject *a = NULL;
  if (a == NULL) {
    int n = 1;
    a = (PyArrayObject *)PyArray_FromDims(0, &n, PyArray_LONG);
    if (a == NULL)
      return -1;
  }
  *(long *)(a->data) = value;
  return PyTrajectory_WriteArray(trajectory, variable, a);
}

/* Set step number */

static int
PyTrajectory_Step(PyTrajectoryObject *trajectory, int step)
{
  if (trajectory->cycle > 0) {
    trajectory->write = 1;
    trajectory->steps = (trajectory->steps % trajectory->cycle) + 1;
    return PyTrajectory_WriteInteger(trajectory,
				     (PyObject *)trajectory->var_step,
				     step);
  }
  else if (step >= trajectory->first_step) {
    trajectory->first_step = 1;
    trajectory->write = 1;
    trajectory->steps++;
    return PyTrajectory_WriteInteger(trajectory,
				     (PyObject *)trajectory->var_step,
				     step);
  }
  else {
    trajectory->write = 0;
    return 0;
  }
}

/* Method table */

static PyMethodDef trajectory_methods[] = {
  {"readParticleVector", (PyCFunction)trajectory_read_particle_vector, 1},
  {"readParticleScalar", (PyCFunction)trajectory_read_particle_scalar, 1},
  {"readParticleTrajectories",
                   (PyCFunction)trajectory_read_particle_trajectories, 1},
  {"close", (PyCFunction)trajectory_close, 1},
  {NULL, NULL}		/* sentinel */
};

/* Attribute access */

static PyObject *
getattr(PyTrajectoryObject *self, char *name)
{
  if (self->file == NULL) {
    PyErr_SetString(PyExc_ValueError, "access to closed trajectory");
    return NULL;
  }
  if (strcmp(name, "file") == 0) {
    Py_INCREF(self->file);
    return (PyObject *)self->file;
  }
  else if (strcmp(name, "nsteps") == 0) {
    return PyInt_FromLong((long)self->steps);
  }
  else if (strcmp(name, "recently_read_box_size") == 0) {
    if (self->box_buffer == NULL) {
    PyErr_SetString(PyExc_AttributeError, "no box size information");
    return NULL;
    }
    Py_INCREF(self->box_buffer);
    return (PyObject *)self->box_buffer;
  }
  else
    return Py_FindMethod(trajectory_methods, (PyObject *)self, name);
}

static int
PyTrajectory_SetAttribute(PyTrajectoryObject *self, char *name,
			  PyObject *value)
{
  return PyNetCDFFile_SetAttribute(self->file, name, value);
}

/* Type definition */

PyTypeObject PyTrajectory_Type = {
  PyObject_HEAD_INIT(NULL)
  0,		/*ob_size*/
  "Trajectory",	/*tp_name*/
  sizeof(PyTrajectoryObject),	/*tp_basicsize*/
  0,		/*tp_itemsize*/
  /* methods */
  (destructor)trajectory_dealloc, /*tp_dealloc*/
  0,			/*tp_print*/
  (getattrfunc)getattr, /*tp_getattr*/
  (setattrfunc)PyTrajectory_SetAttribute, /*tp_setattr*/
  0,			/*tp_compare*/
  0,                    /*tp_repr*/
  0,			/*tp_as_number*/
  0,			/*tp_as_sequence*/
  0,			/*tp_as_mapping*/
  0,			/*tp_hash*/
};

/* Add a time stamp */
static int
PyTrajectory_TimeStamp(PyTrajectoryObject *self, char *text)
{
  time_t now = time(NULL);
  static char time_stamp[200];
  sprintf(time_stamp, text, ctime(&now));
  time_stamp[strlen(time_stamp)-1] = '\0';
  return PyNetCDFFile_AddHistoryLine(self->file, time_stamp);
}

static char *
skip_token(char *p)
{
  if (*p == '\'' || *p == '"') {
    char delimiter = *p++;
    while (*p && *p != delimiter) {
      if (*p == '\\')
	p += 2;
      else
	p++;
    }
    if (*p)
      p++;
  }
  else
    p++;
  return p;
}

/* Compare two universe description strings */
static char *
skip_object(char *p)
{
  int parens = 0;
  if (*p == '\'' || *p == '"')
    p = skip_token(p);
  else {
    while (*p && *p != '(')
      p = skip_token(p);
    while (*p) {
      if (*p == '(')
	parens++;
      else if (*p == ')') {
	parens--;
	if (parens == 0)
	  break;
      }
      p = skip_token(p);
    }
    while (*p && *p != ',')
      p = skip_token(p);
  }
  while (*p && (*p == ',' || *p == ' '))
    p = skip_token(p);
  return p;
}

static int
verify_description(char *d1, char *d2)
{
  char *p1 = d1;
  char *p2 = d2;

  while (*p1 && *p1 != '[')
    p1++;
  while (*p1 && (*p1 == '[' || *p1 == ' '))
    p1++;
  while (*p2 && *p2 != '[')
    p2++;
  while (*p2 && (*p2 == '[' || *p2 == ' '))
    p2++;
  while (*p1 && *p2) {
    char *e1, *e2;
    int n1, n2;
    while (*p1 &&
	   ((*p1 == 'o' && *(p1+1) == '(') ||
	    *p1 == '\'' || *p1 == '"'))
      p1 = skip_object(p1);
    while (*p2 &&
	   ((*p2 == 'o' && *(p2+1) == '(') ||
	    *p2 == '\'' || *p2 == '"'))
      p2 = skip_object(p2);
    if (!*p1 || !*p2)
      break;
    if (*p1 == ']' || *p2 == ']')
      return *p1 == ']' && *p2 == ']';
    e1 = skip_object(p1);
    n1 = e1-p1;
    e2 = skip_object(p2);
    n2 = e2-p2;
    if (n1 != n2)
      return 0;
    if (strncmp(p1, p2, n1) != 0)
      return 0;
    p1 = e1;
    p2 = e2;
  }
  return 1;
}

/* Create a trajectory object */

static PyTrajectoryObject *
PyTrajectory_Open(PyObject *universe, PyObject *description,
		  PyArrayObject *index_map,
		  char *filename, char *mode, int floattype, int cycle,
		  int block_size)
{
  PyTrajectoryObject *self;
  PyUniverseSpecObject *universe_spec;
  PyObject *n_ob;
  PyNetCDFVariableObject *description_var;

  if (!PyObject_HasAttrString(universe, "is_universe")) {
    PyErr_SetString(PyExc_TypeError, "not a universe object");
    return NULL;
  }
  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;
  if (!PyString_Check(description)) {
    PyErr_SetString(PyExc_TypeError, "system description not a string");
    return NULL;
  }
  n_ob = PyObject_CallMethod((PyObject *)universe, "numberOfAtoms", NULL);
  if (n_ob == NULL)
    return NULL;

  self = PyObject_NEW(PyTrajectoryObject, &PyTrajectory_Type);
  if (self == NULL)
    return NULL;
  self->natoms = PyInt_AsLong(n_ob);
  Py_DECREF(n_ob);
  self->universe = universe;
  Py_INCREF(self->universe);
  self->index_map = index_map;
  Py_XINCREF(self->index_map);
  self->var_step = NULL;
  self->sbuffer = NULL;
  self->vbuffer = NULL;
  self->box_buffer = NULL;
  self->floattype = floattype;
  self->trajectory_atoms = (index_map == NULL) ? self->natoms
                                               : index_map->dimensions[0];
  self->cycle = cycle;
  if (cycle > 0)
    block_size = 1;
  self->block_size = block_size;
  self->file = PyNetCDFFile_Open(filename, mode);
  if (self->file == NULL) {
    goto error;
  }
  Py_INCREF(self->file);

  if ((self->index_map != NULL || floattype != PyArray_DOUBLE)
       && mode[0] != 'r') {
    int dim[2];
    if (self->index_map != NULL)
      dim[0] = self->index_map->dimensions[0];
    else
      dim[0] = self->natoms;
    dim[1] = 3;
    self->vbuffer = (PyArrayObject *)PyArray_FromDims(2, dim, floattype);
    if (self->vbuffer == NULL)
      goto error;
    self->sbuffer = (PyArrayObject *)PyArray_FromDimsAndData(1, dim,
					 floattype, self->vbuffer->data);
    if (self->sbuffer == NULL)
      goto error;
  }

  self->first_step = 0;
  self->var_step = PyNetCDFFile_GetVariable(self->file, "step");
  if (self->var_step == NULL) {
    if (mode[0] == 'r') {
      PyErr_SetString(PyExc_TypeError, "not a trajectory file");
      goto error;
    }
    else {
      char *description_length = "description_length";
      Py_ssize_t len;
      int ret;
      PyNetCDFFile_SetAttributeString(self->file, "Conventions", "MMTK/Trajectory");
      if (block_size == 1)
	PyNetCDFFile_SetAttribute(self->file, "trajectory_type",
				  PyInt_FromLong(0));
      else
	PyNetCDFFile_SetAttribute(self->file, "trajectory_type",
				  PyInt_FromLong(1));
      PyTrajectory_TimeStamp(self, "Created %s");
      if (PyNetCDFFile_CreateDimension(self->file, step_number, cycle) == -1)
	goto error;
      if (block_size > 1
	  && PyNetCDFFile_CreateDimension(self->file, minor_step_number,
					  block_size) == -1)
	goto error;
      if (PyNetCDFFile_CreateDimension(self->file, atom_number,
				       self->trajectory_atoms) == -1)
	goto error;
      if (PyNetCDFFile_CreateDimension(self->file, xyz, 3) == -1)
	goto error;
      if (universe_spec->geometry_data_length > 0
	  && PyNetCDFFile_CreateDimension(self->file, box_size_length,
					  universe_spec->geometry_data_length)
	  == -1)
	goto error;
      len = PyString_Size(description);
      if (len > INT_MAX) {
	PyErr_SetString(PyExc_ValueError, "description string too long");
	goto error;
      }
      if (PyNetCDFFile_CreateDimension(self->file, description_length,
				       (int)len) == -1)
	goto error;
      description_var = PyNetCDFFile_CreateVariable(self->file,
						    "description", 'c',
						    &description_length, 1);
      if (description_var == NULL)
	goto error;
      ret = PyNetCDFVariable_WriteString(description_var,
					 (PyStringObject *)description);
      Py_DECREF(description_var);
      if (ret == -1)
	goto error;
      if (self->block_size > 1) {
	char *dim_names[2];
	dim_names[0] = step_number;
	dim_names[1] = minor_step_number;
	self->var_step = PyNetCDFFile_CreateVariable(self->file, "step", 'l',
						     dim_names, 2);
      }
      else
	self->var_step = PyNetCDFFile_CreateVariable(self->file, "step", 'l',
						     &step_number, 1);
      if (self->var_step == NULL)
	goto error;
      if (self->cycle > 0
	  && PyTrajectory_SetAttribute(self, "last_step", PyInt_FromLong(-1))
	     == -1)
	  goto error;
      self->steps = 0;
    }
  }
  else {
    size_t *step_shape = PyNetCDFVariable_GetShape(self->var_step);
    Py_INCREF(self->var_step);
    if (PyNetCDFVariable_GetRank(self->var_step) == 2) {
      PyArrayObject *step_array;
      self->block_size = step_shape[1];
      if (step_shape[0] == 0)
	self->steps = 0;
      else {
	PyNetCDFIndex *indices;
	int i;
	self->steps = (step_shape[0]-1)*step_shape[1];
	indices = PyNetCDFVariable_Indices(self->var_step);
	indices[0].start = step_shape[0]-1;
	indices[0].stop = step_shape[0];
	indices[0].item = 1;
	step_array = PyNetCDFVariable_ReadAsArray(self->var_step, indices);
	if (step_array == NULL)
	  return NULL;
	for (i = 0; i < step_shape[1]; i++) {
	  if (((long *)step_array->data)[i] == NC_FILL_INT)
	    break;
	  self->steps++;
	}
	Py_DECREF(step_array);
      }
    }
    else {
      self->block_size = 1;
      self->steps = step_shape[0];
    }
    if (mode[0] == 'a') {
      PyObject *name;
      PyNetCDFVariableObject *variable;
      Py_ssize_t pos = 0;
      while (PyDict_Next(self->file->variables, &pos, &name,
			 (PyObject **)&variable))
	if (variable->type == PyArray_FLOAT
	    || variable->type == PyArray_DOUBLE) {
	  self->floattype = variable->type;
	  break;
	}
      self->first_step = 1;
      if (self->cycle > 0) {
	PyObject *last = PyNetCDFFile_GetAttribute(self->file, "last_step");
	if (last != NULL) {
	  int step = *(int *)((PyArrayObject *)last)->data;
	  if (step >= 0)
	    self->steps = (step+1)%self->cycle;
	}
      }
    }
    n_ob = PyDict_GetItemString(self->file->dimensions, atom_number);
    if (n_ob == NULL)
      goto error;
    self->trajectory_atoms = PyInt_AsLong(n_ob);
    description_var = PyNetCDFFile_GetVariable(self->file, "description");
    if (description_var != NULL) {
      PyStringObject *traj_description =
	           PyNetCDFVariable_ReadAsString(description_var);
      if (traj_description == NULL)
	goto error;
      if (!verify_description(PyString_AsString(description),
			      PyString_AsString((PyObject *)
						    traj_description))) {
	PyErr_SetString(PyExc_ValueError,
			"trajectory file not compatible with universe");
	goto error;
      }
      Py_DECREF(traj_description);
    }
    else
      PyErr_Clear();
  }
  self->last_flush = clock();
  return self;

error:
  trajectory_dealloc(self);
  return NULL;
}

static PyObject *
Trajectory(PyObject *self, PyObject *args)
{
  PyObject *universe, *description, *index_map;
  char *filename;
  char *mode = "r";
  int dpflag = 0;
  int cycle = 0;
  int block_size = 1;

  if (!PyArg_ParseTuple(args, "OO!Os|siii:Trajectory",
			&universe, &PyString_Type, &description, &index_map, 
			&filename, &mode, &dpflag, &cycle, &block_size))
    return NULL;
  if (index_map == Py_None)
    index_map = NULL;
  else if (!PyArray_Check(index_map)) {
    PyErr_SetString(PyExc_TypeError, "index map must be an array");
    return NULL;
  }
  return (PyObject *)PyTrajectory_Open(universe, description,
				       (PyArrayObject *)index_map,
				       filename, mode,
				       dpflag?PyArray_DOUBLE:PyArray_FLOAT,
				       cycle, block_size);
}

/*
 * Handle regular output during trajectory generation.
 */

enum PySpec_TYPE {PySpec_None, PySpec_Trajectory, PySpec_Print,
		  PySpec_Function};

/* Preprocess an output specification */

static int
get_spec(PyObject *universe, PyObject *spec,
	 PyTrajectoryOutputSpec *output, int type, char *description,
	 PyTrajectoryVariable *data, int nvar)
{
  static char text[200];
  int i;

  output->first = PyInt_AsLong(PyTuple_GetItem(spec, (Py_ssize_t)1));
  output->last = PyInt_AsLong(PyTuple_GetItem(spec, (Py_ssize_t)2));
  output->frequency = PyInt_AsLong(PyTuple_GetItem(spec, (Py_ssize_t)3));
  output->close = 0;
  output->type = type;
  output->destination = NULL;
  output->parameters = NULL;
  output->scratch = NULL;

  if (type != PySpec_Function) {

    PyObject *what = PyTuple_GetItem(spec, (Py_ssize_t)5);
    Py_ssize_t n = PyObject_Length(what);
    output->destination = PyTuple_GetItem(spec, (Py_ssize_t)4);
    if (output->destination == Py_None)
      return 0;
    output->what = 0;
    while (n-- > 0) {
      PyObject *item = PyObject_GetItem(what, PyInt_FromSsize_t(n));
      PyTrajectory_DataClassName *class = class_names;
      char *s;
      if (!PyString_Check(item)) {
	PyErr_SetString(PyExc_TypeError, "output item not a string");
	Py_DECREF(item);
	return -1;
      }
      s = PyString_AsString(item);
      while (class->name != NULL) {
	if (strcmp(s, class->name) == 0)
	  output->what |= class->number;
	class++;
      }
      Py_DECREF(item);
    }

  }

  if (type == PySpec_Trajectory) {

    PyTrajectoryObject *trajectory;

    output->variables = (PyObject **)malloc(nvar*sizeof(PyObject *));
    if (output->variables == NULL)
      return -1;
    trajectory = (PyTrajectoryObject *)output->destination;
    while (!PyTrajectory_Check(output->destination)) {
      output->destination = PyObject_GetAttrString(output->destination,
						   "trajectory");
      trajectory = (PyTrajectoryObject *)output->destination;
      if (output->destination == NULL) {
	PyErr_SetString(PyExc_TypeError, "not a trajectory");
	return -1;
      }
    }
    Py_INCREF(output->destination);
    if (!PyTrajectory_Check(output->destination)) {
      PyErr_SetString(PyExc_TypeError, "not a trajectory");
      Py_DECREF(output->destination);
      return -1;
    }
    if (description != NULL) {
      if (strlen(description) < 150)
	strcpy(text, description);
      else
	strcpy(text, "Trajectory");
      strcat(text, " started %s");
      if (PyTrajectory_TimeStamp(trajectory, text) == -1)
	return -1;
    }
    for (i = 0; i < nvar; i++) {
      if (output->what & data[i].class) {
	output->variables[i] =
	  PyTrajectory_GetVariable(trajectory, data[i].name,
				   data[i].type, 0, data[i].unit, 1);
	if (output->variables[i] == NULL) {
	  Py_DECREF(output->destination);
	  free(output->variables);
	  return -1;
	}
      }
    }

  }

  if (type == PySpec_Print) {

    if (PyString_Check(output->destination)) {
      PyObject *file =
	PyFile_FromString(PyString_AsString(output->destination), "a");
      if (file == NULL)
	return -1;
      output->destination = file;
      output->close = 1;
    }
    else
      Py_INCREF(output->destination);
    if (!PyObject_HasAttrString(output->destination, "write")) {
      PyErr_SetString(PyExc_TypeError, "not a file");
      Py_DECREF(output->destination);
      return -1;
    }

  }

  if (type == PySpec_Function) {

    output->function = (trajectory_fn *)
                       PyCObject_AsVoidPtr(PyTuple_GetItem(spec,
							   (Py_ssize_t)4));
    output->parameters = PyTuple_GetItem(spec, (Py_ssize_t)5);
    Py_INCREF(output->parameters);
    if (!output->function(data, output->parameters, -1, &output->scratch))
      return -1;

  }

  return 1;
}

static PyTrajectoryOutputSpec *
PyTrajectory_OutputSpecification(PyObject *universe,
				 PyListObject *spec_list, char *description,
				 PyTrajectoryVariable *data)
{
  PyTrajectoryOutputSpec *output;
  PyTrajectoryVariable *v;
  Py_ssize_t nspecs = PyList_Size((PyObject *)spec_list);
  int nvar = 0;
  for (v = data; v->name != NULL; v++)
    nvar++;
  output = (PyTrajectoryOutputSpec *)
    malloc((nspecs+1)*sizeof(PyTrajectoryOutputSpec));
  if (output != NULL) {
    int n = 0;
    int ret;
    int i;
    for (i = 0; i < nspecs; i++) {
      PyObject *spec = PyList_GetItem((PyObject *)spec_list, (Py_ssize_t)i);
      PyObject *which;
      char *which_str;
      int type;
      if (!PyTuple_Check(spec)) {
	PyErr_SetString(PyExc_TypeError, "must be a tuple");
	free(output);
	return NULL;
      }
      which = PyTuple_GetItem(spec, (Py_ssize_t)0);
      if (!PyString_Check(which)) {
	PyErr_SetString(PyExc_TypeError, "must be a string");
	free(output);
	return NULL;
      }
      which_str = PyString_AsString(which);
      if (strcmp(which_str, "print") == 0)
	type = PySpec_Print;
      else if (strcmp(which_str, "trajectory") == 0)
	type = PySpec_Trajectory;
      else if (strcmp(which_str, "function") == 0)
	type = PySpec_Function;
      else {
	PyErr_SetString(PyExc_TypeError, "illegal specification id");
	free(output);
	return NULL;
      }
      ret = get_spec(universe, spec, output+n, type, description, data, nvar);
      if (ret == -1)
	return NULL;
      if (ret == 1)
	n++;
    }
    output[n].type = PySpec_None;
  }
  return output;
}

/* Clean up */
static void
PyTrajectory_OutputFinish(PyTrajectoryOutputSpec *spec, int step,
			  int error_flag, int time_stamp_flag,
			  PyTrajectoryVariable *data)
{
  PyTrajectoryOutputSpec *s = spec;
  PyTrajectory_Output(spec, -step, data, NULL);
  while (spec->type != PySpec_None) {
    if (spec->type == PySpec_Trajectory) {
      char *text;
      PyTrajectory_Flush((PyTrajectoryObject *)spec->destination);
      if (error_flag) {
	if (PyErr_CheckSignals())
	  text = "Trajectory interrupted %s";
	else
	  text = "Trajectory terminated by error %s";
      }
      else
	text = "Trajectory finished %s";
      if (time_stamp_flag || error_flag)
	PyTrajectory_TimeStamp((PyTrajectoryObject *)spec->destination, text);
      PyTrajectory_Flush((PyTrajectoryObject *)spec->destination);
      free(spec->variables);
    }
    if (spec->type == PySpec_Function) {
      spec->function(data, spec->parameters, -2, &spec->scratch);
    }
    if (spec->close) {
      if (spec->type == PySpec_Trajectory)
	PyTrajectory_Close((PyTrajectoryObject *)spec->destination);
      else
	PyObject_CallMethod(spec->destination, "close", NULL);
    }
    Py_XDECREF(spec->destination);
    Py_XDECREF(spec->parameters);
    spec++;
  }
  free(s);
}

/* Do output for one step */
static int
PyTrajectory_Output(PyTrajectoryOutputSpec *spec, int step,
		    PyTrajectoryVariable *data, PyThreadState **thread)
{
  PyTrajectoryVariable *var;
  double *array_data;
  char buffer[100];
  int interrupt = 0;
  int i;

  for (var = data; var->name != NULL; var++)
    var->modified = 0;
  while (spec->type != PySpec_None) {
    if ((step >= spec->first && step < spec->last &&
	 (step-spec->first) % spec->frequency == 0) || step < 0) {
      if (spec->type == PySpec_Print && step >= 0) {
	if (thread != NULL)
	  PyEval_RestoreThread(*thread);
	if (step > 0)
	  PyFile_WriteString("\n", spec->destination);
	sprintf(buffer, "Step %d\n", step);
	PyFile_WriteString(buffer, spec->destination);

	for (var = data; var->name != NULL; var++)
	  if (spec->what & var->class) {
	    switch (var->type) {
	    case PyTrajectory_Scalar:
	      sprintf(buffer, var->text, *var->value.dp);
	      PyFile_WriteString(buffer, spec->destination);
	      break;
	    case PyTrajectory_ParticleScalar:
	      sprintf(buffer, var->text);
	      PyFile_WriteString(buffer, spec->destination);
	      array_data = (double *)var->value.array->data;
	      for (i = 0; i < var->value.array->dimensions[0]; i++) {
		sprintf(buffer, "  %5d: %g\n", i, *array_data++);
		PyFile_WriteString(buffer, spec->destination);
	      }
	      break;
	    case PyTrajectory_ParticleVector:
	      sprintf(buffer, var->text);
	      PyFile_WriteString(buffer, spec->destination);
	      array_data = (double *)var->value.array->data;
	      for (i = 0; i < var->value.array->dimensions[0]; i++) {
		sprintf(buffer, "  %5d: %g, %g, %g\n", i,
			array_data[0], array_data[1], array_data[2]);
		PyFile_WriteString(buffer, spec->destination);
		array_data += 3;
	      }
	      break;
	    case PyTrajectory_BoxSize:
	      sprintf(buffer, var->text);
	      PyFile_WriteString(buffer, spec->destination);
	      array_data = var->value.dp;
	      for (i = 0; i < var->length; i++) {
		sprintf(buffer, " %g", *array_data++);
		PyFile_WriteString(buffer, spec->destination);
		if (i == var->length-1)
		  sprintf(buffer, "\n");
		else
		  sprintf(buffer, ",");
		PyFile_WriteString(buffer, spec->destination);
	      }
	      break;
	    }
	  }

	if (PyObject_HasAttrString(spec->destination, "flush"))
	  PyObject_CallMethod(spec->destination, "flush", NULL);
	if (PyErr_CheckSignals())
	  interrupt = -1;
	if (thread != NULL)
	  *thread = PyEval_SaveThread();
      }
      if (spec->type == PySpec_Trajectory) {
	PyTrajectoryObject *trajectory =
	  (PyTrajectoryObject *)spec->destination;
	PyObject **tvar;
	clock_t cpu_time;
	if (trajectory->cycle > 0 && step < 0)
	  step = -step;
	if (step >= 0) {
	  if (thread != NULL)
	    PyEval_RestoreThread(*thread);
	  if (PyTrajectory_Step(trajectory, step) == -1)
	    return -1;
	  for (var = data, tvar = spec->variables; var->name != NULL;
	       var++, tvar++)
	    if (spec->what & var->class) {
	      PyArrayObject *array;
	      switch (var->type) {
	      case PyTrajectory_Scalar:
		if (PyTrajectory_WriteFloats(trajectory, *tvar,
					     var->value.dp, 1) == -1)
		  return -1;
		break;
	      case PyTrajectory_ParticleScalar:
		if (trajectory->index_map != NULL) {
		  int i;
		  long *indices = (long *)trajectory->index_map->data;
		  double *source = (double *)var->value.array->data;
		  if (trajectory->floattype == PyArray_DOUBLE) {
		    double *dest = (double *)trajectory->sbuffer->data;
		    for (i = 0; i < trajectory->index_map->dimensions[0]; i++)
		      dest[i] = source[indices[i]];
		  }
		  else {
		    float *dest = (float *)trajectory->sbuffer->data;
		    for (i = 0; i < trajectory->index_map->dimensions[0]; i++)
		      dest[i] = (float)source[indices[i]];
		  }
		  array = trajectory->sbuffer;
		}
		else {
		  if (trajectory->floattype == PyArray_DOUBLE)
		    array = var->value.array;
		  else {
		    double *source = (double *)var->value.array->data;
		    float *dest = (float *)trajectory->sbuffer->data;
		    int natoms = var->value.array->dimensions[0];
		    int i;
		    for (i = 0; i < natoms; i++)
		      dest[i] = (float)source[i];
		    array = trajectory->sbuffer;
		  }
		}
		if (PyTrajectory_WriteArray(trajectory, *tvar, array) == -1)
		  return -1;
		break;
	      case PyTrajectory_ParticleVector:
		if (trajectory->index_map != NULL) {
		  int i;
		  long *indices = (long *)trajectory->index_map->data;
		  vector3 *source = (vector3 *)var->value.array->data;
		  if (trajectory->floattype == PyArray_DOUBLE) {
		    vector3 *dest = (vector3 *)trajectory->vbuffer->data;
		    for (i = 0; i < trajectory->index_map->dimensions[0]; i++){
		      dest[i][0] = source[indices[i]][0];
		      dest[i][1] = source[indices[i]][1];
		      dest[i][2] = source[indices[i]][2];
		    }
		  }
		  else {
		    float *dest = (float *)trajectory->vbuffer->data;
		    for (i = 0; i < trajectory->index_map->dimensions[0]; i++){
		      dest[3*i] = (float)source[indices[i]][0];
		      dest[3*i+1] = (float)source[indices[i]][1];
		      dest[3*i+2] = (float)source[indices[i]][2];
		    }
		  }
		  array = trajectory->vbuffer;
		}
		else {
		  if (trajectory->floattype == PyArray_DOUBLE)
		    array = var->value.array;
		  else {
		    vector3 *source = (vector3 *)var->value.array->data;
		    int natoms = var->value.array->dimensions[0];
		    float *dest = (float *)trajectory->vbuffer->data;
		    int i;
		    for (i = 0; i < natoms; i++) {
		      dest[3*i] = (float)source[i][0];
		      dest[3*i+1] = (float)source[i][1];
		      dest[3*i+2] = (float)source[i][2];
		    }
		    array = trajectory->vbuffer;
		  }
		}
		if (PyTrajectory_WriteArray(trajectory, *tvar, array) == -1)
		  return -1;
		break;
	      case PyTrajectory_BoxSize:
		if (PyTrajectory_WriteFloats(trajectory, *tvar,
					     var->value.dp, var->length) == -1)
		  return -1;
		break;
	      }
	    }
	  cpu_time = clock();
	  if (trajectory->cycle > 0) {
	    if (PyTrajectory_SetAttribute(trajectory, "last_step",
					  PyInt_FromLong(trajectory->steps-1))
		== -1)
	      return -1;
	    if (PyTrajectory_Flush(trajectory) == -1)
	      return -1;
	  }
	  else if (cpu_time-trajectory->last_flush > 900*CLOCKS_PER_SEC) {
	    if (PyTrajectory_Flush(trajectory) == -1)
	      return -1;
	    trajectory->last_flush = cpu_time;
	  }
	  if (PyErr_CheckSignals())
	    interrupt = -1;
	  if (thread != NULL)
	    *thread = PyEval_SaveThread();
	}
      }
      if (spec->type == PySpec_Function) {
	if (!spec->function(data, spec->parameters, step, &spec->scratch))
	  return -1;
      }
    }
    spec++;
  }
  return interrupt;
}


/* Transformation from/to box coordinates */

static PyObject *
boxTransformation(PyObject *dummy, PyObject *args)
{
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *pt_in, *pt_out, *box_size;
  int to_box;
  vector3 *in, *out;
  double *box;

  if (!PyArg_ParseTuple(args, "O!O!O!O!i",
			&PyUniverseSpec_Type, &universe_spec,
			&PyArray_Type, &pt_in,
			&PyArray_Type, &pt_out,
			&PyArray_Type, &box_size,
			&to_box))
    return NULL;
  in = (vector3 *)pt_in->data;
  out = (vector3 *)pt_out->data;
  box = (double *)box_size->data;
  universe_spec->trajectory_function(in, out, pt_in->dimensions[0],
				     box, to_box);
  Py_INCREF(Py_None);
  return Py_None;
}

/* Snapshot generator: process the current situation as one step of
   a trajectory */

static PyObject *
snapshot(PyObject *dummy, PyObject *args)
{
  PyObject *universe;
  PyListObject *spec_list;
  PyObject *data_dict, *data_key, *data_value;
  PyTrajectoryVariable *vars, *v, *vtemp;
  PyTrajectoryOutputSpec *output;
  char string_buffer[80];
  char *name;
  int energy_terms;
  Py_ssize_t dict_pos;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "OO!O!i", &universe,
			&PyDict_Type, &data_dict,
			&PyList_Type, &spec_list,
			&energy_terms))
    return NULL;
  
  vars = (PyTrajectoryVariable *)malloc((9+energy_terms)
					* sizeof(PyTrajectoryVariable));
  if (vars == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  vars[0].name = "temperature";
  vars[0].text = "Temperature: %lf\n";
  vars[0].unit = temperature_unit_name;
  vars[0].type = PyTrajectory_Scalar;
  vars[0].class = PyTrajectory_Thermodynamic;
  vars[0].value.dp = NULL;
  vars[1].name = "pressure";
  vars[1].text = "Pressure: %lf\n";
  vars[1].unit = pressure_unit_name;
  vars[1].type = PyTrajectory_Scalar;
  vars[1].class = PyTrajectory_Thermodynamic;
  vars[1].value.dp = NULL;
  vars[2].name = "configuration";
  vars[2].text = "Configuration:\n";
  vars[2].unit = length_unit_name;
  vars[2].type = PyTrajectory_ParticleVector;
  vars[2].class = PyTrajectory_Configuration;
  vars[2].value.array = NULL;
  vars[3].name = "velocities";
  vars[3].text = "Velocities:\n";
  vars[3].unit = velocity_unit_name;
  vars[3].type = PyTrajectory_ParticleVector;
  vars[3].class = PyTrajectory_Velocities;
  vars[3].value.array = NULL;
  vars[4].name = "gradients";
  vars[4].text = "Energy gradients:\n";
  vars[4].unit = energy_gradient_unit_name;
  vars[4].type = PyTrajectory_ParticleVector;
  vars[4].class = PyTrajectory_Gradients;
  vars[4].value.array = NULL;
  vars[5].name = "gradient_norm";
  vars[5].text = "Gradient norm: %lf\n";
  vars[5].unit = energy_gradient_unit_name;
  vars[5].type = PyTrajectory_Scalar;
  vars[5].class = PyTrajectory_Energy;
  vars[5].value.dp = NULL;
  vars[6].name = "box_size";
  vars[6].text = "Box size:";
  vars[6].unit = length_unit_name;
  vars[6].type = PyTrajectory_BoxSize;
  vars[6].class = PyTrajectory_Configuration;
  vars[6].value.dp = NULL;
  vars[7].name = "time";
  vars[7].text = "Time: %lf\n";
  vars[7].unit = time_unit_name;
  vars[7].type = PyTrajectory_Scalar;
  vars[7].class = PyTrajectory_Time;
  vars[7].value.dp = NULL;
  vars[8].name = NULL;

  for (v = vars; v->name != NULL; ) {
    PyObject *value = PyDict_GetItemString(data_dict, v->name);
    if (value == NULL) {
      for (vtemp = v; vtemp->name != NULL; vtemp++)
	*vtemp = *(vtemp+1);
    }
    else {
      if (v->type == PyTrajectory_Scalar) {
	v->value.dp = (double *)malloc(sizeof(double));
	if (v->value.dp == NULL) {
	  PyErr_NoMemory();
	  goto error2;
	}
	*(v->value.dp) = PyFloat_AsDouble(value);
      }
      else if (v->type == PyTrajectory_BoxSize) {
	v->value.dp = (double *)((PyArrayObject *)value)->data;
	v->length = ((PyArrayObject *)value)->dimensions[0];
      }
      else
	v->value.array = (PyArrayObject *)value;
      v++;
    }
  }

  dict_pos = 0;
  while (PyDict_Next(data_dict, &dict_pos, &data_key, &data_value)) {
    name = PyString_AsString(data_key);
    if (strcmp(name+strlen(name)-7, "_energy") == 0) {
      char *s;
      strcpy(string_buffer, name);
      for (s = string_buffer; *s; s++)
	if (*s == '_')
	  *s = ' ';
      strcpy(string_buffer + strlen(string_buffer), ": %lf\n");
      v->name = name;
      v->text = string_buffer;
      v->unit = energy_unit_name;
      v->type = PyTrajectory_Scalar;
      v->class = PyTrajectory_Energy;
      v->value.dp = (double *)malloc(sizeof(double));
      if (v->value.dp == NULL) {
	PyErr_NoMemory();
	goto error2;
      }
      *(v->value.dp) = PyFloat_AsDouble(data_value);
      v++;
    }
  }
  v->name = NULL;

  output = PyTrajectory_OutputSpecification(universe, spec_list, NULL, vars);
  if (output == NULL)
    goto error2;
  if (PyTrajectory_Output(output, 1, vars, NULL) == -1)
    goto error;

  PyTrajectory_OutputFinish(output, 1, 0, 0, vars);
  for (v = vars; v->name != NULL; v++)
    if (v->type == PyTrajectory_Scalar)
      free(v->value.dp);
  free(vars);
  Py_INCREF(Py_None);
  return Py_None;

error:
  PyTrajectory_OutputFinish(output, 1, 1, 0, vars);
error2:
  for (v = vars; v->name != NULL; v++)
    if (v->type == PyTrajectory_Scalar)
      free(v->value.dp);
  free(vars);
  return NULL;
}

/* Trajectory reader */

static PyObject *
readTrajectory(PyObject *dummy, PyObject *args)
{
  PyObject *universe;
  PyTrajectoryObject *input;
  PyListObject *spec_list;

  PyObject *input_vars;
  Py_ssize_t ninvars;
  PyTrajectoryVariable *vars;
  PyObject *name;
  PyNetCDFVariableObject *var;
  Py_ssize_t i;
  int j;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "OO!O!", &universe,
			&PyTrajectory_Type, &input,
			&PyList_Type, &spec_list))
    return NULL;

  /* Create temporary arrays for all variables */
  input_vars = input->file->variables;
  ninvars = PyDict_Size(input_vars);
  vars = (PyTrajectoryVariable *)
              malloc((ninvars+1)*sizeof(PyTrajectoryVariable));
  if (vars == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  i = 0; j = 0;
  while (PyDict_Next(input_vars, &i, &name, (PyObject **)&var)) {
    char *n = PyString_AsString(name);
    if (var->unlimited && strcmp(n, "step") != 0) {
      if (var->nd == 3) { /* particle array */
	int shape[2];
	vars[j].type = PyTrajectory_ParticleVector;
	shape[0] = input->natoms;
	shape[1] = 3;
	vars[j].value.array = (PyArrayObject *)
	                           PyArray_FromDims(2, shape,input->floattype);
	if (vars[j].value.array == NULL) {
	  PyErr_NoMemory();
	  goto error;
	}
      }
      else if (var->nd == 2) { /* particle scalar */
	/* not implemented */
	continue;
      }
      else { /* scalar */
	vars[j].type = PyTrajectory_Scalar;
	vars[j].value.dp = (double *)malloc(sizeof(double));
	if (vars[j].value.dp == NULL) {
	  PyErr_NoMemory();
	  goto error;
	}
      }
      vars[j].name = n;
      vars[j].unit = PyString_AsString(PyNetCDFVariable_GetAttribute(var,
								     "units"));
      vars[j].text = "";
      vars[j].class = 0;
      j += 1;
    }
  }
  vars[j].name = NULL;

  /* Loop over all steps */

  /* Error exit */
error:
  return NULL;
}


/* Table of functions defined in the module */

static PyMethodDef module_methods[] = {
  {"Trajectory", Trajectory, 1},
  {"boxTransformation", boxTransformation, 1},
  {"snapshot", snapshot, 1},
  {"readTrajectory", readTrajectory, 1},
  {NULL, NULL}		/* sentinel */
};


/* Module initialization */

DL_EXPORT(void)
initMMTK_trajectory(void)
{
  PyObject *module, *dict;
  PyObject *netcdf;
  static void *PyTrajectory_API[PyTrajectory_API_pointers];

  /* Patch object type */
#ifdef EXTENDED_TYPES
  if (PyType_Ready(&PyTrajectory_Type) < 0)
    return;
#else
  PyTrajectory_Type.ob_type = &PyType_Type;
#endif

  /* Create the module and add the type object */
  module = Py_InitModule("MMTK_trajectory", module_methods);
  dict = PyModule_GetDict(module);
  PyDict_SetItemString(dict, "trajectory_type",(PyObject *)&PyTrajectory_Type);

  /* Import the array module */
#ifdef import_array
  import_array();
#endif

  /* Import the universe module */
  module = PyImport_ImportModule("MMTK_universe");
  if (module != NULL) {
    PyObject *module_dict = PyModule_GetDict(module);
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API");
    if (PyCObject_Check(c_api_object))
      PyUniverse_API = (void **)PyCObject_AsVoidPtr(c_api_object);
  }

  /* Initialize C API pointer array and store in module */
  PyTrajectory_API[PyTrajectory_Type_NUM] = (void *)&PyTrajectory_Type;
  PyTrajectory_API[PyTrajectory_Open_NUM] = (void *)&PyTrajectory_Open;
  PyTrajectory_API[PyTrajectory_Close_NUM] = (void *)&PyTrajectory_Close;
  PyTrajectory_API[PyTrajectory_OutputSpecification_NUM] =
    (void *)&PyTrajectory_OutputSpecification;
  PyTrajectory_API[PyTrajectory_OutputFinish_NUM] =
    (void *)&PyTrajectory_OutputFinish;
  PyTrajectory_API[PyTrajectory_Output_NUM] = (void *)&PyTrajectory_Output;
  PyDict_SetItemString(dict, "_C_API",
		       PyCObject_FromVoidPtr((void *)PyTrajectory_API, NULL));

  /* Define maxint */
  PyDict_SetItemString(dict, "maxint", PyInt_FromLong(INT_MAX));

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_trajectory");

  /* Import netcdf and retrieve its C API address array */
  netcdf = PyImport_ImportModule("Scientific.IO.NetCDF");
  if (netcdf != NULL) {
    PyObject *module_dict = PyModule_GetDict(netcdf);
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API");
    fflush(stdout);
    if (PyCObject_Check(c_api_object)) {
      PyNetCDF_API = (void **)PyCObject_AsVoidPtr(c_api_object);
    }
  }
  else
    PyErr_Clear();
}
