/* I/O in DCD format
 *
 * Written by Lutz Ehrlich
 * Adapted to MMTK conventions by Konrad Hinsen
 */

#include "MMTK/universe.h"
#include "MMTK/trajectory.h"
#include "MMTK/readdcd.h"

/* Global variables */

double angstrom_factor;
double akma_time_factor;

/* Allocate and initialize Output variable descriptors */

static  PyTrajectoryVariable *
get_data_descriptors(PyArrayObject *configuration, double *time,
		     double *box_size, int box_size_length)
{
  static PyTrajectoryVariable vars[4];
  if (vars != NULL) {
    vars[0].name = "time";
    vars[0].text = "Time: %lf\n";
    vars[0].unit = time_unit_name;
    vars[0].type = PyTrajectory_Scalar;
    vars[0].class = PyTrajectory_Time;
    vars[0].value.dp = time;
    vars[1].name = "configuration";
    vars[1].text = "Configuration:\n";
    vars[1].unit = length_unit_name;
    vars[1].type = PyTrajectory_ParticleVector;
    vars[1].class = PyTrajectory_Configuration;
    vars[1].value.array = configuration;
    if (box_size_length > 0) {
      vars[2].name = "box_size";
      vars[2].text = "Box size:";
      vars[2].unit = length_unit_name;
      vars[2].type = PyTrajectory_BoxSize;
      vars[2].class = PyTrajectory_Configuration;
      vars[2].value.dp = box_size;
      vars[2].length = box_size_length;
      vars[3].name = NULL;
    }
    else
      vars[2].name = NULL;
  }
  return vars;
}


static  PyObject *
readDCD(PyObject *dummy, PyObject *args)
{
  PyObject *universe;
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *configuration;
  PyListObject *spec_list;
  PyTrajectoryOutputSpec *output;
  vector3 *x;
  int atoms;

  char buffer[100];
  PyTrajectoryVariable *data_descriptors;
  int j;

  /* DCD transformation variables */
  char * dcdFileName;
  FILE * dcdFile;
  int *  dcdFreeatoms;
  int    dcdAtoms;
  int    dcdNFrames=0;
  int    dcdFrameStart=0;
  int    dcdFrameSkip=0;
  float  dcdTimeStep;
  int    dcdNamnf;
  float *dcdX = NULL;
  float *dcdY = NULL;
  float *dcdZ = NULL;
  int    dcdErrcode;
  int    currFrames=0;
  double time;

  dcdFreeatoms = NULL;
    

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "OO!O!s", &universe,
			&PyArray_Type, &configuration,
			&PyList_Type, &spec_list, &dcdFileName))
    return NULL;

  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;
  atoms = configuration->dimensions[0];
  x = (vector3 *)configuration->data;

  /* Prepare output data descriptors */
  data_descriptors =
    get_data_descriptors(configuration, &time,
			 universe_spec->geometry_data,
			 universe_spec->geometry_data_length);

  /* Initialize output */
  output = PyTrajectory_OutputSpecification(universe, spec_list,
					    dcdFileName,
					    data_descriptors);
  if (output == NULL)
    return NULL;

  /* open DCD file  */
  dcdFile = open_dcd_read( dcdFileName );
  if ( dcdFile == NULL ){
    PyErr_SetString(PyExc_IOError, "Cannot open file");
    goto error;
  }

  /* read header */
  dcdErrcode = read_dcdheader(dcdFile, &dcdAtoms, &dcdNFrames, 
			      &dcdFrameStart,  &dcdFrameSkip, &dcdTimeStep,
			      &dcdNamnf, &dcdFreeatoms);
  if ( dcdErrcode == DCD_BADFORMAT ) {
    PyErr_SetString(PyExc_IOError, "Not a DCD file");
    goto error;
  }
  else if ( dcdErrcode != 0 ) {
    PyErr_SetString(PyExc_IOError, "DCD reading error");
    goto error;
  }
  if( atoms != dcdAtoms ){
    snprintf(buffer, sizeof(buffer),
             "number of atoms in DCD file (%d) doesn't "
	    "match universe (%d)", dcdAtoms, atoms);
    PyErr_SetString(PyExc_ValueError, buffer);
    goto error;
  }
  if ( dcdNamnf != 0 ){
    PyErr_SetString(PyExc_ValueError, "Can't read DCD files with free atoms");
    goto error;
  }

  /* allocate the dcd{X,Y,Z} arrays */
  dcdX = ( float*) malloc(dcdAtoms * sizeof(float) );
  dcdY = ( float*) malloc(dcdAtoms * sizeof(float) );
  dcdZ = ( float*) malloc(dcdAtoms * sizeof(float) );
  if( (dcdX==NULL) || (dcdY==NULL) || (dcdZ==NULL) ){
    PyErr_NoMemory();
    goto error;
  }
  
  /* read in the frames one after the other */
  currFrames = 0;
  time = 0.;
  while ( 1 ){
    int err_code = read_dcdstep(dcdFile, dcdAtoms, dcdX, dcdY, dcdZ,
				dcdNamnf, (currFrames == 0), dcdFreeatoms);
    if (err_code == -1)
      break;
    if (err_code < 0) {
      PyErr_SetString(PyExc_IOError, "DCD read error");
      goto error;
    }
    for (j = 0; j < dcdAtoms; j++) {
       x[j][0] = angstrom_factor * dcdX[j];
       x[j][1] = angstrom_factor * dcdY[j];
       x[j][2] = angstrom_factor * dcdZ[j];
    }
    if (PyTrajectory_Output(output, currFrames, data_descriptors, NULL) == -1)
      goto error;
    currFrames++;
    time += dcdFrameSkip*dcdTimeStep*akma_time_factor;
  }
  close_dcd_read(dcdFile,0,dcdFreeatoms);

  /* Clean up and return None */
  if( dcdX != NULL )
    free(dcdX);
  if( dcdY != NULL )
    free(dcdY);
  if( dcdZ != NULL )
    free(dcdZ);
  PyTrajectory_OutputFinish(output, currFrames-1, 0, 1, data_descriptors);
  Py_INCREF(Py_None);
  return Py_None;

  /* Clean up and return error */
error:
  if( dcdX != NULL )
    free(dcdX);
  if( dcdY != NULL )
    free(dcdY);
  if( dcdZ != NULL )
    free(dcdZ);
  close_dcd_read(dcdFile,0,dcdFreeatoms);
  PyTrajectory_OutputFinish(output, currFrames, 1, 1, data_descriptors);
  return NULL;
}


static PyObject *
writeOpenDCD(PyObject *dummy, PyObject *args)
{
  /* DCD transformation variables */
  char * dcdFileName;
  int    dcdAtoms;
  int    dcdNFrames=0;
  int    dcdFrameStart=0;
  int    dcdNSavc;
  double time;
  FILE * fd;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "siiiid", &dcdFileName,&dcdAtoms, &dcdNFrames, 
			&dcdFrameStart, &dcdNSavc, &time))
    return NULL;

  /* open DCD file  */
  fd = open_dcd_write( dcdFileName );
  if ( fd == NULL ){
    PyErr_SetString(PyExc_IOError, "Cannot open file");
    return NULL;
  }

  /* write header */
  write_dcdheader(fd , dcdFileName, dcdAtoms, dcdNFrames, 
                  dcdFrameStart,  dcdNSavc,
                  time/akma_time_factor);
  return (PyObject*)PyCObject_FromVoidPtr((void *) fd, NULL);
}

static  PyObject *
writeDCDStep(PyObject *dummy, PyObject *args)
{
  int    dcdAtoms;
  float *dcdX = NULL;
  float *dcdY = NULL;
  float *dcdZ = NULL;
  int err;
  FILE * fd;
  PyObject *fd_cobj;
  PyArrayObject *xconfig, *yconfig, *zconfig;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "O!O!O!O!",
                        &PyCObject_Type, &fd_cobj,
			&PyArray_Type, &xconfig,
			&PyArray_Type, &yconfig,
			&PyArray_Type, &zconfig))
    return NULL;

  fd = PyCObject_AsVoidPtr(fd_cobj);
  dcdAtoms = xconfig->dimensions[0];
  dcdX = (float *)xconfig->data;
  dcdY = (float *)yconfig->data;
  dcdZ = (float *)zconfig->data;

  err = write_dcdstep( fd, dcdAtoms, dcdX, dcdY, dcdZ);
  if (err != 1) {
    PyErr_SetString(PyExc_IOError, "Couldn't write DCD step");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static  PyObject *
writeCloseDCD(PyObject *dummy, PyObject *args)
{
  PyObject *fd_cobj;
  FILE * fd;

  /* Parse and check arguments */
  if (!PyArg_ParseTuple(args, "O!",
                        &PyCObject_Type, &fd_cobj))
    return NULL;

  fd = PyCObject_AsVoidPtr(fd_cobj);
  close_dcd_read(fd, 0, NULL);

  Py_INCREF(Py_None);
  return Py_None;

}


static PyMethodDef DCD_methods[] = {
  {"readDCD", readDCD, 1},
  {"writeOpenDCD", writeOpenDCD, 1},
  {"writeDCDStep", writeDCDStep, 1},
  {"writeCloseDCD", writeCloseDCD, 1},
  {NULL, NULL}		/* sentinel */
};


/* Initialization function for the module */

DL_EXPORT(void)
initMMTK_DCD(void)
{
  PyObject *units;

  /* Create the module and add the functions */
  Py_InitModule("MMTK_DCD", DCD_methods);

  /* Import the array module */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules */
  import_MMTK_trajectory();

  /* Get the length and time unit conversion factor from Units */
  units = PyImport_ImportModule("MMTK.Units");
  if (units != NULL) {
    PyObject *module_dict = PyModule_GetDict(units);
    PyObject *factor = PyDict_GetItemString(module_dict, "Ang");
    angstrom_factor = PyFloat_AsDouble(factor);
    factor = PyDict_GetItemString(module_dict, "akma_time");
    akma_time_factor = PyFloat_AsDouble(factor);
  }

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_DCD");
}
