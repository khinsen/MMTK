/* Molecular surface code
 *
 * Copyright 2000 by Peter McCluskey (pcm@rahul.net).
 * You may do anything you want with it, provided this notice is kept intact.
 */

#include "Python.h"
#include <math.h>
#include <time.h>
#include "MMTK/core.h"

typedef double TPoint[3];
typedef struct {
  int index;
  double dist2;
} NeighborData;

static void
normalize(TPoint p)
{
  double length = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  p[0] /= length;
  p[1] /= length;
  p[2] /= length;
}

static int
add_point(const TPoint p, TPoint *r, int r_index, PyObject *pt_dict)
{
  PyObject *p_orig;
  char key[512];
  snprintf(key, sizeof(key), "%.6f,%.6f,%.6f", p[0], p[1], p[2]);
  p_orig = PyDict_GetItemString(pt_dict, key);
  if(p_orig)
    Py_DECREF(p_orig);
  else
  {
    PyObject *py_one = PyInt_FromLong(1);
    PyDict_SetItemString(pt_dict, key, py_one);
    r[r_index][0] = p[0];
    r[r_index][1] = p[1];
    r[r_index][2] = p[2];
    ++r_index;
  }
  return r_index;
}

static int
tess_triangle(const TPoint v0, const TPoint v1, const TPoint v2, int npts,
	      TPoint *r, int r_index, PyObject *pt_dict)
{
    int n4 = npts/4;
    TPoint new_pt0;
    TPoint new_pt1;
    TPoint new_pt2;
    new_pt0[0] = v0[0] + v1[0];
    new_pt0[1] = v0[1] + v1[1];
    new_pt0[2] = v0[2] + v1[2];
    new_pt1[0] = v1[0] + v2[0];
    new_pt1[1] = v1[1] + v2[1];
    new_pt1[2] = v1[2] + v2[2];
    new_pt2[0] = v2[0] + v0[0];
    new_pt2[1] = v2[1] + v0[1];
    new_pt2[2] = v2[2] + v0[2];
    normalize(new_pt0);
    normalize(new_pt1);
    normalize(new_pt2);
    if(n4 <= 3)
    {
      /* other 2 points in pts will be created as pts[0] of */
      /* another triangle, except for top level points */
      r_index = add_point(v0, r, r_index, pt_dict);
      r_index = add_point(new_pt0, r, r_index, pt_dict);
      r_index = add_point(new_pt1, r, r_index, pt_dict);
      r_index = add_point(new_pt2, r, r_index, pt_dict);
      return r_index;
    }
    else
    {
      int i = tess_triangle(v0, new_pt0, new_pt2, n4, r, r_index, pt_dict);
      i = tess_triangle(new_pt0, v1, new_pt1, n4, r, i, pt_dict);
      i = tess_triangle(new_pt0, new_pt1, new_pt2, n4, r, i, pt_dict);
      return tess_triangle(new_pt1, v2, new_pt2, n4, r, i, pt_dict);
    }
}

static TPoint *
tesselate(int num_points)
{
    static const TPoint north = {0.0, 0.0, 1.0};
    static const TPoint south = {0.0, 0.0, -1.0};
    static const TPoint noon  = {1.0, 0.0, 0.0};
    static const TPoint night = {-1.0, 0.0, 0.0};
    static const TPoint dawn  = {0.0, 1.0, 0.0};
    static const TPoint dusk  = {0.0, -1.0, 0.0};
    int npts = (num_points - 2)/4;
    int i = 0;
    PyObject *pt_dict = PyDict_New();
    TPoint *r = (TPoint *)malloc(num_points * sizeof(TPoint));
    i = add_point(north, r, i, pt_dict);
    i = tess_triangle(north, dawn, noon, npts, r, i, pt_dict);
    i = add_point(noon,  r, i, pt_dict);
    i = tess_triangle(north, noon, dusk, npts, r, i, pt_dict);
    i = add_point(dusk,  r, i, pt_dict);
    i = tess_triangle(north, dusk, night, npts, r, i, pt_dict);
    i = add_point(night, r, i, pt_dict);
    i = tess_triangle(north, night, dawn, npts, r, i, pt_dict);
    i = add_point(dawn,  r, i, pt_dict);
    i = add_point(south, r, i, pt_dict);
    i = tess_triangle(south, dawn, night, npts, r, i, pt_dict);
    i = tess_triangle(south, night, dusk, npts, r, i, pt_dict);
    i = tess_triangle(south, dusk, noon, npts, r, i, pt_dict);
    i = tess_triangle(south, noon, dawn, npts, r, i, pt_dict);
    Py_DECREF(pt_dict);
    if(i != num_points)
    {
      free(r);
      return NULL;
    }
    return r;
}

static int
nbor_data_1_atom(PyObject *nbors, int i, PyObject *atom_data,
		 NeighborData *nbor_data)
{
  double max_dist_2;
  PyObject *boxes;
  PyObject *bsize;
  int j;
  double box_size;
  int n_nbors1 = 0;

  boxes = PyObject_GetAttrString(nbors, "boxes");
  bsize = PyObject_GetAttrString(nbors, "box_size");
  box_size = PyFloat_AsDouble(bsize);
  Py_DECREF(bsize);
  max_dist_2 = box_size*box_size;
  {
    PyObject *pos1 = PyList_GetItem(atom_data, (Py_ssize_t)i);
    double posx = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)0));
    double posy = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)1));
    double posz = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)2));
    int boxn;
    static const int nbor_boxes[27][3] = {
      { -1, -1, -1},
      { -1, -1,  0},
      { -1, -1,  1},
      { -1,  0, -1},
      { -1,  0,  0},
      { -1,  0,  1},
      { -1,  1, -1},
      { -1,  1,  0},
      { -1,  1,  1},
      {  0, -1, -1},
      {  0, -1,  0},
      {  0, -1,  1},
      {  0,  0, -1},
      {  0,  0,  0},
      {  0,  0,  1},
      {  0,  1, -1},
      {  0,  1,  0},
      {  0,  1,  1},
      {  1, -1, -1},
      {  1, -1,  0},
      {  1, -1,  1},
      {  1,  0, -1},
      {  1,  0,  0},
      {  1,  0,  1},
      {  1,  1, -1},
      {  1,  1,  0},
      {  1,  1,  1},
    };
    int key0 = (int)floor(posx/box_size);
    int key1 = (int)floor(posy/box_size);
    int key2 = (int)floor(posz/box_size);
    for(boxn = 0; boxn < 27; ++boxn)
    {
      char key[128];
      PyObject *alist;
      snprintf(key, sizeof(key), "%d %d %d", key0+nbor_boxes[boxn][0],
               key1+nbor_boxes[boxn][1], key2+nbor_boxes[boxn][2]);
      alist = PyDict_GetItemString(boxes, key);
      if(!alist && i == -1) printf("none in list at %s, atom %d\n", key, i);
      if(alist)
      {
	int n2 = (int)PyObject_Length(alist);
	if(i == -1) printf("%3d in list at %s\n", n2, key);
	for(j = 0; j < n2; ++j)
	{
	  PyObject * py_i2 = PyList_GetItem(alist, (Py_ssize_t)j);
	  int i2 = PyInt_AsLong(py_i2);
	  if(i2 != i)
	  {
	    PyObject *apos = PyList_GetItem(atom_data, (Py_ssize_t)i2);
	    double vaax = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)0)) - posx;
	    double vaay = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)1)) - posy;
	    double vaaz = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)2)) - posz;
	    double d2 = vaax*vaax + vaay*vaay + vaaz*vaaz;
	    if(d2 <= max_dist_2)
	    {
	      nbor_data[n_nbors1].index = i2;
	      nbor_data[n_nbors1++].dist2 = d2;
	    }
	  }
	}
      }
    }
  }
  Py_DECREF(boxes);
  return n_nbors1;
}

/*
 * returns a tuple of length 2. The first member of the tuple is a list of
 * positions of surface points. The second member is None or a list of
 * unit vectors.
 */

static PyObject *
surface1atom(PyObject *dummy, PyObject *args)
{
  PyObject *nbors;
  NeighborData *nbor_data;
  PyObject *atom_data;
  PyObject *ai_pos;
  PyObject *points;
  double radius, rad1sq, rad1_2;
  int point_density;
  int return_unit_pts;	/* should I return a list of unit points? */
  PyObject *unit_points = NULL;
  PyObject *ret_tup;
  int n_tess, n_nbors;
  Py_ssize_t n_atoms;
  int i, j, aindex;
  double *vx, *vy, *vz, *thresh;
  double aposx, aposy, aposz;
  int last_index = 0;
  static TPoint *tesselations = NULL;
  static int last_point_density = -1;

  if (!PyArg_ParseTuple(args, "OiOOdii", &nbors, &aindex, &atom_data,
			&ai_pos, &radius, &point_density, &return_unit_pts))
    return NULL;
  if(PyObject_Length(ai_pos) < 3)
  {
    PyErr_SetString(PyExc_TypeError, "3rd argument must be a tuple of 3 floats");
    return NULL;
  }
  if(point_density != last_point_density)
  {
    if(tesselations)
      free(tesselations);
    last_point_density = point_density;
    tesselations = tesselate(point_density);
    if(!tesselations)
    {
      PyErr_SetString(PyExc_ValueError,
		      "point_density invalid, must be 2**(2*N) + 2, where N > 1");
      return NULL;
    }
  }
  ret_tup = PyTuple_New((Py_ssize_t)2);
  if(return_unit_pts)
    unit_points = PyList_New((Py_ssize_t)0);
  n_atoms = PyObject_Length(atom_data);
  nbor_data = (NeighborData*)malloc(n_atoms*sizeof(nbor_data[0]));
  n_nbors = nbor_data_1_atom(nbors, aindex, atom_data, nbor_data);
  aposx = PyFloat_AsDouble(PyTuple_GetItem(ai_pos, (Py_ssize_t)0));
  aposy = PyFloat_AsDouble(PyTuple_GetItem(ai_pos, (Py_ssize_t)1));
  aposz = PyFloat_AsDouble(PyTuple_GetItem(ai_pos, (Py_ssize_t)2));
  points = PyList_New((Py_ssize_t)0);
  n_tess = point_density;
  vx = (double *)malloc(n_nbors*sizeof(*vx));
  vy = (double *)malloc(n_nbors*sizeof(*vy));
  vz = (double *)malloc(n_nbors*sizeof(*vz));
  thresh = (double *)malloc(n_nbors*sizeof(*thresh));

  rad1sq = radius*radius;
  rad1_2 = 2*radius;
  for(i = 0; i < n_nbors; ++i)
  {
    double r2;
    PyObject *apos;
    apos = PyList_GetItem(atom_data, (Py_ssize_t)nbor_data[i].index);
    if(!apos) return NULL;
    vx[i] = PyFloat_AsDouble(PyTuple_GetItem(apos, (Py_ssize_t)0)) - aposx;
    vy[i] = PyFloat_AsDouble(PyTuple_GetItem(apos, (Py_ssize_t)1)) - aposy;
    vz[i] = PyFloat_AsDouble(PyTuple_GetItem(apos, (Py_ssize_t)2)) - aposz;
    r2 = PyFloat_AsDouble(PyTuple_GetItem(apos, (Py_ssize_t)3));
    thresh[i] = (nbor_data[i].dist2 + rad1sq - r2*r2) / rad1_2;
  }
  for(i = 0; i < n_tess; ++i)
  {
    int buried = 0;
    double ptx = tesselations[i][0];
    double pty = tesselations[i][1];
    double ptz = tesselations[i][2];
    for(j = last_index; j < n_nbors; ++j)
    {
      if(ptx*vx[j] + pty*vy[j] + ptz*vz[j] > thresh[j])
      {
	buried = 1;
	last_index = j;
	break;
      }
#if 0
      else
      {
	double l2 = sqrt(ptx*ptx + pty*pty + ptz*ptz)
	  *sqrt(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]);
	double cosa = (ptx*vx[j] + pty*vy[j] + ptz*vz[j])/l2;
	double dist41 = sqrt((radius*ptx+0.0933)*(radius*ptx+0.0933)
			     + (radius*pty-0.2487)*(radius*pty-0.2487)
			     + (radius*ptz+0.0903)*(radius*ptz+0.0903));
	double dist28 = sqrt(radius*ptx*radius*ptx + radius*pty*radius*pty + radius*ptz*radius*ptz);
	if(cosa > 0.75 && fabs(ptx) < 0.02 && pty > 0.93 && fabs(ptz - 0.195) < 0.05 && ++cnt0 < 1000)
	  printf("%2d %2d %5.2f*%5.2f + %5.2f*%5.2f + %5.2f*%5.2f > %4.2f (%4.2f) %5.2f %6.2f %6.2f\n",
	       i,j,ptx*10, vx[j]*10, pty*10, vy[j]*10,
	       ptz*10, vz[j]*10, thresh[j]*10, 10*(ptx*vx[j] + pty*vy[j] + ptz*vz[j]), cosa, dist28, dist41);
      }
#endif
    }
    if(!buried)
    {
      for(j = 0; j < last_index; ++j)
      {
	if(ptx*vx[j] + pty*vy[j] + ptz*vz[j] > thresh[j])
	{
	  buried = 1;
	  last_index = j;
	  break;
	}
#if 0
      else
      {
	double l2 = sqrt(ptx*ptx + pty*pty + ptz*ptz)
	  *sqrt(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]);
	double cosa = (ptx*vx[j] + pty*vy[j] + ptz*vz[j])/l2;
	if(cosa > 0.75 && fabs(ptx) < 0.02 && pty > 0.93 && fabs(ptz - 0.195) < 0.05 && ++cnt0 < 1000)
	  printf("%4d %3d %5.2f*%6.2f + %5.2f*%6.2f + %5.2f*%6.2f > %5.2f (%5.2f) %6.2f\n",
	       i,j,ptx*10, vx[j]*10, pty*10, vy[j]*10,
	       ptz*10, vz[j]*10, thresh[j]*10, 10*(ptx*vx[j] + pty*vy[j] + ptz*vz[j]), cosa);
      }
#endif
      }
    }
    if(!buried)
    {
      PyObject *new_pt = PyTuple_New((Py_ssize_t)3);
      PyTuple_SetItem(new_pt, (Py_ssize_t)0,
		      PyFloat_FromDouble(radius*ptx + aposx));
      PyTuple_SetItem(new_pt, (Py_ssize_t)1,
		      PyFloat_FromDouble(radius*pty + aposy));
      PyTuple_SetItem(new_pt, (Py_ssize_t)2,
		      PyFloat_FromDouble(radius*ptz + aposz));
      PyList_Append(points, new_pt);
      Py_DECREF(new_pt);
      if(return_unit_pts)
      {
	PyObject *new_pt = PyTuple_New((Py_ssize_t)3);
	PyTuple_SetItem(new_pt, (Py_ssize_t)0, PyFloat_FromDouble(ptx));
	PyTuple_SetItem(new_pt, (Py_ssize_t)1, PyFloat_FromDouble(pty));
	PyTuple_SetItem(new_pt, (Py_ssize_t)2, PyFloat_FromDouble(ptz));
	PyList_Append(unit_points, new_pt);
	Py_DECREF(new_pt);
      }
    }
  }
  free(vx);
  free(vy);
  free(vz);
  free(thresh);
  free(nbor_data);

  /*
  surface_(&total_area, areas, dareas, doublereal *radius, doublereal *weight, doublereal *probe);
  */
  PyTuple_SetItem(ret_tup, (Py_ssize_t)0, points);
  if(return_unit_pts)
    PyTuple_SetItem(ret_tup, (Py_ssize_t)1, unit_points);
  else
  {
    Py_INCREF(Py_None);
    PyTuple_SetItem(ret_tup, (Py_ssize_t)1, Py_None);
  }
  return ret_tup;
}

/* not sure whether this is still usefull; may be obsoleted by the next
   function, which uses less memory due to lazy evaluation.
 */
static PyObject*
FindNeighbors(PyObject *dummy, PyObject *args)
{
  PyObject *atoms;
  double radius;
  double max_rad;
  PyObject *atom_data;
  double max_dist_2;
  PyObject *nbors, *boxes;
  Py_ssize_t n_atoms;
  int i, j;
  double box_size;
  PyObject **nlist3, *nlist;

  if (!PyArg_ParseTuple(args, "OddOd", &atoms, &radius,
			&max_rad, &atom_data, &max_dist_2))
    return NULL;
  n_atoms = PyObject_Length(atoms);
  nbors = PyTuple_New(n_atoms);
  nlist3 = (PyObject **)malloc(n_atoms*sizeof(PyObject*));
  boxes = PyDict_New();
  box_size = 2*(max_rad + radius);
  printf("box_size %.2f %.2f %.2f\n", box_size*10, max_rad*10, radius*10);
  for(i = 0; i < n_atoms; ++i)
  {
    PyObject *v, *tmp;
    char key[128];
    PyObject *pos = PyList_GetItem(atom_data, (Py_ssize_t)i);
    snprintf(key, sizeof(key), "%d %d %d",
             (int)floor(PyFloat_AsDouble(PyTuple_GetItem(pos,
						     (Py_ssize_t)0))/box_size),
             (int)floor(PyFloat_AsDouble(PyTuple_GetItem(pos,
						     (Py_ssize_t)1))/box_size),
             (int)floor(PyFloat_AsDouble(PyTuple_GetItem(pos,
						     (Py_ssize_t)2))/box_size));
    v = PyDict_GetItemString(boxes, key);
    if(!v)
      PyDict_SetItemString(boxes, key, v = PyList_New((Py_ssize_t)0));
    PyList_Append(v, tmp = PyInt_FromLong(i));
    Py_DECREF(tmp);
    if(0)
      printf("key %12.12s %3d %6.2f %6.2f %6.2f\n",
	     key, (int)PyObject_Length(v),
	    PyFloat_AsDouble(PyTuple_GetItem(pos, (Py_ssize_t)0)),
	    PyFloat_AsDouble(PyTuple_GetItem(pos, (Py_ssize_t)1)),
	    PyFloat_AsDouble(PyTuple_GetItem(pos, (Py_ssize_t)2)));
  }
  for(i = 0; i < n_atoms; ++i)
  {
    PyObject *pos1 = PyList_GetItem(atom_data, (Py_ssize_t)i);
    double posx = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)0));
    double posy = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)1));
    double posz = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)2));
    int n_nbors1 = 0;
    int boxn;
    static const int nbor_boxes[27][3] = {
      { -1, -1, -1},
      { -1, -1,  0},
      { -1, -1,  1},
      { -1,  0, -1},
      { -1,  0,  0},
      { -1,  0,  1},
      { -1,  1, -1},
      { -1,  1,  0},
      { -1,  1,  1},
      {  0, -1, -1},
      {  0, -1,  0},
      {  0, -1,  1},
      {  0,  0, -1},
      {  0,  0,  0},
      {  0,  0,  1},
      {  0,  1, -1},
      {  0,  1,  0},
      {  0,  1,  1},
      {  1, -1, -1},
      {  1, -1,  0},
      {  1, -1,  1},
      {  1,  0, -1},
      {  1,  0,  0},
      {  1,  0,  1},
      {  1,  1, -1},
      {  1,  1,  0},
      {  1,  1,  1},
    };
    int key0 = (int)floor(posx/box_size);
    int key1 = (int)floor(posy/box_size);
    int key2 = (int)floor(posz/box_size);
    for(boxn = 0; boxn < 27; ++boxn)
    {
      char key[128];
      PyObject *alist;
      snprintf(key, sizeof(key), "%d %d %d", key0+nbor_boxes[boxn][0],
               key1+nbor_boxes[boxn][1], key2+nbor_boxes[boxn][2]);
      alist = PyDict_GetItemString(boxes, key);
      if(!alist && i == -1) printf("none in list at %s\n", key);
      if(alist)
      {
	int n2 = (int)PyObject_Length(alist);
	if(i == -1) printf("%3d in list at %s\n", n2, key);
	for(j = 0; j < n2; ++j)
	{
	  PyObject * py_i2 = PyList_GetItem(alist, (Py_ssize_t)j);
	  int i2 = PyInt_AsLong(py_i2);
	  if(i2 != i)
	  {
	    PyObject *apos = PyList_GetItem(atom_data, (Py_ssize_t)i2);
	    double vaax = PyFloat_AsDouble(PyTuple_GetItem(apos,
                                                         (Py_ssize_t)0)) - posx;
	    double vaay = PyFloat_AsDouble(PyTuple_GetItem(apos,
                                                         (Py_ssize_t)1)) - posy;
	    double vaaz = PyFloat_AsDouble(PyTuple_GetItem(apos,
                                                         (Py_ssize_t)2)) - posz;
	    double d2 = vaax*vaax + vaay*vaay + vaaz*vaaz;
	    if(d2 <= max_dist_2)
	    {
	      PyObject *tup1 = PyTuple_New((Py_ssize_t)2);
	      Py_INCREF(py_i2);
	      PyTuple_SetItem(tup1, (Py_ssize_t)0, py_i2);
	      PyTuple_SetItem(tup1, (Py_ssize_t)1, PyFloat_FromDouble(d2));
	      nlist3[n_nbors1++] = tup1;
	    }
	  }
	}
      }
    }
    nlist = PyTuple_New((Py_ssize_t)n_nbors1);
    for(j = 0; j < n_nbors1; ++j)
    {
      PyTuple_SetItem(nlist, (Py_ssize_t)j, nlist3[j]);
    }
    PyTuple_SetItem(nbors, (Py_ssize_t)i, nlist);
  }
  free(nlist3);
  Py_DECREF(boxes);
  return nbors;
}

static PyObject*
FindNeighborsOfAtom(PyObject *dummy, PyObject *args)
{
  PyObject *atoms;
  PyObject *atom_data;
  double max_dist_2;
  PyObject *boxes;
  Py_ssize_t n_atoms;
  int i, j;
  double box_size;
  PyObject **nlist3, *nlist;

  if (!PyArg_ParseTuple(args, "OiOdO", &atoms, &i, &boxes, &box_size, &atom_data))
    return NULL;
  n_atoms = PyObject_Length(atoms);
  nlist3 = (PyObject **)malloc(n_atoms*sizeof(PyObject*));
  max_dist_2 = box_size*box_size;
  {
    PyObject *pos1 = PyList_GetItem(atom_data, (Py_ssize_t)i);
    double posx = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)0));
    double posy = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)1));
    double posz = PyFloat_AsDouble(PyTuple_GetItem(pos1, (Py_ssize_t)2));
    int n_nbors1 = 0;
    int boxn;
    static const int nbor_boxes[27][3] = {
      { -1, -1, -1},
      { -1, -1,  0},
      { -1, -1,  1},
      { -1,  0, -1},
      { -1,  0,  0},
      { -1,  0,  1},
      { -1,  1, -1},
      { -1,  1,  0},
      { -1,  1,  1},
      {  0, -1, -1},
      {  0, -1,  0},
      {  0, -1,  1},
      {  0,  0, -1},
      {  0,  0,  0},
      {  0,  0,  1},
      {  0,  1, -1},
      {  0,  1,  0},
      {  0,  1,  1},
      {  1, -1, -1},
      {  1, -1,  0},
      {  1, -1,  1},
      {  1,  0, -1},
      {  1,  0,  0},
      {  1,  0,  1},
      {  1,  1, -1},
      {  1,  1,  0},
      {  1,  1,  1},
    };
    int key0 = (int)floor(posx/box_size);
    int key1 = (int)floor(posy/box_size);
    int key2 = (int)floor(posz/box_size);
    for(boxn = 0; boxn < 27; ++boxn)
    {
      char key[128];
      PyObject *alist;
      snprintf(key, sizeof(key),
               "%d %d %d", key0+nbor_boxes[boxn][0],
               key1+nbor_boxes[boxn][1], key2+nbor_boxes[boxn][2]);
      alist = PyDict_GetItemString(boxes, key);
      if(!alist && i == -1) printf("none in list at %s\n", key);
      if(alist)
      {
	int n2 = (int)PyObject_Length(alist);
	if(i == -1) printf("%3d in list at %s\n", n2, key);
	for(j = 0; j < n2; ++j)
	{
	  PyObject * py_i2 = PyList_GetItem(alist, (Py_ssize_t)j);
	  int i2 = PyInt_AsLong(py_i2);
	  if(i2 != i)
	  {
	    PyObject *apos = PyList_GetItem(atom_data, (Py_ssize_t)i2);
	    double vaax = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)0)) - posx;
	    double vaay = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)1)) - posy;
	    double vaaz = PyFloat_AsDouble(PyTuple_GetItem(apos,
							(Py_ssize_t)2)) - posz;
	    double d2 = vaax*vaax + vaay*vaay + vaaz*vaaz;
	    if(d2 <= max_dist_2)
	    {
	      PyObject *tup1 = PyTuple_New((Py_ssize_t)2);
	      Py_INCREF(py_i2);
	      PyTuple_SetItem(tup1, (Py_ssize_t)0, py_i2);
	      PyTuple_SetItem(tup1, (Py_ssize_t)1, PyFloat_FromDouble(d2));
	      nlist3[n_nbors1++] = tup1;
	    }
	  }
	}
      }
    }
    nlist = PyTuple_New((Py_ssize_t)n_nbors1);
    for(j = 0; j < n_nbors1; ++j)
    {
      PyTuple_SetItem(nlist, (Py_ssize_t)j, nlist3[j]);
    }
  }
  free(nlist3);
  return nlist;
}

/*
 * List of functions defined in the module
 */

static PyMethodDef surface_methods[] = {
  {"surface1atom", surface1atom, METH_VARARGS},
  {"FindNeighbors", FindNeighbors, METH_VARARGS},
  {"FindNeighborsOfAtom", FindNeighborsOfAtom, METH_VARARGS},
  {NULL, NULL}		/* sentinel */
};

/* Initialization function for the module */

DL_EXPORT(void)
#ifdef IBMPC
__declspec(dllexport)
#endif
initMMTK_surface(void)
{
  /* Create the module and add the functions */
  Py_InitModule("MMTK_surface", surface_methods);

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_surface");
}
