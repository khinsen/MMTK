/* Low-level force field calculations: non-bonded interactions
 *
 * Written by Konrad Hinsen
 * last revision: 2007-4-24
 */

#define NO_IMPORT
#define _FORCEFIELD_MODULE
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_MMTKFF_API

#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

#define THREAD_DEBUG 0

/* PMTA definitions */

#ifdef WITH_DPMTA
#include "dpmta.h"
#endif

/* Utility functions */

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))


/* Nonbonded list update */

/* Sort all atoms into small subboxes such that atoms in any box
 * can interact only with particles in a small shell of neighbouring boxes */
int
nblist_update(PyNonbondedListObject *nblist, int natoms,
	      double *coordinates, double *geometry_data)
{
  vector3 *x = (vector3 *)coordinates;
  long *subset = (long *)((PyArrayObject *)nblist->atom_subset)->data;
  int n_sub = ((PyArrayObject *)nblist->atom_subset)->dimensions[0];
  vector3 box1, box2;
  double box_size[3];
  int *p;
  int i, ix, iy, iz, minx, miny, minz, maxx, maxy, maxz, n;

  if (nblist->box_number == NULL) {
    nblist->box_number = (int *)malloc(2*natoms*sizeof(int));
    if (nblist->box_number == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    nblist->box_atoms = nblist->box_number + natoms;
  }
  
  nblist->universe_spec->correction_function(x, natoms, geometry_data);
  nblist->universe_spec->bounding_box_function(&box1, &box2, x, natoms,
					       geometry_data);
  nblist->lastx = x;
#if 0
  printf("box1: %lf, %lf, %lf\n", box1[0], box1[1], box1[2]);
  printf("box2: %lf, %lf, %lf\n", box2[0], box2[1], box2[2]);
  printf("cutoff: %lf\n", nblist->cutoff);
#endif
  if (nblist->cutoff > 0. && (!nblist->universe_spec->is_periodic
			      || nblist->universe_spec->is_orthogonal)) {
    int done = 0;
    double factor = 1.;
    while (!done) {
      int nboxes;
      nblist->box_count[0] = (int)(NBLIST_NEIGHBORS*(box2[0]-box1[0])
				   /(factor*nblist->cutoff));
      nblist->box_count[1] = (int)(NBLIST_NEIGHBORS*(box2[1]-box1[1])
				   /(factor*nblist->cutoff));
      nblist->box_count[2] = (int)(NBLIST_NEIGHBORS*(box2[2]-box1[2])
				   /(factor*nblist->cutoff));
      if (nblist->box_count[0] == 0) nblist->box_count[0] = 1;
      if (nblist->box_count[1] == 0) nblist->box_count[1] = 1;
      if (nblist->box_count[2] == 0) nblist->box_count[2] = 1;
      nboxes = nblist->box_count[0]*nblist->box_count[1]*nblist->box_count[2];
      if (nboxes > 2*natoms)
	factor *= 1.1;
      else
	done = 1;
    }
  }
  else
    nblist->box_count[0] = nblist->box_count[1] = nblist->box_count[2] = 1;
  box_size[0] = (box2[0]-box1[0])/nblist->box_count[0];
  box_size[1] = (box2[1]-box1[1])/nblist->box_count[1];
  box_size[2] = (box2[2]-box1[2])/nblist->box_count[2];
  if (box_size[0] == 0.) box_size[0] = 1.;
  if (box_size[1] == 0.) box_size[1] = 1.;
  if (box_size[2] == 0.) box_size[2] = 1.;
#if 0
  printf("division: %d/%d/%d\n", nblist->box_count[0],
	 nblist->box_count[1], nblist->box_count[2]);
  printf("cell size: %lf, %lf, %lf\n", box_size[0], box_size[1], box_size[2]);
#endif

  nblist->neighbors[0][0] = 0;
  nblist->neighbors[0][1] = 0;
  nblist->neighbors[0][2] = 0;
  i = 1;
  minx = -(int)((nblist->cutoff+box_size[0])/box_size[0]);
  miny = -(int)((nblist->cutoff+box_size[1])/box_size[1]);
  minz = -(int)((nblist->cutoff+box_size[2])/box_size[2]);
  maxx = -minx+1;
  maxy = -miny+1;
  maxz = -minz+1;
  if (nblist->universe_spec->is_periodic) {
    maxx = min(maxx, (nblist->box_count[0]+1)/2);
    maxy = min(maxy, (nblist->box_count[1]+1)/2);
    maxz = min(maxz, (nblist->box_count[2]+1)/2);
    minx = max(minx, maxx-nblist->box_count[0]);
    miny = max(miny, maxy-nblist->box_count[1]);
    minz = max(minz, maxz-nblist->box_count[2]);
  }
  else {
    maxx = min(maxx, nblist->box_count[0]);
    maxy = min(maxy, nblist->box_count[1]);
    maxz = min(maxz, nblist->box_count[2]);
    minx = max(minx, 1-nblist->box_count[0]);
    miny = max(miny, 1-nblist->box_count[1]);
    minz = max(minz, 1-nblist->box_count[2]);
  }
#if 0
  printf("Box neighbor list:\n");
  printf("  minx: %d\tmaxx: %d\n", minx, maxx);
  printf("  miny: %d\tmaxy: %d\n", miny, maxy);
  printf("  minz: %d\tmaxz: %d\n", minz, maxz);
#endif
  for (ix = minx; ix < maxx; ix++)
    for (iy = miny; iy < maxy; iy++)
      for (iz = minz; iz < maxz; iz++)
	if (!(ix == 0 && iy == 0 && iz == 0)) {
	  double dx = (abs(ix)-1.)*box_size[0];
	  double dy = (abs(iy)-1.)*box_size[1];
	  double dz = (abs(iz)-1.)*box_size[2];
	  if (dx < 0.) dx = 0.;
	  if (dy < 0.) dy = 0.;
	  if (dz < 0.) dz = 0.;
	  if (dx*dx+dy*dy+dz*dz <= sqr(nblist->cutoff)) {
	    nblist->neighbors[i][0] = ix;
	    nblist->neighbors[i][1] = iy;
	    nblist->neighbors[i][2] = iz;
#if 0
	    printf(" %d: %d/%d/%d\n", i, ix, iy, iz);
#endif
	    i++;
	  }
#if 0
	  else {
	    printf(" %d/%d/%d -> %f/%f/%f -> %f\n",
		   ix, iy, iz, dx, dy, dz, sqrt(dx*dx+dy*dy+dz*dz));
	  }
#endif
	}
  nblist->nneighbors = i;
#if 0
  printf("Box neighbors: %d\n", i);
#endif

  nblist->nboxes =
         nblist->box_count[0]*nblist->box_count[1]*nblist->box_count[2];
  if (nblist->nboxes > nblist->allocated_boxes) {
    free(nblist->boxes);
    nblist->boxes = malloc(nblist->nboxes*sizeof(nbbox));
    if (nblist->boxes == NULL) {
      PyErr_NoMemory();
      return 0;
    }
    nblist->allocated_boxes = nblist->nboxes;
  }

  i = 0;
  for (iz = 0; iz < nblist->box_count[2]; iz++)
    for (iy = 0; iy < nblist->box_count[1]; iy++)
      for (ix = 0; ix < nblist->box_count[0]; ix++) {
	nblist->boxes[i].ix = ix;
	nblist->boxes[i].iy = iy;
	nblist->boxes[i].iz = iz;
	nblist->boxes[i].n = 0;
	nblist->boxes[i].i = 0;
	i++;
      }

  n = (n_sub == 0) ? natoms : n_sub;
  for (i = 0; i < n; i++) {
    int box = 0;
    int ai, n;
    ai = (n_sub == 0) ? i : subset[i];
    box = (int)((x[ai][0]-box1[0])/box_size[0]);
    if (box == nblist->box_count[0]) box--;
    n = (int)((x[ai][1]-box1[1])/box_size[1]);
    if (n == nblist->box_count[1]) n--;
    box += nblist->box_count[0]*n;
    n = (int)((x[ai][2]-box1[2])/box_size[2]);
    if (n == nblist->box_count[2]) n--;
    box += nblist->box_count[0]*nblist->box_count[1]*n;
    nblist->box_number[ai] = box;
#if 0
    printf("Atom %d in box %d (%d/%d/%d)\n", ai, box,
	   nblist->boxes[box].ix, nblist->boxes[box].iy, nblist->boxes[box].iz);
#endif 
    nblist->boxes[box].n++;
  }
  p = nblist->box_atoms;
  for (i = 0; i < nblist->nboxes; i++) {
    nblist->boxes[i].atoms = p;
    nblist->boxes[i].i = 0;
    p += nblist->boxes[i].n;
  }
  for (i = 0; i < n; i++) {
    int ai = (n_sub == 0) ? i : subset[i];
    nbbox *box = nblist->boxes + nblist->box_number[ai];
    box->atoms[box->i++] = ai;
  }  
  for (i = 0; i < nblist->nboxes; i++)
    nblist->boxes[i].i = 0;

#if 0
  ix = 0;
  for (i = 0; i < nblist->nboxes; i++) {
    int j;
    for (j = 0; j < nblist->boxes[i].n; j++) {
      int atom = nblist->boxes[i].atoms[j];
      ix++;
      if (nblist->box_number[atom] != i)
	printf("Box number error: %d: %d - %d\n",
	       atom, i, nblist->box_number[atom]);
    }
  }
  if (ix != ((n_sub == 0) ? natoms : n_sub))
    printf("Atom number error: %d - %d\n", ix, natoms);
#endif

  nblist->iterator.state = nblist_start;
  nblist->iterator.n = -1;

  return 1;
}

/* Iterator over nonbonded list */
int
nblist_iterate(PyNonbondedListObject *nblist, struct nblist_iterator *iterator)
{
  long *excluded, *one_four;
  int n_ex, n_14;
  int ix, iy, iz, use_box;

  switch (iterator->state) {

  case nblist_start:
    iterator->n = -1;
    iterator->ibox = -1;
    iterator->jbox = -1;
    iterator->ineighbor = nblist->nneighbors-1;
    iterator->box1 = nblist->boxes;
    iterator->box2 = nblist->boxes;
    iterator->i = iterator->box1->n-1;
    iterator->j = iterator->box2->n-1;
    iterator->state = nblist_continue;
  case nblist_continue:
    iterator->j++;
    if (iterator->j == iterator->box2->n) {
      int lasti = iterator->box1->n - ((iterator->box1==iterator->box2)?1:0);
      iterator->i++;
      if (iterator->i >= lasti) {
	do {
	  iterator->ineighbor++;
	  if (iterator->ineighbor == nblist->nneighbors) {
	    do {
	      iterator->ibox++;
	      if (iterator->ibox == nblist->nboxes) {
		iterator->state = nblist_finished;
		return 0;
	      }
	      iterator->box1 = nblist->boxes+iterator->ibox;
	      iterator->ineighbor = 0;
	    } while (iterator->box1->n == 0);
	  }
	  ix = nblist->neighbors[iterator->ineighbor][0]
	       + iterator->box1->ix;
	  iy = nblist->neighbors[iterator->ineighbor][1]
               + iterator->box1->iy;
	  iz = nblist->neighbors[iterator->ineighbor][2]
               + iterator->box1->iz;
	  use_box = 1;
	  if (nblist->universe_spec->is_periodic) {
	    if (ix < 0) ix += nblist->box_count[0];
	    if (iy < 0) iy += nblist->box_count[1];
	    if (iz < 0) iz += nblist->box_count[2];
	    if (ix >= nblist->box_count[0])
	      ix -= nblist->box_count[0];
	    if (iy >= nblist->box_count[1])
	      iy -= nblist->box_count[1];
	    if (iz >= nblist->box_count[2])
	      iz -= nblist->box_count[2];
	  }
	  else if (ix < 0 || iy < 0 || iz < 0
		   || ix >= nblist->box_count[0]
		   || iy >= nblist->box_count[1]
		   || iz >= nblist->box_count[2]) {
	    use_box = 0;
	  }
	  iterator->jbox = ix + nblist->box_count[0]
	                   *(iy + nblist->box_count[1]*iz);
	  if (iterator->jbox < iterator->ibox)
	    use_box = 0;
	  if (use_box) {
	    if (nblist->boxes[iterator->jbox].n == 0) {
	      use_box = 0;
	    }
	    if (iterator->ibox == iterator->jbox &&
		nblist->boxes[iterator->jbox].n == 1) {
	      use_box = 0;
	    }
	  }
	} while (!use_box);
	iterator->box2 = nblist->boxes+iterator->jbox;
	iterator->i = 0;
      }
      if (iterator->ibox == iterator->jbox)
	iterator->j = iterator->i + 1;
      else
	iterator->j = 0;
    }
    iterator->a1 = iterator->box1->atoms[iterator->i];
    iterator->a2 = iterator->box2->atoms[iterator->j];
    iterator->n++;
    return 1;
    break;

  case nblist_start_excluded:
    iterator->i = -2;
    iterator->state = nblist_continue_excluded;
  case nblist_continue_excluded:
    excluded = (long *)((PyArrayObject *)nblist->excluded_pairs)->data;
    n_ex = 2*((PyArrayObject *)nblist->excluded_pairs)->dimensions[0];
    iterator->i += 2;
    if (iterator->i == n_ex) {
      iterator->state = nblist_finished;
      return 0;
    }
    iterator->a1 = excluded[iterator->i];
    iterator->a2 = excluded[iterator->i+1];
    return 1;
    break;

  case nblist_start_14:
    iterator->i = -2;
    iterator->state = nblist_continue_14;
  case nblist_continue_14:
    one_four = (long *)((PyArrayObject *)nblist->one_four_pairs)->data;
    n_14 = 2*((PyArrayObject *)nblist->one_four_pairs)->dimensions[0];
    iterator->i += 2;
    if (iterator->i == n_14) {
      iterator->state = nblist_finished;
      return 0;
    }
    iterator->a1 = one_four[iterator->i];
    iterator->a2 = one_four[iterator->i+1];
    return 1;
    break;

  case nblist_finished:
    return 0;
    break;

  }
}

/* Evaluator for all non-bonded interactions */

#ifdef GRADIENTFN
#define pair_term(ljfactor, esfactor, ewaldfactor) \
{ \
  vector3 rij; \
  double r_sq, r; \
  double deriv = 0., deriv2 = 0.; \
  (*d_fn)(rij, x[a2], x[a1], distance_data); \
  r_sq = vector_length_sq(rij); \
  if (es_flag || ewald_flag) \
    r = sqrt(r_sq); \
  \
  if (lj_flag && (lj_cutoff_sq == 0. || r_sq <= lj_cutoff_sq)) { \
    int type1 = lj_type[a1]; \
    int type2 = lj_type[a2]; \
    if (type1 >= 0 && type2 >= 0) { \
      double eps = ljfactor*eps_sigma[2*(ntypes*type1+type2)]; \
      double sigma = eps_sigma[2*(ntypes*type1+type2)+1]; \
      double sr2 = sqr(sigma)/r_sq; \
      double sr6 = cube(sr2); \
      double sr12 = sqr(sr6); \
      double e = 4.*eps*(sr12-sr6); \
      double v = 24.*eps*(2.*sr12-sr6); \
      if (r_sq < 0.02) { \
        lj_energy2 += e; \
        lj_virial2 += v; \
      } \
      else { \
        lj_energy1 += e; \
        lj_virial1 += v; \
      } \
      deriv += -24.*eps*(2.*sr12-sr6)/r_sq; \
      deriv2 += 24.*eps*(26.*sr12-7.*sr6)/r_sq; \
    } \
  } \
  \
  if (es_flag && (es_cutoff_sq == 0. || r_sq <= es_cutoff_sq)) { \
    double qiqj = esfactor*charge[a1]*charge[a2]*electrostatic_energy_factor; \
    double term = qiqj*(1./r-es_inv_cutoff); \
    es_energy += term; \
    deriv += -qiqj*(1./r_sq-sqr(es_inv_cutoff))/r; \
    deriv2 += 2.*qiqj*(1./(r*r_sq)-es_inv_cutoff*sqr(es_inv_cutoff)); \
  } \
  \
  if (ewald_flag && (ewald_cutoff_sq == 0. || r_sq<=ewald_cutoff_sq)) { \
    double f1 = erfc(beta*r); \
    double qiqj = ewaldfactor*charge[a1]*charge[a2] \
                                * electrostatic_energy_factor; \
    double ef = 2./sqrt(M_PI); \
    ewald_energy += qiqj*(f1/r-erfc_cutoff*ewald_inv_cutoff); \
    if (energy->gradients != NULL) \
      deriv += -qiqj*(f1/r_sq-erfc_cutoff*ewald_inv_cutoff \
		   +ef*beta*(exp(-beta*beta*r_sq)/r \
		     -ewald_inv_cutoff*exp(-beta*beta*ewald_cutoff_sq)))/r; \
    if (energy->force_constants != NULL) { \
      deriv2 += 2.*qiqj*(f1/(r*r_sq) - erfc_cutoff*cube(ewald_inv_cutoff) \
			 + ef*beta*exp(-beta*beta*r_sq) \
			 * (beta*beta+1./r_sq)); \
      if (ewald_inv_cutoff > 0.) \
	deriv2 -= 2.*qiqj*ef*beta*exp(-beta*beta*ewald_cutoff_sq) \
	  * (beta*beta+sqr(ewald_inv_cutoff)); \
    } \
  } \
  \
  if (energy->gradients != NULL) { \
    vector3 grad; \
    grad[0] = deriv*rij[0]; \
    grad[1] = deriv*rij[1]; \
    grad[2] = deriv*rij[2]; \
    if (energy->gradient_fn != NULL) { \
      (*energy->gradient_fn)(energy, a1, grad); \
      vector_changesign(grad); \
      (*energy->gradient_fn)(energy, a2, grad); \
    } \
    else { \
      vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data; \
      f[a1][0] += grad[0]; \
      f[a1][1] += grad[1]; \
      f[a1][2] += grad[2]; \
      f[a2][0] -= grad[0]; \
      f[a2][1] -= grad[1]; \
      f[a2][2] -= grad[2]; \
    } \
  } \
  if (energy->force_constants != NULL) { \
    add_pair_fc(energy, a1, a2, rij, r_sq, deriv, deriv2); \
  } \
}
#else
#define pair_term(ljfactor, esfactor, ewaldfactor) \
{ \
  vector3 rij; \
  double r_sq, r; \
  double deriv = 0., deriv2 = 0.; \
  (*d_fn)(rij, x[a2], x[a1], distance_data); \
  r_sq = vector_length_sq(rij); \
  if (es_flag || ewald_flag) \
    r = sqrt(r_sq); \
  \
  if (lj_flag && (lj_cutoff_sq == 0. || r_sq <= lj_cutoff_sq)) { \
    int type1 = lj_type[a1]; \
    int type2 = lj_type[a2]; \
    if (type1 >= 0 && type2 >= 0) { \
      double eps = ljfactor*eps_sigma[2*(ntypes*type1+type2)]; \
      double sigma = eps_sigma[2*(ntypes*type1+type2)+1]; \
      double sr2 = sqr(sigma)/r_sq; \
      double sr6 = cube(sr2); \
      double sr12 = sqr(sr6); \
      double e = 4.*eps*(sr12-sr6); \
      double v = 24.*eps*(2.*sr12-sr6); \
      if (r_sq < 0.02) { \
        lj_energy2 += e; \
        lj_virial2 += v; \
      } \
      else { \
        lj_energy1 += e; \
        lj_virial1 += v; \
      } \
      deriv += -24.*eps*(2.*sr12-sr6)/r_sq; \
      deriv2 += 24.*eps*(26.*sr12-7.*sr6)/r_sq; \
    } \
  } \
  \
  if (es_flag && (es_cutoff_sq == 0. || r_sq <= es_cutoff_sq)) { \
    double qiqj = esfactor*charge[a1]*charge[a2]*electrostatic_energy_factor; \
    double term = qiqj*(1./r-es_inv_cutoff); \
    es_energy += term; \
    deriv += -qiqj*(1./r_sq-sqr(es_inv_cutoff))/r; \
    deriv2 += 2.*qiqj*(1./(r*r_sq)-es_inv_cutoff*sqr(es_inv_cutoff)); \
  } \
  \
  if (ewald_flag && (ewald_cutoff_sq == 0. || r_sq<=ewald_cutoff_sq)) { \
    double f1 = erfc(beta*r); \
    double qiqj = ewaldfactor*charge[a1]*charge[a2] \
                                * electrostatic_energy_factor; \
    double ef = 2./sqrt(M_PI); \
    ewald_energy += qiqj*(f1/r-erfc_cutoff*ewald_inv_cutoff); \
    if (energy->gradients != NULL) \
      deriv += -qiqj*(f1/r_sq-erfc_cutoff*ewald_inv_cutoff \
		   +ef*beta*(exp(-beta*beta*r_sq)/r \
		     -ewald_inv_cutoff*exp(-beta*beta*ewald_cutoff_sq)))/r; \
    if (energy->force_constants != NULL) { \
      deriv2 += 2.*qiqj*(f1/(r*r_sq) - erfc_cutoff*cube(ewald_inv_cutoff) \
			 + ef*beta*exp(-beta*beta*r_sq) \
			 * (beta*beta+1./r_sq)); \
      if (ewald_inv_cutoff > 0.) \
	deriv2 -= 2.*qiqj*ef*beta*exp(-beta*beta*ewald_cutoff_sq) \
	  * (beta*beta+sqr(ewald_inv_cutoff)); \
    } \
  } \
  \
  if (energy->gradients != NULL) { \
    vector3 grad; \
    grad[0] = deriv*rij[0]; \
    grad[1] = deriv*rij[1]; \
    grad[2] = deriv*rij[2]; \
    { \
      vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data; \
      f[a1][0] += grad[0]; \
      f[a1][1] += grad[1]; \
      f[a1][2] += grad[2]; \
      f[a2][0] -= grad[0]; \
      f[a2][1] -= grad[1]; \
      f[a2][2] -= grad[2]; \
    } \
  } \
  if (energy->force_constants != NULL) { \
    add_pair_fc(energy, a1, a2, rij, r_sq, deriv, deriv2); \
  } \
}
#endif

void
nonbonded_evaluator(PyFFEnergyTermObject *self,
		    PyFFEvaluatorObject *eval,
		    energy_spec *input,
		    energy_data *energy)
{
  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self->data[0];
  long *excluded = (long *)((PyArrayObject *)nblist->excluded_pairs)->data;
  int n_ex = 2*((PyArrayObject *)nblist->excluded_pairs)->dimensions[0];
  long *one_four = (long *)((PyArrayObject *)nblist->one_four_pairs)->data;
  int n_14 = 2*((PyArrayObject *)nblist->one_four_pairs)->dimensions[0];
  distance_fn *d_fn = nblist->universe_spec->distance_function;
  double *distance_data = nblist->universe_spec->geometry_data;
  vector3 *x = (vector3 *)input->coordinates->data;

  double *eps_sigma;  int ntypes; long *lj_type;
  double lj_cutoff_sq, lj_one_four;
  double *charge; double es_cutoff_sq, es_inv_cutoff, es_one_four;
  double ewald_cutoff_sq, ewald_inv_cutoff, beta, erfc_cutoff, ewald_one_four;
  double lj_energy1, lj_energy2, lj_virial1, lj_virial2;
  double es_energy, ewald_energy;
  int lj_flag, es_flag, ewald_flag, ewald_recip_flag;
  int ibox, jbox, slicecounter, k;
#if THREAD_DEBUG
  int paircount = 0;
#endif

  PyFFEnergyTermObject *lj_ev, *es_ev, *ewald_ev;
  lj_ev = (PyFFEnergyTermObject *)self->data[1];
  es_ev = (PyFFEnergyTermObject *)self->data[2];
  ewald_ev = (PyFFEnergyTermObject *)self->data[3];
  lj_flag = lj_ev != NULL;
  es_flag = es_ev != NULL;
  ewald_flag = ewald_ev != NULL;

  lj_energy1 = lj_energy2 = lj_virial1 = lj_virial2 = 0.;
  es_energy = 0.;
  ewald_energy = 0.;

  if (lj_flag) {
    eps_sigma = (double *)((PyArrayObject *)lj_ev->data[1])->data;
    ntypes = ((PyArrayObject *)lj_ev->data[1])->dimensions[0];
    lj_type = (long *)((PyArrayObject *)lj_ev->data[2])->data;
    lj_cutoff_sq = sqr(lj_ev->param[0]);
    lj_one_four = lj_ev->param[1]-1.;
  }

  if (es_flag) {
    charge = (double *)((PyArrayObject *)es_ev->data[1])->data;
    es_cutoff_sq = sqr(es_ev->param[0]);
    es_inv_cutoff = (es_ev->param[0] == 0.) ? 0. : 1./es_ev->param[0];
    es_one_four = es_ev->param[1]-1.;
  }

  if (ewald_flag) {
    charge = (double *)((PyArrayObject *)ewald_ev->data[1])->data;
    ewald_cutoff_sq = sqr(ewald_ev->param[0]);
    ewald_inv_cutoff = (ewald_ev->param[0] == 0.) ? 0. : 1./ewald_ev->param[0];
    if (ewald_ev->param[0] > 0. && ewald_ev->param[3] == 0.)
      ewald_inv_cutoff = 1./ewald_ev->param[0];
    else
      ewald_inv_cutoff = 0.;
    beta = ewald_ev->param[2];
    erfc_cutoff = erfc(beta*ewald_ev->param[0]);
    ewald_recip_flag = (ewald_ev->param[3] > 0.);
    ewald_one_four = ewald_ev->param[1]-1.;
  }

  if (input->thread_id == 0)
    nblist_update(nblist, input->natoms, (double *)x, distance_data);
#ifdef WITH_THREAD
  barrier(eval->binfo+self->barrier_index, input->thread_id, input->nthreads);
#endif

  if (!(lj_flag || es_flag || ewald_flag))
    return;

  slicecounter = input->nslices-input->slice_id;
  for (ibox = 0; ibox < nblist->nboxes; ibox++) {
    nbbox *box1 = &nblist->boxes[ibox];
    int ineighbor;
    
    for (ineighbor = 0; ineighbor < nblist->nneighbors; ineighbor++) {
      int ix = nblist->neighbors[ineighbor][0]+box1->ix;
      int iy = nblist->neighbors[ineighbor][1]+box1->iy;
      int iz = nblist->neighbors[ineighbor][2]+box1->iz;
      nbbox *box2;
      int i, j;
      if (nblist->universe_spec->is_periodic) {
	if (ix < 0) ix += nblist->box_count[0];
	if (iy < 0) iy += nblist->box_count[1];
	if (iz < 0) iz += nblist->box_count[2];
	if (ix >= nblist->box_count[0]) ix -= nblist->box_count[0];
	if (iy >= nblist->box_count[1]) iy -= nblist->box_count[1];
	if (iz >= nblist->box_count[2]) iz -= nblist->box_count[2];
      }
      else if (ix < 0 || iy < 0 || iz < 0
	       || ix >= nblist->box_count[0]
	       || iy >= nblist->box_count[1]
	       || iz >= nblist->box_count[2])
	continue;
      jbox = ix + nblist->box_count[0]*(iy + nblist->box_count[1]*iz);
      if (jbox < ibox)
	continue;
      box2 = &nblist->boxes[jbox];
      if (ibox == jbox) {
	for (i = 0; i < box1->n; i++) {
	  int a1 = box1->atoms[i];
	  for (j = i+1; j < box2->n; j++) {
	    int a2 = box2->atoms[j];
	    if (--slicecounter == 0) {
	      slicecounter = input->nslices;
	      pair_term(1., 1., 1.);
#if THREAD_DEBUG
	      paircount++;
#endif
	    }
	  }
	}
      }
      else {
	for (i = 0; i < box1->n; i++) {
	  int a1 = box1->atoms[i];
	  for (j = 0; j < box2->n; j++) {
	    int a2 = box2->atoms[j];
	    if (--slicecounter == 0) {
	      slicecounter = input->nslices;
	      pair_term(1., 1., 1.);
#if THREAD_DEBUG
	      paircount++;
#endif
	    }
	  }
	}
      }
    }
  }
#if THREAD_DEBUG
    printf("Slice %d: %d nonbonded pairs\n", input->slice_id, paircount);
#endif

  es_inv_cutoff = 0.;
  ewald_inv_cutoff = 0.;
  erfc_cutoff = 0.;
  beta = 0.;
  for (k = 2*input->slice_id; k < n_ex; k += 2*input->nslices) {
    int a1 = excluded[k];
    int a2 = excluded[k+1];
    pair_term(-1., -1., -1.);
  }
  for (k = 2*input->slice_id; k < n_14; k += 2*input->nslices) {
    int a1 = one_four[k];
    int a2 = one_four[k+1];
    pair_term(lj_one_four, es_one_four, ewald_one_four);
  }

  energy->energy_terms[self->index] = lj_energy1 + lj_energy2;
  energy->energy_terms[self->index+1] = es_energy;
  energy->energy_terms[self->index+2] = ewald_energy;
  energy->energy_terms[self->virial_index] +=
    lj_virial1 + lj_virial2 + es_energy + ewald_energy;
}


/* Lennard-Jones and electrostatic interactions
 *
 * The real work is done in nonbonded_evaluator, these two functions
 * just return the energy values under the proper name.
 */

void
lennard_jones_evaluator(PyFFEnergyTermObject *self,
			PyFFEvaluatorObject *eval,
			energy_spec *input,
			energy_data *energy)
{
}

void
electrostatic_evaluator(PyFFEnergyTermObject *self,
			PyFFEvaluatorObject *eval,
			energy_spec *input,
			energy_data *energy)
{
  PyNonbondedListObject *nblist = (PyNonbondedListObject *)self->data[0];
  long *subset = (long *)((PyArrayObject *)nblist->atom_subset)->data;
  int n_sub = ((PyArrayObject *)nblist->atom_subset)->dimensions[0];
  double *charge = (double *)((PyArrayObject *)self->data[1])->data;
  double cutoff_sq = sqr(self->param[0]);
  double inv_cutoff = (self->param[0] == 0.) ? 0. : 1./self->param[0];
  double e = 0., v = 0.;
  int k;

  if (cutoff_sq > 0.) {
    int n = (n_sub == 0) ? input->natoms : n_sub;
    double sum = 0.;
    for (k = 0; k < n; k++) {
      int i = (n_sub == 0) ? k : subset[k];
      sum += sqr(charge[i]);
    }
    e -= 0.5*inv_cutoff*sum*electrostatic_energy_factor;
    v -= 0.5*inv_cutoff*sum*electrostatic_energy_factor;
  }
  energy->energy_terms[self->index] = e;
  energy->energy_terms[self->virial_index] += v;
}


#ifdef WITH_DPMTA
/* Electrostatic interactions via multipoles */

#ifdef GRADIENTFN
#define pair_es_energy(qiqj, distance_vector_function) \
{ \
  vector3 rij; \
  double r_sq; \
  double deriv = 0.; \
  double deriv2 = 0.; \
  distance_vector_function(rij, x[j], x[i], distance_data); \
  r_sq = vector_length_sq(rij); \
  if ((qiqj) != 0.) { \
    double r = sqrt(r_sq); \
    double d = -(qiqj)/(r*r_sq); \
    e += (qiqj)/r; \
    v += (qiqj)/r; \
    deriv += d; \
    deriv2 += 2.*(qiqj)/(r*r_sq); \
  } \
  if (energy->gradients != NULL) { \
    vector3 grad; \
    grad[0] = deriv*rij[0]; \
    grad[1] = deriv*rij[1]; \
    grad[2] = deriv*rij[2]; \
    if (energy->gradient_fn != NULL) { \
      (*energy->gradient_fn)(energy, i, grad); \
      vector_changesign(grad); \
      (*energy->gradient_fn)(energy, j, grad); \
    } \
    else { \
      vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data; \
      f[i][0] += grad[0]; \
      f[i][1] += grad[1]; \
      f[i][2] += grad[2]; \
      f[j][0] -= grad[0]; \
      f[j][1] -= grad[1]; \
      f[j][2] -= grad[2]; \
    } \
  } \
  if (energy->force_constants != NULL) \
    add_pair_fc(energy, i, j, rij, r_sq, deriv, deriv2); \
}
#else
#define pair_es_energy(qiqj, distance_vector_function) \
{ \
  vector3 rij; \
  double r_sq; \
  double deriv = 0.; \
  double deriv2 = 0.; \
  distance_vector_function(rij, x[j], x[i], distance_data); \
  r_sq = vector_length_sq(rij); \
  if ((qiqj) != 0.) { \
    double r = sqrt(r_sq); \
    double d = -(qiqj)/(r*r_sq); \
    e += (qiqj)/r; \
    v += (qiqj)/r; \
    deriv += d; \
    deriv2 += 2.*(qiqj)/(r*r_sq); \
  } \
  if (energy->gradients != NULL) { \
    vector3 grad; \
    grad[0] = deriv*rij[0]; \
    grad[1] = deriv*rij[1]; \
    grad[2] = deriv*rij[2]; \
    { \
      vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data; \
      f[i][0] += grad[0]; \
      f[i][1] += grad[1]; \
      f[i][2] += grad[2]; \
      f[j][0] -= grad[0]; \
      f[j][1] -= grad[1]; \
      f[j][2] -= grad[2]; \
    } \
  } \
  if (energy->force_constants != NULL) \
    add_pair_fc(energy, i, j, rij, r_sq, deriv, deriv2); \
}
#endif

void
es_mp_evaluator(PyFFEnergyTermObject *self,
		PyFFEvaluatorObject *eval,
		energy_spec *input,
		energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  distance_fn *d_fn = self->universe_spec->distance_function;
  double *distance_data = self->universe_spec->geometry_data;
  PmtaInitData *initdata;
  static PmtaInitData *last_initdata = NULL;
  PmtaParticle *particle;
  PmtaPartInfo *force;
  PmtaVector virial;
  double potential;
  PyObject *pair_lists = self->data[0];
  PyArrayObject *array_ex = 
                  (PyArrayObject *)PyList_GetItem(pair_lists, (Py_ssize_t)0);
  PyArrayObject *array_14 =
                  (PyArrayObject *)PyList_GetItem(pair_lists, (Py_ssize_t)1);
  PyArrayObject *array_subset =
                  (PyArrayObject *)PyList_GetItem(pair_lists, (Py_ssize_t)2);
  long *excluded = (long *)array_ex->data;
  long n_ex = 2*array_ex->dimensions[0];
  long *one_four = (long *)array_14->data;
  long n_14 = 2*array_14->dimensions[0];
  long *subset = (long *)array_subset->data;
  long n_sub = array_subset->dimensions[0];
  double *charge = (double *)((PyArrayObject *)self->data[1])->data;
  double es_one_four_factor = self->param[0];
  double e, v;
  int k, n;

  n = (n_sub == 0) ? input->natoms : n_sub;
  initdata = (PmtaInitData *)self->scratch;
  particle = (PmtaParticle *)(initdata+1);
  force = (PmtaPartInfo *)(particle+n);
  if (initdata != last_initdata) {
    PMTAinit(initdata, NULL);
    last_initdata = initdata;
  }
  for (k = 0; k < n; k++) {
    int i = (n_sub == 0) ? k : subset[k];
    particle[k].p.x = x[i][0];
    particle[k].p.y = x[i][1];
    particle[k].p.z = x[i][2];
    particle[k].q = charge[i];
  }
  if (self->param[1] == 0.) {
    PmtaVector v1, v2, v3, center;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = ymin = zmin = undefined;
    xmax = ymax = zmax = -undefined;
    for (k = 0; k < n; k++) {
      int i = (n_sub == 0) ? k : subset[k];
      if (x[i][0] < xmin) xmin = x[i][0];
      if (x[i][1] < ymin) ymin = x[i][1];
      if (x[i][2] < zmin) zmin = x[i][2];
      if (x[i][0] > xmax) xmax = x[i][0];
      if (x[i][1] > ymax) ymax = x[i][1];
      if (x[i][2] > zmax) zmax = x[i][2];
    }
    v1.x = xmax-xmin; v1.y = 0.; v1.z = 0.;
    v2.x = 0.; v2.y = ymax-ymin; v2.z = 0.;
    v3.x = 0.; v3.y = 0.; v3.z = zmax-zmin;
    center.x = 0.5*(xmin+xmax);
    center.y = 0.5*(ymin+ymax);
    center.z = 0.5*(zmin+zmax);
    PMTAresize(&v1, &v2, &v3, &center);
  }
  else {
    int change = 0;
    if (distance_data[0] != initdata->v1.x) {
      initdata->v1.x = distance_data[0];
      change = 1;
    }
    if (distance_data[1] != initdata->v2.y) {
      initdata->v2.y = distance_data[1];
      change = 1;
    }
    if (distance_data[2] != initdata->v3.z) {
      initdata->v3.z = distance_data[2];
      change = 1;
    }
    if (change)
      PMTAresize(&initdata->v1, &initdata->v2, &initdata->v3,
		 &initdata->cellctr);
  }
  PMTAforce(n, particle, force, NULL);
  PMTAvirial(&potential, &virial, NULL, NULL);
  for (k = 0; k < n; k++) {
    int i = (n_sub == 0) ? k : subset[k];
    vector3 grad;
    grad[0] = -electrostatic_energy_factor*force[k].f.x;
    grad[1] = -electrostatic_energy_factor*force[k].f.y;
    grad[2] = -electrostatic_energy_factor*force[k].f.z;
    if (energy->gradients != NULL) {
#ifdef GRADIENTFN
      if (energy->gradient_fn != NULL)
	(*energy->gradient_fn)(energy, i, grad);
      else
#endif
      {
	vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data;
	f[i][0] += grad[0];
	f[i][1] += grad[1];
	f[i][2] += grad[2];
      }
    }
  }
  e = potential*electrostatic_energy_factor;
  v = -electrostatic_energy_factor*(virial.x+virial.y+virial.z);
  if (energy->force_constants != NULL) {
    PyErr_SetString(PyExc_ValueError,
		    "no second derivatives in multipole evaluator");
    energy->error = 1;
  }
  if (d_fn == distance_vector_pointer) {
#   define distance_vector_function distance_vector_1
#   include "nonbonded1.i"
#   undef distance_vector_function
  }
  else if (d_fn == orthorhombic_distance_vector_pointer) {
#   define distance_vector_function distance_vector_2
#   include "nonbonded1.i"
#   undef distance_vector_function
  }
  else if (d_fn == parallelepipedic_distance_vector_pointer) {
#   define distance_vector_function distance_vector_3
#   include "nonbonded1.i"
#   undef distance_vector_function
  }
  else {
#   define distance_vector_function (*d_fn)
#   include "nonbonded1.i"
#   undef distance_vector_function
  }
  energy->energy_terms[self->index] = e;
  energy->energy_terms[self->virial_index] += v;
}
#endif
