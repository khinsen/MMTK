/* Ewald method for electrostatic interactions
 *
 * Written by Konrad Hinsen
 * last revision: 2002-12-10
 */


#define NO_IMPORT
#define _FORCEFIELD_MODULE
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_MMTKFF_API

#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"

#define DEBUG 0

/*
 * Include erfc code, unless user wants to use an erfc from libm.
 */
#ifndef LIBM_HAS_ERFC
#include "polevl.c"
#include "ndtr.c"
#endif
#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif

/*
 * Plain Ewald sum, for orthorhombic systems only.
 */
typedef struct {
   double real;
   double imag;
} complex;

static complex c_1 = {1., 0.};

#define c_mult_real(a, b) ((a).real*(b).real - (a).imag*(b).imag)
#define c_mult_imag(a, b) ((a).real*(b).imag + (a).imag*(b).real)

#if 0
static inline complex
c_mult(complex a, complex b)
{
  complex r;
  r.real = a.real*b.real - a.imag*b.imag;
  r.imag = a.real*b.imag + a.imag*b.real;
  return r;
}
#endif

static double
c_abssq(complex z)
{
  return z.real*z.real + z.imag*z.imag;
}

static double
reciprocal_sum(energy_spec *input, energy_data *energy,
	       double volume, double *charge, double beta,
	       long *kmax, double cutoff_sq,
	       box_fn *box_transformation_fn, double *universe_data,
	       void *scratch, PyFFEvaluatorObject *eval,
	       PyFFEnergyTermObject *term)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  vector3 *xb = (vector3 *)scratch;
  complex *eikx = (complex *)(xb + input->natoms);
  complex *eiky = eikx + input->natoms*(kmax[0]+1);
  complex *eikz = eiky + input->natoms*(2*kmax[1]+1);
  complex *eikr = eikz + input->natoms*(2*kmax[2]+1);
  int *nkvect = (int *)(eikr+input->natoms);
  int *kx = nkvect + 1;
  int *ky = kx + *nkvect;
  int *kz = ky + *nkvect;
  vector3 r1, r2, r3;
  int nx, ny, nz, negy, negz, xk, yk, zk;
  int atoms_per_slice, first_atom, last_atom;
  int vectors_per_slice, first_vector, last_vector;
  double e;
  int i, nk;

  if (input->thread_id == 0) {
    (*box_transformation_fn)(x, xb, input->natoms, universe_data, 1);
    for (i = 0; i < input->natoms; i++) {
      xb[i][0] *= 2.*M_PI;
      xb[i][1] *= 2.*M_PI;
      xb[i][2] *= 2.*M_PI;
    }
  }
  else {
    eikr = (complex *)malloc(2*input->natoms*sizeof(double));
    if (eikr == NULL) {
      energy->error = 1;
      return 0.;
    }
  }
#ifdef WITH_THREAD
  barrier(eval->binfo+term->barrier_index, input->thread_id, input->nthreads);
#endif

  r1[0] = 2.*M_PI;  r1[1] = r1[2] = 0.;
  r2[1] = 2.*M_PI;  r2[2] = r2[0] = 0.;
  r3[2] = 2.*M_PI;  r3[0] = r3[1] = 0.;
  (*box_transformation_fn)(&r1, &r1, 1, universe_data, 1);
  (*box_transformation_fn)(&r2, &r2, 1, universe_data, 1);
  (*box_transformation_fn)(&r3, &r3, 1, universe_data, 1);

  nx = kmax[0]+1;
  ny = 2*kmax[1]+1;
  nz = 2*kmax[2]+1;
  negy = kmax[1];
  negz = kmax[2];

  atoms_per_slice = (input->natoms+input->nslices-1)/input->nslices;
  first_atom = input->slice_id*atoms_per_slice;
  last_atom = (input->slice_id+1)*atoms_per_slice;
  if (last_atom > input->natoms)
    last_atom = input->natoms;

  for (i = first_atom; i < last_atom; i++) {
    eikx[nx*i] = c_1;
    eiky[ny*i] = c_1;
    eikz[nz*i] = c_1;
    eikx[nx*i+1].real = cos(xb[i][0]);
    eikx[nx*i+1].imag = sin(xb[i][0]);
    eiky[ny*i+1].real = cos(xb[i][1]);
    eiky[ny*i+1].imag = sin(xb[i][1]);
    eikz[nz*i+1].real = cos(xb[i][2]);
    eikz[nz*i+1].imag = sin(xb[i][2]);
    eiky[ny*i+negy+1].real = eiky[ny*i+1].real;
    eiky[ny*i+negy+1].imag = -eiky[ny*i+1].imag;
    eikz[nz*i+negz+1].real = eikz[nz*i+1].real;
    eikz[nz*i+negz+1].imag = -eikz[nz*i+1].imag;
    for (xk = 2; xk <= kmax[0]; xk++) {
      eikx[nx*i+xk].real = c_mult_real(eikx[nx*i+xk-1], eikx[nx*i+1]);
      eikx[nx*i+xk].imag = c_mult_imag(eikx[nx*i+xk-1], eikx[nx*i+1]);
    }
    for (yk = 2; yk <= kmax[1]; yk++) {
      eiky[ny*i+yk].real = c_mult_real(eiky[ny*i+yk-1], eiky[ny*i+1]);
      eiky[ny*i+yk].imag = c_mult_imag(eiky[ny*i+yk-1], eiky[ny*i+1]);
      eiky[ny*i+negy+yk].real = eiky[ny*i+yk].real;
      eiky[ny*i+negy+yk].imag = -eiky[ny*i+yk].imag;
    }
    for (zk = 2; zk <= kmax[2]; zk++) {
      eikz[nz*i+zk].real = c_mult_real(eikz[nz*i+zk-1], eikz[nz*i+1]);
      eikz[nz*i+zk].imag = c_mult_imag(eikz[nz*i+zk-1], eikz[nz*i+1]);
      eikz[nz*i+negz+zk].real = eikz[nz*i+zk].real;
      eikz[nz*i+negz+zk].imag = -eikz[nz*i+zk].imag;
    }
  }
#ifdef WITH_THREAD
  barrier(eval->binfo+term->barrier_index+1,
	  input->thread_id, input->nthreads);
#endif

  vectors_per_slice = (*nkvect+input->nslices-1)/input->nslices;
  first_vector = input->slice_id*vectors_per_slice;
  last_vector = (input->slice_id+1)*vectors_per_slice;
  if (last_vector > *nkvect)
    last_vector = *nkvect;

  e = 0.;
  for (nk = first_vector; nk < last_vector; nk++) {
    int xk = kx[nk];
    int yk = ky[nk];
    int zk = kz[nk];
    int iyk = (yk < 0) ? negy-yk : yk;
    int izk = (zk < 0) ? negz-zk : zk;
    vector3 k = {0., 0., 0.};
    double ksq, expfactor;
    complex sum;

    vector_add(k, r1, xk);
    vector_add(k, r2, yk);
    vector_add(k, r3, zk);
    ksq = vector_length_sq(k);
    expfactor = exp(-ksq/(4.*beta*beta))/ksq;
    if (xk != 0)
      expfactor *= 2;

    sum.real = sum.imag = 0.;
    for (i = 0; i < input->natoms; i++) {
      complex temp;
      temp.real = c_mult_real(eikx[nx*i+xk], eiky[ny*i+iyk]);
      temp.imag = c_mult_imag(eikx[nx*i+xk], eiky[ny*i+iyk]);
      eikr[i].real = c_mult_real(temp, eikz[nz*i+izk]);
      eikr[i].imag = c_mult_imag(temp, eikz[nz*i+izk]);
      sum.real += charge[i]*eikr[i].real;
      sum.imag += charge[i]*eikr[i].imag;
    }
    e += c_abssq(sum)*expfactor;
    if (energy->gradients != NULL) {
      for (i = 0; i < input->natoms; i++) {
	vector3 grad;
	double a = 4.*M_PI*expfactor*electrostatic_energy_factor *
	  charge[i]*(sum.imag*eikr[i].real-sum.real*eikr[i].imag)/volume;
	grad[0] = a*k[0];
	grad[1] = a*k[1];
	grad[2] = a*k[2];
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
    if (energy->force_constants != NULL) {
      tensor3 fij;
      int j;
      tensor_product(fij, k, k,
		     4.*M_PI*expfactor*electrostatic_energy_factor/volume);
#if 0
      for (i = 0; i < input->natoms; i++)
	for (j = i; j < input->natoms; j++) {
	  double f = charge[i]*charge[j];
	}
      if (energy->fc_fn != NULL) {
      }
      else {
      }
#endif
      PyErr_SetString(PyExc_ValueError, "not yet implemented");
      energy->error = 1;
    }
  }
  if (input->thread_id != 0)
    free(eikr);
  return 2.*M_PI*e*electrostatic_energy_factor/volume;
}

/*
 * Count k vectors for memory allocation and initialize scratch area.
 */
int
init_kvectors(box_fn *box_transformation_fn, double *universe_data, int natoms,
	      long *kmax, double cutoff_sq, void *scratch, int nvect)
{
  vector3 *xb = (vector3 *)scratch;  /* Size natoms */
  complex *eikx = (complex *)(xb + natoms);  /* Size (kmax[0]+1)*natoms */
  complex *eiky = eikx + natoms*(kmax[0]+1);  /* Size (2*kmax[1]+1)*natoms */
  complex *eikz = eiky + natoms*(2*kmax[1]+1);  /* Size (2*kmax[2]+1)*natoms */
  complex *eikr = eikz + natoms*(2*kmax[2]+1);  /* Size natoms */
  int *nkvect = (int *)(eikr+natoms);  /* Size 1 */
  int *kx = nkvect + 1;  /* Size nvect */
  int *ky = kx + nvect;  /* Size nvect */
  int *kz = ky + nvect;  /* Size nvect */

  vector3 r1, r2, r3;
  int x, y, z;
  int nk = 0;

  r1[0] = 2.*M_PI;  r1[1] = r1[2] = 0.;
  r2[1] = 2.*M_PI;  r2[2] = r2[0] = 0.;
  r3[2] = 2.*M_PI;  r3[0] = r3[1] = 0.;
  (*box_transformation_fn)(&r1, &r1, 1, universe_data, 1);
  (*box_transformation_fn)(&r2, &r2, 1, universe_data, 1);
  (*box_transformation_fn)(&r3, &r3, 1, universe_data, 1);

  if (scratch != NULL)
    *nkvect = nvect;
  for (x = 0; x <= kmax[0]; x++) {
    for (y = -kmax[1]; y <= kmax[1]; y++) {
      for (z = -kmax[2]; z <= kmax[2]; z++) {
	vector3 k = {0., 0., 0.};
	double ksq;
	vector_add(k, r1, x);
	vector_add(k, r2, y);
	vector_add(k, r3, z);
	ksq = vector_length_sq(k);
	if (ksq < cutoff_sq && ksq > 0) {
	  if (scratch != NULL) {
	    kx[nk] = x; ky[nk] = y; kz[nk] = z;
	  }
	  nk++;
	}
      }
    }
  }
#if 0
  printf("%d k vectors\n", nk);
#endif
  return nk;
}

/* The Ewald evaluator does not calculate the real-space sum,
   which is evaluated in nonbonded_evaluator! */
void
es_ewald_evaluator(PyFFEnergyTermObject *self,
		   PyFFEvaluatorObject *eval,
		   energy_spec *input,
		   energy_data *energy)
{
  box_fn *box_transformation_fn = self->universe_spec->box_function;
  volume_fn *v_fn = self->universe_spec->volume_function;
  double *universe_data = self->universe_spec->geometry_data;
  double volume = (*v_fn)(1., universe_data);

  double *charge = (double *)((PyArrayObject *)self->data[1])->data;
  long *kmax = (long *)((PyArrayObject *)self->data[2])->data;
  double inv_cutoff = (self->param[0] == 0.) ? 0. : 1./self->param[0];
  double beta = self->param[2];
  double reciprocal_cutoff_sq = self->param[3];
  double charge_sum;
  double e = 0.;
  int k;

  if (reciprocal_cutoff_sq > 0.)
    inv_cutoff = 0.;

  if (input->slice_id == 0) {
    charge_sum = 0.;
    for (k = 0; k < input->natoms; k++)
      charge_sum += charge[k]*charge[k];
    e -= charge_sum*electrostatic_energy_factor * 
         (beta/sqrt(M_PI)+0.5*inv_cutoff*erfc(beta*self->param[0]));
  }
  energy->energy_terms[self->index] = e;

  if (reciprocal_cutoff_sq > 0.)
    energy->energy_terms[self->index+1] =
      reciprocal_sum(input, energy, volume, charge, beta,
		     kmax, reciprocal_cutoff_sq,
		     box_transformation_fn, universe_data,
		     self->scratch, eval, self);

  energy->energy_terms[self->virial_index] +=
    energy->energy_terms[self->index] + energy->energy_terms[self->index+1];
}
