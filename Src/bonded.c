/* Low-level force field calculations: bonded interactions
 *
 * Written by Konrad Hinsen
 * last revision: 2009-4-22
 */

#define NO_IMPORT
#define _FORCEFIELD_MODULE
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_MMTKFF_API

#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"


/* Harmonic bond potential */

void
harmonic_bond_evaluator(PyFFEnergyTermObject *self,
			PyFFEvaluatorObject *eval,
			energy_spec *input,
			energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  long *index = (long *)((PyArrayObject *)self->data[0])->data;
  double *param = (double *)((PyArrayObject *)self->data[1])->data;

  int nterms = (self->n+input->nslices-1)/(input->nslices);
  int term = input->slice_id*nterms;
  int last_term = (input->slice_id+1)*nterms;
  double e = 0., v = 0.;

  if (last_term > self->n)
    last_term = self->n;
  index += 2*term;
  param += 2*term;

  /* Loop over bond terms */
  while (term++ < last_term) {
    /*
     * E = k (r-r0)^2
     * r = | R_i-R_j |
     *
     * index[0], index[1]: particle indices i, j
     * param[0]: r0
     * param[1]: k
     */
    int i = index[0];
    int j = index[1];
    vector3 rij;
    double lrij, dr;
    self->universe_spec->distance_function(rij, x[j], x[i],
					   self->universe_spec->geometry_data);
    lrij = vector_length(rij);
    dr = lrij-param[0];
    e += param[1]*sqr(dr);
    v += -2.*param[1]*dr*lrij;
    if (energy->gradients != NULL) {
      double deriv = 2.*param[1]*dr/lrij;
      vector3 grad;
      grad[0] = deriv*rij[0];
      grad[1] = deriv*rij[1];
      grad[2] = deriv*rij[2];
#ifdef GRADIENTFN
      if (energy->gradient_fn != NULL) {
	(*energy->gradient_fn)(energy, i, grad);
	vector_changesign(grad);
	(*energy->gradient_fn)(energy, j, grad);
      }
      else
#endif
      {
	vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data;
	f[i][0] += grad[0];
	f[i][1] += grad[1];
	f[i][2] += grad[2];
	f[j][0] -= grad[0];
	f[j][1] -= grad[1];
	f[j][2] -= grad[2];
      }
    }
    if (energy->force_constants != NULL) {
      double f1 = 2.*param[1]*dr/lrij;
      double f2 = 2.*param[1];
      add_pair_fc(energy, i, j, rij, sqr(lrij), f1, f2);
    }
    index += 2;
    param += 2;
  }
  energy->energy_terms[self->index] = e;
  energy->energy_terms[self->virial_index] += v;
}

/* Harmonic angle potential */

void
harmonic_angle_evaluator(PyFFEnergyTermObject *self,
			 PyFFEvaluatorObject *eval,
			 energy_spec *input,
			 energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;
  long *index = (long *)((PyArrayObject *)self->data[0])->data;
  double *param = (double *)((PyArrayObject *)self->data[1])->data;

  int nterms = (self->n+input->nslices-1)/(input->nslices);
  int term = input->slice_id*nterms;
  int last_term = (input->slice_id+1)*nterms;
  double e = 0.;

  if (last_term > self->n)
    last_term = self->n;
  index += 3*term;
  param += 2*term;

  /* Loop over bond angle terms */
  while (term++ < last_term) {
    /*
     * E = k (theta-theta0)^2
     * cos theta = (R_i-R_j)*(R_k-R_j)
     *
     * index[1]: index of central particle j
     * index[0], index[2]: indices of outer particles i, k
     * param[0]: theta0 (radians)
     * param[1]: k
     */
    int i = index[0];
    int j = index[1];
    int k = index[2];
    vector3 rij, rkj;
    double lrij, lrkj;
    double cos_theta, sin_theta, theta, dtheta;
    self->universe_spec->distance_function(rij, x[j], x[i],
					   self->universe_spec->geometry_data);
    lrij = vector_length(rij);
    vector_scale(rij, 1./lrij);
    self->universe_spec->distance_function(rkj, x[j], x[k],
					   self->universe_spec->geometry_data);
    lrkj = vector_length(rkj);
    vector_scale(rkj, 1./lrkj);
    cos_theta = dot(rij, rkj);
    if (cos_theta > 1.) cos_theta = 1.;
    if (cos_theta < -1.) cos_theta = -1.;
    sin_theta = sqrt(1.-sqr(cos_theta));
    theta = acos(cos_theta);
    dtheta = (theta-param[0]);
    e += param[1]*sqr(dtheta);
    if (energy->gradients != NULL || energy->force_constants != NULL) {
      /* First derivative of angle potential */
      double deriv = -2.*param[1]*dtheta/sin_theta;
      vector3 di, dk, dj;
      di[0] = (rkj[0]-cos_theta*rij[0])/lrij;
      di[1] = (rkj[1]-cos_theta*rij[1])/lrij;
      di[2] = (rkj[2]-cos_theta*rij[2])/lrij;
      dk[0] = (rij[0]-cos_theta*rkj[0])/lrkj;
      dk[1] = (rij[1]-cos_theta*rkj[1])/lrkj;
      dk[2] = (rij[2]-cos_theta*rkj[2])/lrkj;
      dj[0] = -di[0]-dk[0];
      dj[1] = -di[1]-dk[1];
      dj[2] = -di[2]-dk[2];
      if (energy->gradients != NULL) {
#ifdef GRADIENTFN
	if (energy->gradient_fn != NULL) {
	  vector3 grad;
	  vector_copy(grad, di);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, i, grad);
	  vector_copy(grad, dj);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, j, grad);
	  vector_copy(grad, dk);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, k, grad);
	}
	else
#endif
        {
	  vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data;
	  f[i][0] += deriv*di[0];
	  f[i][1] += deriv*di[1];
	  f[i][2] += deriv*di[2];
	  f[j][0] += deriv*dj[0];
	  f[j][1] += deriv*dj[1];
	  f[j][2] += deriv*dj[2];
	  f[k][0] += deriv*dk[0];
	  f[k][1] += deriv*dk[1];
	  f[k][2] += deriv*dk[2];
	}
      }
      if (energy->force_constants != NULL) {
	/* Second derivative of angle potential */
	double deriv2 = 2.*param[1] *
	      (1.-(cos_theta/sin_theta)*dtheta)/sqr(sin_theta);
	double min_r = (lrij < lrkj) ? sqr(lrij) : sqr(lrkj);
	int swapij = (i > j);
	int swapik = (i > k);
	int swapjk = (j > k);
	if (energy->fc_fn != NULL) {
	  int l, m;
	  tensor3 fcii, fcjj, fckk, fcij, fcik, fcjk;
	  for (l = 0; l < 3; l++)
	    for (m = 0; m < 3; m++) {
	      double a, b, ab, ba;
	      a = 3.*cos_theta*rij[l]*rij[m]-rij[l]*rkj[m]-rkj[l]*rij[m];
	      b = 3.*cos_theta*rkj[l]*rkj[m]-rkj[l]*rij[m]-rij[l]*rkj[m];
	      ab = rij[l]*rkj[m]*cos_theta-rij[l]*rij[m]-rkj[l]*rkj[m];
	      ba = rkj[l]*rij[m]*cos_theta-rkj[l]*rkj[m]-rij[l]*rij[m];
	      if (l == m) {
		a -= cos_theta;
		b -= cos_theta;
		ab += 1.;
		ba += 1.;
	      }
	      a *= deriv/sqr(lrij);
	      b *= deriv/sqr(lrkj);
	      ab *= deriv/(lrij*lrkj);
	      ba *= deriv/(lrij*lrkj);
	      fcii[l][m] = deriv2*di[l]*di[m] + a;
	      fcjj[l][m] = deriv2*dj[l]*dj[m] + a + b + ab + ba;
	      fckk[l][m] = deriv2*dk[l]*dk[m] + b;
	      fcij[l][m] = deriv2*di[l]*dj[m] - a - ab;
	      fcik[l][m] = deriv2*di[l]*dk[m] + ab;
	      fcjk[l][m] = deriv2*dj[l]*dk[m] - ab - b;
	    }
	  (*energy->fc_fn)(energy, i, i, fcii, min_r);
	  (*energy->fc_fn)(energy, j, j, fcjj, min_r);
	  (*energy->fc_fn)(energy, k, k, fckk, min_r);
	  if (swapij) {
	    tensor_transpose(fcij);
	    (*energy->fc_fn)(energy, j, i, fcij, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, i, j, fcij, min_r);
	  if (swapik) {
	    tensor_transpose(fcik);
	    (*energy->fc_fn)(energy, k, i, fcik, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, i, k, fcik, min_r);
	  if (swapjk) {
	    tensor_transpose(fcjk);
	    (*energy->fc_fn)(energy, k, j, fcjk, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, j, k, fcjk, min_r);
	}
	else {
	  double *fc_data = 
	    (double *)((PyArrayObject *)energy->force_constants)->data;
	  double *fcii = fc_data + 9*input->natoms*i+3*i;
	  double *fcjj = fc_data + 9*input->natoms*j+3*j;
	  double *fckk = fc_data + 9*input->natoms*k+3*k;
	  double *fcij, *fcik, *fcjk;
	  int l, m;
	  if (swapij)
	    fcij = fc_data + 9*input->natoms*j+3*i;
	  else
	    fcij = fc_data + 9*input->natoms*i+3*j;
	  if (swapik)
	    fcik = fc_data + 9*input->natoms*k+3*i;
	  else
	    fcik = fc_data + 9*input->natoms*i+3*k;
	  if (swapjk)
	    fcjk = fc_data + 9*input->natoms*k+3*j;	
	  else
	    fcjk = fc_data + 9*input->natoms*j+3*k;
	  for (l = 0; l < 3; l++)
	    for (m = 0; m < 3; m++) {
	      int o = 3*input->natoms*l + m;
	      double a, b, ab, ba;
	      fcii[o] += deriv2*di[l]*di[m];
	      fcjj[o] += deriv2*dj[l]*dj[m];
	      fckk[o] += deriv2*dk[l]*dk[m];
	      if (swapij)
		fcij[o] += deriv2*dj[l]*di[m];
	      else
		fcij[o] += deriv2*di[l]*dj[m];
	      if (swapik)
		fcik[o] += deriv2*dk[l]*di[m];
	      else
		fcik[o] += deriv2*di[l]*dk[m];
	      if (swapjk)
		fcjk[o] += deriv2*dk[l]*dj[m];
	      else
		fcjk[o] += deriv2*dj[l]*dk[m];
	      a = 3.*cos_theta*rij[l]*rij[m]-rij[l]*rkj[m]-rkj[l]*rij[m];
	      b = 3.*cos_theta*rkj[l]*rkj[m]-rkj[l]*rij[m]-rij[l]*rkj[m];
	      ab = rij[l]*rkj[m]*cos_theta-rij[l]*rij[m]-rkj[l]*rkj[m];
	      ba = rkj[l]*rij[m]*cos_theta-rkj[l]*rkj[m]-rij[l]*rij[m];
	      if (l == m) {
		a -= cos_theta;
		b -= cos_theta;
		ab += 1.;
		ba += 1.;
	      }
	      a *= deriv/sqr(lrij);
	      b *= deriv/sqr(lrkj);
	      ab *= deriv/(lrij*lrkj);
	      ba *= deriv/(lrij*lrkj);
	      fcii[o] += a;
	      fcjj[o] += a + b + ab + ba;
	      fckk[o] += b;
	      if (swapij)
		fcij[o] -= a + ba;
	      else
		fcij[o] -= a + ab;
	      if (swapik)
		fcik[o] += ba;
	      else
		fcik[o] += ab;
	      if (swapjk)
		fcjk[o] -= ba + b;
	      else
		fcjk[o] -= ab + b;
	    }
	}
      }
    }
    index += 3;
    param += 2;
  }
  energy->energy_terms[self->index] = e;
}

/* Cosine dihedral potential */

static void
add_fc_tensor(double *fc, int n, int swap, tensor3 t, double f)
{
  int i, j;
  if (swap) {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) {
	int o = 3*n*i+j;
	fc[o] += f*t[j][i];
      }
  }
  else {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) {
	int o = 3*n*i+j;
	fc[o] += f*t[i][j];
      }
  }
}

void
cosine_dihedral_evaluator(PyFFEnergyTermObject *self,
			  PyFFEvaluatorObject *eval,
			  energy_spec *input,
			  energy_data *energy)
{
  vector3 *x = (vector3 *)input->coordinates->data;

  long *index = (long *)((PyArrayObject *)self->data[0])->data;
  double *param = (double *)((PyArrayObject *)self->data[1])->data;

  int nterms = (self->n+input->nslices-1)/(input->nslices);
  int term = input->slice_id*nterms;
  int last_term = (input->slice_id+1)*nterms;
  double e = 0.;

  if (last_term > self->n)
    last_term = self->n;
  index += 4*term;
  param += 4*term;

  /* Loop over dihedral angle terms */
  while (term++ < last_term) {
    /*
     * E = V [1 + cos(n phi-gamma)]
     * phi = angle between plane1 and plane2 using the IUPAC sign convention
     * plane1 defined by R_i, R_j, R_k
     * plane2 defined by R_j, R_k, R_l
     *
     * index[1], index[1]: indices of the axis particles j, k
     * index[0], index[3]: indices of outer particles i, l
     * param[0]: n = int(param[0])
     * param[1]: cos gamma
     * param[2]: sin gamma
     * param[3]: V
     */
    int i = index[0];
    int j = index[1];
    int k = index[2];
    int l = index[3];
    int n;
    vector3 rij, rkj, rlk, rkj_cross_rkl, rij_cross_rkj, r, s;
    double lrij, lrkj, lrlk, lm, ln, lr, ls;
    double dot_rij_rkj, dot_rlk_rkj;
    double cos_phi, sqr_cos_phi, cos_n_phi;
    double sin_phi, sin_n_phi_ratio;
    double phi, dphi, sign_phi, cos_phase, sin_phase;
    self->universe_spec->distance_function(rij, x[j], x[i],
					   self->universe_spec->geometry_data);
    lrij = vector_length(rij);
    self->universe_spec->distance_function(rkj, x[j], x[k],
					   self->universe_spec->geometry_data);
    lrkj = vector_length(rkj);
    self->universe_spec->distance_function(rlk, x[k], x[l],
					   self->universe_spec->geometry_data);
    lrlk = vector_length(rlk);
    cross(rkj_cross_rkl, rlk, rkj);
    ln = vector_length(rkj_cross_rkl);
    cross(rij_cross_rkj, rij, rkj);
    lm = vector_length(rij_cross_rkj);
    sign_phi = 1.;
    if (dot(rij, rkj_cross_rkl) < 0.) sign_phi = -1.;
    vector_scale(rkj, 1./lrkj);
    dot_rij_rkj = dot(rij, rkj);
    r[0] = rij[0]-dot_rij_rkj*rkj[0];
    r[1] = rij[1]-dot_rij_rkj*rkj[1];
    r[2] = rij[2]-dot_rij_rkj*rkj[2];
    lr = vector_length(r);
    vector_scale(r, 1./lr);
    dot_rlk_rkj = dot(rlk, rkj);
    s[0] = rlk[0]-dot_rlk_rkj*rkj[0];
    s[1] = rlk[1]-dot_rlk_rkj*rkj[1];
    s[2] = rlk[2]-dot_rlk_rkj*rkj[2];
    ls = vector_length(s);
    vector_scale(s, 1./ls);
    cos_phi = dot(r, s);
    if (cos_phi > 1.) cos_phi = 1.;
    if (cos_phi < -1.) cos_phi = -1.;
    n = (int)param[0];
    if (n == 0) {
      phi = acos(cos_phi)*sign_phi;
      dphi = phi-param[1];
      dphi = fmod(dphi+2.*M_PI, 2.*M_PI);
      if (dphi > M_PI)
	dphi -= 2*M_PI;
      e += param[3]*sqr(dphi);
    }
    else {
      cos_phase = param[1];
      sin_phase = param[2];
      sqr_cos_phi = sqr(cos_phi);
      sin_phi = 2.;
      switch (n) {
      case 1:
	cos_n_phi = cos_phi;
	sin_n_phi_ratio = 1.;
	break;
      case 2:
	cos_n_phi = 2.*sqr_cos_phi-1.;
	sin_n_phi_ratio = 2.*cos_phi;
	break;
      case 3:
	cos_n_phi = (4.*sqr_cos_phi-3.)*cos_phi;
	sin_n_phi_ratio = 4.*sqr_cos_phi - 1.;
	break;
      case 4:
	cos_n_phi = 8.*(sqr_cos_phi-1.)*sqr_cos_phi+1.;
	sin_n_phi_ratio = 4.*(2.*sqr_cos_phi-1.)*cos_phi;
	break;
#if 0
	/* n=5 and n=6 don't occur in the Amber force field and have
	   therefore not been tested. Use this section at your own risk. */
      case 5:
	cos_n_phi = ((16.*sqr_cos_phi-20.)*sqr_cos_phi+5.)*cos_phi;
	sin_n_phi_ratio =4.*(4.*sqr_cos_phi-3.)*sqr_cos_phi+1.;
	break;
      case 6:
	cos_n_phi = ((32.*sqr_cos_phi-48.)*sqr_cos_phi+18.)*sqr_cos_phi-1.;
	sin_n_phi_ratio = (32.*(sqr_cos_phi-1.)*sqr_cos_phi+6.)*cos_phi;
	break;
#endif
      default:
	phi = acos(cos_phi)*sign_phi;
	sin_phi = sqrt(1.-sqr_cos_phi)*sign_phi;
	cos_n_phi = cos(n*phi);
	if (fabs(sin_phi) > 1.e-4)
	  sin_n_phi_ratio = sin(n*phi)/sin_phi;
	else
	  sin_n_phi_ratio = n;
      }
      if (fabs(sin_phase) < 1.e-8)
	e += param[3]*(1.+cos_n_phi*cos_phase);
      else {
	if (sin_phi == 2.)
	  sin_phi = sqrt(1.-sqr_cos_phi)*sign_phi;
	e += param[3]*(1.+cos_n_phi*cos_phase
		       + sin_n_phi_ratio*sin_phi*sin_phase);
      }
    }
    if (energy->gradients != NULL || energy->force_constants != NULL) {
      double deriv;
      vector3 di, dj, dk, dl, ds;
      if (n == 0)
	deriv = 2.*param[3]*dphi;
      else {
	if (sin_phi == 2.)
	  sin_phi = sqrt(1.-sqr_cos_phi)*sign_phi;
	deriv = n*param[3]*(-sin_n_phi_ratio*sin_phi*cos_phase
			    +cos_n_phi*sin_phase);
      }
      vector_copy(di, rij_cross_rkj);
      vector_scale(di, lrkj/sqr(lm));
      vector_copy(dl, rkj_cross_rkl);
      vector_scale(dl, -lrkj/sqr(ln));
      ds[0] = (dot_rij_rkj*di[0]+dot_rlk_rkj*dl[0])/lrkj;
      ds[1] = (dot_rij_rkj*di[1]+dot_rlk_rkj*dl[1])/lrkj;
      ds[2] = (dot_rij_rkj*di[2]+dot_rlk_rkj*dl[2])/lrkj;
      dj[0] = ds[0]-di[0];
      dj[1] = ds[1]-di[1];
      dj[2] = ds[2]-di[2];
      dk[0] = -ds[0]-dl[0];
      dk[1] = -ds[1]-dl[1];
      dk[2] = -ds[2]-dl[2];
      if (energy->gradients != NULL) {
#ifdef GRADIENTFN
	if (energy->gradient_fn != NULL) {
	  vector3 grad;
	  vector_copy(grad, di);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, i, grad);
	  vector_copy(grad, dj);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, j, grad);
	  vector_copy(grad, dk);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, k, grad);
	  vector_copy(grad, dl);
	  vector_scale(grad, deriv);
	  (*energy->gradient_fn)(energy, l, grad);
	}
	else
#endif
        {
	  vector3 *f = (vector3 *)((PyArrayObject *)energy->gradients)->data;
	  f[i][0] += deriv*di[0];
	  f[i][1] += deriv*di[1];
	  f[i][2] += deriv*di[2];
	  f[j][0] += deriv*dj[0];
	  f[j][1] += deriv*dj[1];
	  f[j][2] += deriv*dj[2];
	  f[k][0] += deriv*dk[0];
	  f[k][1] += deriv*dk[1];
	  f[k][2] += deriv*dk[2];
	  f[l][0] += deriv*dl[0];
	  f[l][1] += deriv*dl[1];
	  f[l][2] += deriv*dl[2];
	}
      }
      if (energy->force_constants != NULL) {
	double deriv2;
	int swapij = (i > j);
	int swapik = (i > k);
	int swapil = (i > l);
	int swapjk = (j > k);
	int swapjl = (j > l);
	int swapkl = (k > l);
	vector3 ga, fa, gb, hb;
	tensor3 aga, afa, bgb, bhb, gg, fg, hg, temp;
	double ff, gg1, gg2, gg3, gg4, hh;
	int i1, i2;

	if (n == 0)
	  deriv2 = 2.*param[3];
	else
	  deriv2 = -n*n*param[3]*(cos_n_phi*cos_phase 
				  + sin_n_phi_ratio*sin_phi*sin_phase);

	/************/
	cross(ga, rkj, rij_cross_rkj);
	cross(fa, rij_cross_rkj, rij);
	cross(gb, rkj, rkj_cross_rkl);
	cross(hb, rkj_cross_rkl, rlk);
	symmetric_tensor_product(aga, rij_cross_rkj, ga, -lrkj);
	symmetric_tensor_product(afa, rij_cross_rkj, fa, -1.);
	symmetric_tensor_product(bgb, rkj_cross_rkl, gb, -lrkj);
	symmetric_tensor_product(bhb, rkj_cross_rkl, hb, -1.);
	/************/
	ff = lrkj/sqr(sqr(lm));
	gg1 = 0.5/(cube(lrkj)*sqr(lm));
	gg2 = -dot_rij_rkj/sqr(sqr(lm));
	gg3 = -0.5/(cube(lrkj)*sqr(ln));
	gg4 = dot_rlk_rkj/sqr(sqr(ln));
	hh = -lrkj/sqr(sqr(ln));
	tensor_copy(gg, aga);
	tensor_scale(gg, gg1);
	tensor_add(gg, afa, gg2);
	tensor_add(gg, bgb, gg3);
	tensor_add(gg, bhb, gg4);
	/************/
	tensor_product(fg, fa, rij_cross_rkj, lrkj);
	tensor_product(temp, rij_cross_rkj, ga, -dot_rij_rkj*lrkj);
	tensor_add(fg, temp, 1.);
	tensor_scale(fg, 1./sqr(sqr(lm)));
	tensor_product(hg, hb, rkj_cross_rkl, lrkj);
	tensor_product(temp, rkj_cross_rkl, gb, -dot_rlk_rkj*lrkj);
	tensor_add(hg, temp, 1.);
	tensor_scale(hg, -1./sqr(sqr(ln)));
	/************/

	if (energy->fc_fn != NULL) {
	  tensor3 fcii, fcjj, fckk, fcll, fcij, fcik, fcil, fcjk, fcjl, fckl;
	  double min_r = lrij;
	  if (lrkj < min_r) min_r = lrkj;
	  if (lrlk < min_r) min_r = lrlk;
	  min_r = sqr(min_r);
	  for (i1 = 0; i1 < 3; i1++)
	    for (i2 = 0; i2 < 3; i2++) {
	      fcii[i1][i2] = deriv2*di[i1]*di[i2];
	      fcjj[i1][i2] = deriv2*dj[i1]*dj[i2];
	      fckk[i1][i2] = deriv2*dk[i1]*dk[i2];
	      fcll[i1][i2] = deriv2*dl[i1]*dl[i2];
	      fcij[i1][i2] = deriv2*di[i1]*dj[i2];
	      fcik[i1][i2] = deriv2*di[i1]*dk[i2];
	      fcil[i1][i2] = deriv2*di[i1]*dl[i2];
	      fcjk[i1][i2] = deriv2*dj[i1]*dk[i2];
	      fcjl[i1][i2] = deriv2*dj[i1]*dl[i2];
	      fckl[i1][i2] = deriv2*dk[i1]*dl[i2];
	    }
	  /************/
	  tensor_add(fcii, aga, ff*deriv);
	  tensor_add(fcjj, aga, ff*deriv);
	  tensor_add(fcjj, gg, deriv);
	  tensor_add(fcjj, fg, -deriv);
	  tensor_add(fckk, bgb, hh*deriv);
	  tensor_add(fckk, gg, deriv);
	  tensor_add(fckk, hg, deriv);
	  tensor_add(fcll, bgb, hh*deriv);
	  tensor_add(fcij, aga, -ff*deriv);
	  tensor_add(fcij, fg, deriv);
	  tensor_add(fcik, fg, -deriv);
	  tensor_add(fcjk, gg, -deriv);
	  tensor_add(fcjk, fg, deriv);
	  tensor_add(fckl, bgb, -hh*deriv);
	  /************/
	  tensor_transpose(fg);
	  tensor_transpose(hg);
	  /************/
	  tensor_add(fcjj, fg, -deriv);
	  tensor_add(fckk, hg, deriv);
	  tensor_add(fcjk, hg, -deriv);
	  tensor_add(fcjl, hg, deriv);
	  tensor_add(fckl, hg, -deriv);
	  /************/
	  (*energy->fc_fn)(energy, i, i, fcii, min_r);	  
	  (*energy->fc_fn)(energy, j, j, fcjj, min_r);
	  (*energy->fc_fn)(energy, k, k, fckk, min_r);
	  (*energy->fc_fn)(energy, l, l, fcll, min_r);
	  if (swapij) {
	    tensor_transpose(fcij);
	    (*energy->fc_fn)(energy, j, i, fcij, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, i, j, fcij, min_r);
	  if (swapik) {
	    tensor_transpose(fcik);
	    (*energy->fc_fn)(energy, k, i, fcik, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, i, k, fcik, min_r);
	  if (swapil) {
	    tensor_transpose(fcil);
	    (*energy->fc_fn)(energy, l, i, fcil, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, i, l, fcil, min_r);
	  if (swapjk) {
	    tensor_transpose(fcjk);
	    (*energy->fc_fn)(energy, k, j, fcjk, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, j, k, fcjk, min_r);
	  if (swapjl) {
	    tensor_transpose(fcjl);
	    (*energy->fc_fn)(energy, l, j, fcjl, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, j, l, fcjl, min_r);
	  if (swapkl) {
	    tensor_transpose(fckl);
	    (*energy->fc_fn)(energy, l, k, fckl, min_r);
	  }
	  else
	    (*energy->fc_fn)(energy, k, l, fckl, min_r);
	}
	else {
	  double *fc_data =
	    (double *)((PyArrayObject *)energy->force_constants)->data;
	  double *fcii = fc_data + 9*input->natoms*i+3*i;
	  double *fcjj = fc_data + 9*input->natoms*j+3*j;
	  double *fckk = fc_data + 9*input->natoms*k+3*k;
	  double *fcll = fc_data + 9*input->natoms*l+3*l;
	  double *fcij, *fcik, *fcil, *fcjk, *fcjl, *fckl;
	  if (swapij)
	    fcij = fc_data + 9*input->natoms*j+3*i;
	  else
	    fcij = fc_data + 9*input->natoms*i+3*j;
	  if (swapik)
	    fcik = fc_data + 9*input->natoms*k+3*i;
	  else
	    fcik = fc_data + 9*input->natoms*i+3*k;
	  if (swapil)
	    fcil = fc_data + 9*input->natoms*l+3*i;
	  else
	    fcil = fc_data + 9*input->natoms*i+3*l;
	  if (swapjk)
	    fcjk = fc_data + 9*input->natoms*k+3*j;	
	  else
	    fcjk = fc_data + 9*input->natoms*j+3*k;
	  if (swapjl)
	    fcjl = fc_data + 9*input->natoms*l+3*j;	
	  else
	    fcjl = fc_data + 9*input->natoms*j+3*l;
	  if (swapkl)
	    fckl = fc_data + 9*input->natoms*l+3*k;	
	  else
	    fckl = fc_data + 9*input->natoms*k+3*l;
	  for (i1 = 0; i1 < 3; i1++)
	    for (i2 = 0; i2 < 3; i2++) {
	      int o = 3*input->natoms*i1 + i2;
	      fcii[o] += deriv2*di[i1]*di[i2];
	      fcjj[o] += deriv2*dj[i1]*dj[i2];
	      fckk[o] += deriv2*dk[i1]*dk[i2];
	      fcll[o] += deriv2*dl[i1]*dl[i2];
	      if (swapij)
		fcij[o] += deriv2*dj[i1]*di[i2];
	      else
		fcij[o] += deriv2*di[i1]*dj[i2];
	      if (swapik)
		fcik[o] += deriv2*dk[i1]*di[i2];
	      else
		fcik[o] += deriv2*di[i1]*dk[i2];
	      if (swapil)
		fcil[o] += deriv2*dl[i1]*di[i2];
	      else
		fcil[o] += deriv2*di[i1]*dl[i2];
	      if (swapjk)
		fcjk[o] += deriv2*dk[i1]*dj[i2];
	      else
		fcjk[o] += deriv2*dj[i1]*dk[i2];
	      if (swapjl)
		fcjl[o] += deriv2*dl[i1]*dj[i2];
	      else
		fcjl[o] += deriv2*dj[i1]*dl[i2];
	      if (swapkl)
		fckl[o] += deriv2*dl[i1]*dk[i2];
	      else
		fckl[o] += deriv2*dk[i1]*dl[i2];
	    }
	  add_fc_tensor(fcii, input->natoms, 0, aga, ff*deriv);
	  add_fc_tensor(fcjj, input->natoms, 0, aga, ff*deriv);
	  add_fc_tensor(fcij, input->natoms, swapij, aga, -ff*deriv);
	  add_fc_tensor(fcjj, input->natoms, 0, gg, deriv);
	  add_fc_tensor(fckk, input->natoms, 0, gg, deriv);
	  add_fc_tensor(fcjk, input->natoms, swapjk, gg, -deriv);
	  add_fc_tensor(fckk, input->natoms, 0, bgb, hh*deriv);
	  add_fc_tensor(fcll, input->natoms, 0, bgb, hh*deriv);
	  add_fc_tensor(fckl, input->natoms, swapkl, bgb, -hh*deriv);
	  add_fc_tensor(fcij, input->natoms, swapij, fg, deriv);
	  add_fc_tensor(fcjj, input->natoms, 0, fg, -deriv);
	  add_fc_tensor(fcjj, input->natoms, 1, fg, -deriv);
	  add_fc_tensor(fcik, input->natoms, swapik, fg, -deriv);
	  add_fc_tensor(fcjk, input->natoms, swapjk, fg, deriv);
	  add_fc_tensor(fcjk, input->natoms, !swapjk, hg, -deriv);
	  add_fc_tensor(fckk, input->natoms, 1, hg, deriv);
	  add_fc_tensor(fckk, input->natoms, 0, hg, deriv);
	  add_fc_tensor(fcjl, input->natoms, !swapjl, hg, deriv);
	  add_fc_tensor(fckl, input->natoms, !swapkl, hg, -deriv);
	}
      }
    }
    index += 4;
    param += 4;
  }
  energy->energy_terms[self->index] = e;
}
