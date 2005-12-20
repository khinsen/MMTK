{
  double qiqj = charge[i]*charge[j]*electrostatic_energy_factor;
  double ef = 2./sqrt(M_PI);
  self->energy += qiqj*(f1/r-f2*inv_cutoff);
  if (gradients != NULL || force_constants != NULL) {
    double deriv = -qiqj*(f1/r_sq-f2*inv_cutoff
	                  +ef*beta*(exp(-beta*beta*r_sq)/r
                                    -inv_cutoff*exp(-beta*beta*cutoff_sq)))/r;
    if (gradients != NULL) {
      vector3 grad;
      grad[0] = deriv*rij[0];
      grad[1] = deriv*rij[1];
      grad[2] = deriv*rij[2];
      if (gradient_fn != NULL) {
	(*gradient_fn)(self, gradients, i, grad);
	vector_changesign(grad);
	(*gradient_fn)(self, gradients, j, grad);
      }
      else {
	vector3 *f = (vector3 *)((PyArrayObject *)gradients)->data;
	f[i][0] += grad[0];
	f[i][1] += grad[1];
	f[i][2] += grad[2];
	f[j][0] -= grad[0];
	f[j][1] -= grad[1];
	f[j][2] -= grad[2];
      }
    }
    if (force_constants != NULL) {
      double deriv2 = 2.*qiqj*(f1/(r*r_sq) - f2*cube(inv_cutoff)
	                       + ef*beta*exp(-beta*beta*r_sq)
                                               * (beta*beta+1./r_sq));
      if (inv_cutoff > 0.)
	deriv2 -= 2.*qiqj*ef*beta*exp(-beta*beta*cutoff_sq)
                         * (beta*beta+sqr(inv_cutoff));
      add_pair_fc(self, force_constants, fc_fn, i, j,
		  rij, r_sq, deriv, deriv2);
    }
  }
}
