for (k = 0; k < nblist->npairs; k++, pair++) {
  if (cutoff_sq == 0. || pair->r_sq <= cutoff_sq) {
    int i = pair->i1;
    int j = pair->i2;
    double r_sq = pair->r_sq;
    double r = sqrt(r_sq);
    double f1 = erfc(beta*r);
    double f2 = erfc_cutoff;
    vector3 rij;
    if (pair->one_four)
      f1 += one_four_factor-1.;
#ifdef NBLIST_WITH_VECTORS
    vector_copy(rij, pair->d);
#else
    distance_vector_function(rij, x[j], x[i], universe_data);
#endif
#   include "ewald2.i"
  }
}

for (k = 0; k < n_ex; k += 2) {
  int i = excluded[k];
  int j = excluded[k+1];
  double r_sq, r, f1, f2;
  vector3 rij;
  distance_vector_function(rij, x[j], x[i], universe_data);
  r_sq = vector_length_sq(rij);
  r = sqrt(r_sq);
  f1 = erfc(beta*r)-1.;
  f2 = erfc_cutoff;
# include "ewald2.i"
}
