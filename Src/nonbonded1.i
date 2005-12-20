for (k = 0; k < n_ex; k += 2) {
  int i = excluded[k];
  int j = excluded[k+1];
  double qiqj = charge[i]*charge[j]*electrostatic_energy_factor;
  pair_es_energy(-qiqj, distance_vector_function);
}
for (k = 0; k < n_14; k += 2) {
  int i = one_four[k];
  int j = one_four[k+1];
  double qiqj = charge[i]*charge[j]*electrostatic_energy_factor;
  pair_es_energy(qiqj*(es_one_four_factor-1.), distance_vector_function);
}
