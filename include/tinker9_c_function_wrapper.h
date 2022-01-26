#ifndef TINKER9_C_FUNCTION_WRAPPER_H_
#define TINKER9_C_FUNCTION_WRAPPER_H_

#include <stdint.h>

void internal_initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile);
int32_t internal_get_n_qm();
int32_t internal_get_n_mm();
void internal_get_qm_atomic_indices(int* qm_atomic_numbers);
void internal_get_qm_mass(double* qm_masses);
void internal_get_mm_mass(double* mm_masses);
void internal_get_qm_xyz(double* qm_coords);
void internal_set_qm_xyz(const double* const qm_coords);
void internal_get_mm_xyz(double* mm_coords);
void internal_set_mm_xyz(const double* const mm_coords);
void internal_get_mm_charge(double* charges);



double internal_get_energy_nonpolar_mm_contribution();
void internal_get_gradients_all_atoms_mm_contribution(double* grad);

#endif // #ifndef TINKER9_C_FUNCTION_WRAPPER_
