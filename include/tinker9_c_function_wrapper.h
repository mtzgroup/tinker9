#ifndef TINKER9_C_FUNCTION_WRAPPER_H_
#define TINKER9_C_FUNCTION_WRAPPER_H_

#include <stdint.h>

void internal_initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile);
int32_t internal_get_n_mm();
double internal_get_energy_nonpolar_mm_contribution();
void internal_get_gradients_all_atoms_mm_contribution(double* grad);

#endif // #ifndef TINKER9_C_FUNCTION_WRAPPER_
