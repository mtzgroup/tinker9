#include <stdio.h>

#include "tinker9_c_function_wrapper.h"
#include "tinkerbox.h"

void initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile)
{
    internal_initialize_tinker(qm_indices, n_qm, xyzfile);
}

int32_t get_n_mm()
{
    return internal_get_n_mm();
}

double get_energy_nonpolar_mm_contribution()
{
    return internal_get_energy_nonpolar_mm_contribution();
}

void get_gradients_all_atoms_mm_contribution(double* grad)
{
    internal_get_gradients_all_atoms_mm_contribution(grad);
}
