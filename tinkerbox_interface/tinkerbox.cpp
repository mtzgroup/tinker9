#include <stdio.h>

#include "tinker9_c_function_wrapper.h"
#include "tinkerbox.h"

void initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile)
{
    internal_initialize_tinker(qm_indices, n_qm, xyzfile);
}

int32_t get_n_qm()
{
    return internal_get_n_qm();
}

int32_t get_n_mm()
{
    return internal_get_n_mm();
}

void get_qm_atomic_indices(int* qm_atomic_numbers)
{
    internal_get_qm_atomic_indices(qm_atomic_numbers);
}

void get_qm_mass(double* qm_masses)
{
    internal_get_qm_mass(qm_masses);
}

void get_mm_mass(double* mm_masses)
{
    internal_get_mm_mass(mm_masses);
}

void get_qm_xyz(double* qm_coords)
{
    internal_get_qm_xyz(qm_coords);
}

void set_qm_xyz(const double* const qm_coords)
{
    internal_set_qm_xyz(qm_coords);
}

void get_mm_xyz(double* mm_coords)
{
    internal_get_mm_xyz(mm_coords);
}

void set_mm_xyz(const double* const mm_coords)
{
    internal_set_mm_xyz(mm_coords);
}

void get_mm_charge(double* charges)
{
    internal_get_mm_charge(charges);
}



double get_energy_nonpolar_mm_contribution()
{
    return internal_get_energy_nonpolar_mm_contribution();
}

void get_gradients_all_atoms_mm_contribution(double* grad)
{
    internal_get_gradients_all_atoms_mm_contribution(grad);
}
