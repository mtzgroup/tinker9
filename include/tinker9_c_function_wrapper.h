#ifndef TINKER9_C_FUNCTION_WRAPPER_H_
#define TINKER9_C_FUNCTION_WRAPPER_H_

#include <stdint.h>

void internal_initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile);
int32_t internal_get_n_qm();
int32_t internal_get_n_mm();
void internal_get_qm_atomic_indices(int32_t* qm_atomic_numbers);
void internal_get_qm_mass(double* qm_masses);
void internal_get_mm_mass(double* mm_masses);
void internal_get_qm_xyz(double* qm_coords);
void internal_set_qm_xyz(const double* const qm_coords);
void internal_get_mm_xyz(double* mm_coords);
void internal_set_mm_xyz(const double* const mm_coords);
void internal_get_mm_charge(double* charges);
void internal_get_mm_static_point_dipole_and_quadrupole(double* dipoles, double* quadrupoles);
void internal_get_mm_polarizability(double* polarizabilities);
double internal_get_energy_nonpolar_mm_contribution();
void internal_get_gradients_all_atoms_mm_contribution(double* grad);
void internal_append_gradient_from_static_dipole_rotation(const double* const mm_torque, double* mm_grad);
void internal_get_electric_field_mm_contribution(double* electric_field_direct_mm, double* electric_field_polarization_mm);
void internal_evaluate_induced_dipole_from_total_electric_field(
        const double* const electric_field_direct_mm,
        const double* const electric_field_polarization_mm,
        double* induced_dipole_direct,
        double* induced_dipole_polarization
    );
void internal_get_mm_induced_dipole(double* mu_d, double* mu_p);

#endif // #ifndef TINKER9_C_FUNCTION_WRAPPER_
