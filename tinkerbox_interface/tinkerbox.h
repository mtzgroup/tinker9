#ifndef TINKERBOX_H_
#define TINKERBOX_H_

#define DSOGLOBAL __attribute__((visibility("default")))

#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

  /**
   * \brief Read the tinker input files (key, prm and xyz) and setup most global variables in tinker.
   * Should be called before you call any other functions in this file.
   *
   * @param qm_indices Should have a size of n_qm.
   *                   Important: The indices must be ONE-indexing!
   * @param n_qm Self-explanatory.
   * @param xyzfile Filename for xyz file. The filename should have less than 240 chars.
   **/
  DSOGLOBAL void initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile);
  
  /**
   * \brief Get number of QM atoms.
   **/
  DSOGLOBAL int32_t get_n_qm();
  
  /**
   * \brief Get number of MM atoms.
   **/
  DSOGLOBAL int32_t get_n_mm();

  /**
   * \brief Get the QM atomic numbers from tinker.
   *
   * @param qm_atomic_numbers Should have n_qm space allocated beforehand.
   **/
  DSOGLOBAL void get_qm_atomic_indices(int* qm_atomic_numbers);

  /**
   * \brief Get the QM atom mass from tinker. In AMU, not in AU.
   *
   * @param qm_masses Should have n_qm space allocated beforehand.
   **/
  DSOGLOBAL void get_qm_mass(double* qm_masses);

  /**
   * \brief Get the MM atom mass from tinker. In AMU, not in AU.
   *
   * @param mm_masses Should have n_mm space allocated beforehand.
   **/
  DSOGLOBAL void get_mm_mass(double* mm_masses);

  /**
   * \brief Get the QM atom coordinates from tinker, in Bohr.
   *
   * @param qm_coords Should have (3 * n_qm) space allocated beforehand.
   *                  Returned coordinates are ordered in x0,y0,z0,x1,y1,z1,...
   **/
  DSOGLOBAL void get_qm_xyz(double* qm_coords);

  /**
   * \brief Set the QM atom coordinates in tinker. qm_coords should be in Bohr.
   *
   * @param qm_coords Should have (3 * n_qm) in size.
   **/
  DSOGLOBAL void set_qm_xyz(const double* const qm_coords);

  /**
   * \brief Get the MM atom coordinates from tinker, in Bohr.
   *
   * @param mm_coords Should have (3 * n_mm) space allocated beforehand.
   *                  Returned coordinates are ordered in x0,y0,z0,x1,y1,z1,...
   **/
  DSOGLOBAL void get_mm_xyz(double* mm_coords);

  /**
   * \brief Set the MM atom coordinates in tinker. mm_coords should be in Bohr.
   *
   * @param mm_coords Should have (3 * n_mm) in size.
   **/
  DSOGLOBAL void set_mm_xyz(const double* const mm_coords);

  /**
   * \brief Get the MM point charges, specified in prm file "charge" or "multipole" keyword, in electron.
   * Important: The "charge" keyword and "multipole" keyword should not be used together in tinker input!
   * Important: In order for tinker to operate properly, for each atom_type specified in prm file,
   *            either none of them should have a "multipole" parameter, or all of them must have "multipole" parameters.
   *            If you run tinker alone, and a "multipole" is missing for one atom_type, tinker will fill it with 0.
   *            This is NOT working for terachem-tinker interface.
   *
   * @param charges Should have n_mm space allocated beforehand.
   **/
  DSOGLOBAL void get_mm_charge(double* charges);

  /**
   * \brief Get the MM static dipoles, specified in prm file "multipole" keyword, in global coordinate, in electron*Bohr.
   * Important: The "multipole" parameters have Bohr in length unit! See "kmpole.f", there's a unit conversion line.
   * Important: If Z-only mode is used for local coordinate definition, make sure there's no local x and y component
   *            of that dipole, otherwise the dipole obtained from this function is inconsistent with
   *            that used for energy calculation.
   *
   * @param dipoles Should have (3 * n_mm) space allocated beforehand.
   *                Returned dipoles are ordered in x0,y0,z0,x1,y1,z1,...
   * @return "npole" variable in tinker memory, suppose to return (n_qm + n_mm). For error checking.
   **/
  DSOGLOBAL int32_t get_mm_static_point_dipole(double* dipoles);

  /**
   * \brief Get the MM polarizability, specified in prm file "polarize" keyword, in Bohr^3.
   *        Only isotropic polarizability is supported now.
   * 
   * @param dipoles Should have n_mm space allocated beforehand.
   **/
  DSOGLOBAL void get_mm_polarizibility(double* polarizabilities);

  /**
   * \brief Get the MM contribution for the total energy, in Hartree.
   *        Polarization energy is NOT included!
   **/
  DSOGLOBAL double get_energy_nonpolar_mm_contribution();

  /**
   * \brief Get the MM contribution for the QM and MM gradient, in Hartree/Bohr.
   *
   * @param grad Should have (3 * (n_qm + n_mm)) space allocated beforehand.
   **/
  DSOGLOBAL void get_gradients_all_atoms_mm_contribution(double* grad);

  /**
   * \brief Append the implicit part of electric field derivatives to the MM gradient.
   *        See the text between equation 15 and 16 of the following paper for detail:
   *        Lipparini, F.; Lagardère, L.; Stamm, B.; Cancès, E.; Schnieders, M.; Ren, P.; Maday, Y.; Piquemal, J.-P.,
   *          Scalable Evaluation of Polarization Energy and Associated Forces in Polarizable Molecular Dynamics: I. Toward Massively Parallel Direct Space Computations.
   *          Journal of Chemical Theory and Computation 2014, 10 (4), 1638-1651.
   * 
   * @param mm_torque tau_i = mu_i \cross E^QM_i at each static dipole position
   *                  Should have the size of (3 * n_mm)
   * @param mm_grad Should have (3 * n_mm) space allocated beforehand.
   **/
  DSOGLOBAL void append_gradient_from_static_dipole_rotation(const double* const mm_torque, double* mm_grad);
  
  /**
   * \brief Get the MM contribution of electric field at each MM atom position.
   *        The electric field will have a unit of electron/Bohr^2.
   * Important: Although not sure if necessary, please make sure a Tinker energy calculation
   *            (get_energy()) was carried out before calling this function.
   * 
   * @param electric_field_direct_mm Should have (n_mm * 3) space allocated beforehand.
   * @param electric_field_polarization_mm Should have (n_mm * 3) space allocated beforehand.
   **/
  DSOGLOBAL void get_electric_field_mm_contribution(double* electric_field_direct_mm, double* electric_field_polarization_mm);

  /**
   * \brief Get the induced dipole based on the electric field at each MM position.
   *        The electric field passed in should have both MM and QM contribution, in electron/Bohr^2.
   *        The obtained induced dipole will have a unit of electron*Bohr.
   * Important: Although not sure if necessary, please make sure a Tinker energy calculation
   *            (get_energy()) was carried out before calling this function.
   * 
   * @param electric_field_direct_mm Should be (n_mm * 3) in size.
   * @param electric_field_polarization_mm Should be (n_mm * 3) in size.
   * @param induced_dipole_direct Should have (n_mm * 3) space allocated beforehand.
   * @param induced_dipole_polarization Should have (n_mm * 3) space allocated beforehand.
   **/
  DSOGLOBAL void evaluate_induced_dipole_from_total_electric_field(
    const double* const electric_field_direct_mm,
    const double* const electric_field_polarization_mm,
    double* induced_dipole_direct,
    double* induced_dipole_polarization
  );

  /**
   * \brief Get the induced dipole stored in tinker memory. Only call this function after polarizable SCF converges.
   *        The obtained induced dipole will have a unit of electron*Bohr.
   * 
   * @param mu_d Should have (n_mm * 3) space allocated beforehand.
   * @param mu_p Should have (n_mm * 3) space allocated beforehand.
   **/
  DSOGLOBAL void get_mm_induced_dipole(double* mu_d, double* mu_p);

#ifdef __cplusplus
}
#endif

#endif // #ifndef TINKERBOX_
