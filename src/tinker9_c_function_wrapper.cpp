#include <stdio.h>

#include "tinker9_c_function_wrapper.h"
#include "tinker8_fortran_function_wrapper.h"
#include "qmmm_global.h"

#include "tinker_rt.h" // initial(), routine.h::tinker_f_*()
#include "md.h" // mdcalc.h::calc::*, mdegv.h::copy_gradient()
#include "energy.h" // energy()
#include "potent.h" // use_potent(), *_term
#include "nblist.h"
#include "elec.h" // chkpole(), rotpole()
#include "field.h" // dfield()
#include "tool/darray.h"

#include <tinker/detail/atoms.hh>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/urey.hh>
#include <tinker/detail/tors.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/units.hh>

#define DIE(...) { printf(__VA_ARGS__); \
    printf("DIE called at line number %d in file %s\n", __LINE__, __FILE__); \
    exit(1); }

template <class T>
static bool if_in_list(const T* const list, const size_t list_length, const T& item)
{
    for (size_t i = 0; i < list_length; i++)
        if (item == list[i])
            return true;
    return false;
}

void internal_initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile)
{
    /**
     * Implemented in src/initial.cpp
     * This function initialize a bunch of global variables to zero-start.
     * It also prints tinker header message (TINKER9_PROMO_STRING).
     */
    tinker::initial();

    /**
     * Implemented in tinker/tinker8_fortran_function_wrapper.f
     * This function reads in the xyz file and key file.
     */
    const int32_t xyzfile_string_length = strlen(xyzfile);
    tinkerbox_getxyz_(xyzfile, &xyzfile_string_length);

    /**
     * Implemented in tinker/mechanic.f
     * This is the function where force field parameters are read in.
     */
    tinker_f_mechanic();

    /**
     * Implemented in src/mechanic.cpp, calling src/osrw.cpp::osrw_mech()
     * Initialize Orthogonal Space Random Walk trash.
     */
    tinker::mechanic2();

    if (n_qm > 0)
        printf("TC anchor: Removing bonding terms of QM atoms.\n");
    // Note: for the following code, the index from fortran side are all one-indexed,
    //       and we expect the index in qm_indices is also one-indexed, so no +-1 is needed.
    /**
     * Zero out bonding interaction between QM atoms.
     * At QM-MM interface, the QM-MM bonds are NOT removed.
     */
    for (size_t i_bond = 0; i_bond < tinker::bndstr::nbond; i_bond++)
    {
        int32_t i_atom_1 = tinker::bndstr::ibnd[i_bond * 2 + 0],
                i_atom_2 = tinker::bndstr::ibnd[i_bond * 2 + 1];
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_1) && if_in_list<int32_t>(qm_indices, n_qm, i_atom_2))
            tinker::bndstr::bk[i_bond] = 0.0;
    }

    /**
     * Zero out bond angle interaction of type QM-QM-QM or QM-QM-MM.
     */
    for (size_t i_angle = 0; i_angle < tinker::angbnd::nangle; i_angle++)
    {
        int32_t i_atom_1 = tinker::angbnd::iang[i_angle * 3 + 0],
                i_atom_2 = tinker::angbnd::iang[i_angle * 3 + 1],
                i_atom_3 = tinker::angbnd::iang[i_angle * 3 + 2];
        int count_in_qm = 0;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_1)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_2)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_3)) count_in_qm++;
        if (count_in_qm > 1)
            tinker::angbnd::ak[i_angle] = 0.0;
    }

    /**
     * Zero out Urey-Bradley interaction the samw way as for bond angle.
     * This is included because Amoeba water model has Urey-Bradley parameters.
     */
    for (size_t i_urey = 0; i_urey < tinker::urey::nurey; i_urey++)
    {
        int32_t i_atom_1 = tinker::urey::iury[i_urey * 3 + 0],
                i_atom_2 = tinker::urey::iury[i_urey * 3 + 1],
                i_atom_3 = tinker::urey::iury[i_urey * 3 + 2];
        int count_in_qm = 0;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_1)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_2)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_3)) count_in_qm++;
        if (count_in_qm > 1)
            tinker::urey::uk[i_urey] = 0.0;
    }
    
    /**
     * Zero out dihedral interaction if, among the 4 atoms, 3 or more of them are QM atoms.
     */
    for (size_t i_torsion = 0; i_torsion < tinker::tors::ntors; i_torsion++)
    {
        int32_t i_atom_1 = tinker::tors::itors[i_torsion * 4 + 0],
                i_atom_2 = tinker::tors::itors[i_torsion * 4 + 1],
                i_atom_3 = tinker::tors::itors[i_torsion * 4 + 2],
                i_atom_4 = tinker::tors::itors[i_torsion * 4 + 3];
        int count_in_qm = 0;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_1)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_2)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_3)) count_in_qm++;
        if (if_in_list<int32_t>(qm_indices, n_qm, i_atom_4)) count_in_qm++;
        if (count_in_qm > 2)
        {
            // The data structure here is: tors[x] is an array of size 4 * n_torsion.
            // For each torsion angle, the 4 fields are: potential, phase offset angle, cosine of phase angle, sine of phase angle.
            // tors[x] is the parameters for the [x]-fold periodicity of the torsion angle.
            tinker::tors::tors1[i_torsion * 4 + 0] = 0.0;
            tinker::tors::tors2[i_torsion * 4 + 0] = 0.0;
            tinker::tors::tors3[i_torsion * 4 + 0] = 0.0;
            tinker::tors::tors4[i_torsion * 4 + 0] = 0.0;
            tinker::tors::tors5[i_torsion * 4 + 0] = 0.0;
            tinker::tors::tors6[i_torsion * 4 + 0] = 0.0;
        }
    }

    /**
     * VDW interactions between QM atoms are removed in src/evdw.cpp::evdw_data() and src/echglj.cpp::echglj_data()
     */

    /**
     * Charges of QM atoms are zeroed out in src/elec.cpp::pchg_data()
     * Notice: Although there's an array pchg in tinker8 source code, that is of the size n_total, but it's not used in tinker9!
     *         Tinker9 recompute pchg array in src/elec.cpp::pchg_data().
     */

    /**
     * Zero out multipoles of QM atoms.
     * Notice: In tinker8 the second dimension of pole array is maxpole=13, in tinker9 it's mpl_total=10.
     */
    if (tinker::use_potent(tinker::mpole_term))
    {
        if (n_qm > 0)
            printf("TC anchor: Removing multipoles of QM atoms.\n");

        for (size_t i_i_qm = 0; i_i_qm < n_qm; i_i_qm++)
        {
            int32_t i_qm = qm_indices[i_i_qm] - 1; // One-indexed to zero-indexed
            for (size_t i_pole_coordinate = 0; i_pole_coordinate < tinker::mpole::maxpole; i_pole_coordinate++)
                tinker::mpole::pole[i_qm * tinker::mpole::maxpole + i_pole_coordinate] = 0.0;
        }
    }

    /**
     * Zero out polarizabilities of QM atoms.
     */
    if (tinker::use_potent(tinker::polar_term))
    {
        if (n_qm > 0)
            printf("TC anchor: Removing polarizabilities of QM atoms.\n");

        for (size_t i_i_qm = 0; i_i_qm < n_qm; i_i_qm++)
        {
            int32_t i_qm = qm_indices[i_i_qm] - 1; // One-indexed to zero-indexed
            tinker::polar::polarity[i_qm] = 0.0;
        }
    }

    /**
     * Stores these variables in tinker global memory.
     * They're used in src/rc_man.cpp::initialize() -> device_data() -> src/energy.cpp::energy_data() -> src/evdw.cpp::evdw_data(),
     * so must be set before initalize() function call.
     */
    QMMMGlobal::n_qm = n_qm;
    QMMMGlobal::qm_indices = new int32_t[n_qm];
    memcpy(QMMMGlobal::qm_indices, qm_indices, n_qm * sizeof(int32_t));

    // tinker::n is defined in include/mdpq.h, set in src/rc_man.cpp::initialize() -> src/mdpq.cpp::n_data(), and it's NOT set already.
    QMMMGlobal::n_mm = tinker::atoms::n - n_qm;
    QMMMGlobal::mm_indices = new int32_t[QMMMGlobal::n_mm];
    for (int32_t i_total = 1, i_mm = 0; i_total <= tinker::atoms::n; i_total++)
        if (!if_in_list<int32_t>(qm_indices, n_qm, i_total))
        {
            QMMMGlobal::mm_indices[i_mm] = i_total;
            i_mm++;
        }

    /**
     * Copied from src/testgrad_x.cpp::x_testgrad()
     * rc_flag defined in include/mdpq.h, calc namespace defined in include/mdcalc.h
     */
    int flags = (tinker::calc::xyz + tinker::calc::mass);
    flags += (tinker::calc::energy + tinker::calc::grad);
    tinker::rc_flag = flags;
    
    /**
     * Implemented in src/rc_man.cpp
     * Initialize Tinker device resource manager.
     * Specifically, it calls src/cudart/gpu_card.cpp::gpu_card_data() to query GPU information,
     * then allocate GPU memory in subsequent rc_man constructor calls.
     */
    tinker::initialize();

    /**
     * Allocate helper arrays for QMMM induced dipole SCF.
     * TODO: not deallocated
     */
    if (tinker::use_potent(tinker::polar_term))
    {
        tinker::darray::allocate(tinker::n, &QMMMGlobal::d_qmmm_electric_field_d);
        tinker::darray::allocate(tinker::n, &QMMMGlobal::d_qmmm_electric_field_p);
    }
    
    /**
     * Print the GPU environment. It seems to use OpenACC as default setting.
     */
#if TINKER_CUDART
    printf("TC anchor: Checking GPU compilation: Tinker9 compiled with CUDA\n");
    if (tinker::vlist_version() & tinker::NBL_SPATIAL)
        printf("TC anchor: Checking GPU runtime: vdw_term running with OpenACC\n");
    else
        printf("TC anchor: Checking GPU runtime: vdw_term running with OpenACC\n");
    if (tinker::clist_version() & tinker::NBL_SPATIAL)
        printf("TC anchor: Checking GPU runtime: charge_term running with OpenACC\n");
    else
        printf("TC anchor: Checking GPU runtime: charge_term running with OpenACC\n");
    if (tinker::mlist_version() & tinker::NBL_SPATIAL)
        printf("TC anchor: Checking GPU runtime: mpole_term, polar_term, chgtrn_term, repuls_term running with OpenACC\n");
    else
        printf("TC anchor: Checking GPU runtime: mpole_term, polar_term, chgtrn_term, repuls_term running with OpenACC\n");
    if (tinker::ulist_version() & tinker::NBL_SPATIAL)
        printf("TC anchor: Checking GPU runtime: polar_term running with OpenACC\n");
    else
        printf("TC anchor: Checking GPU runtime: polar_term running with OpenACC\n");
#else
    printf("TC anchor: Checking GPU compilation: Tinker9 compiled with only OpenACC\n");
    printf("TC anchor: Checking GPU runtime: All running with OpenACC\n");
#endif
    printf("\n\n");
}

int32_t internal_get_n_qm()
{
    return QMMMGlobal::n_qm;
}

int32_t internal_get_n_mm()
{
    return QMMMGlobal::n_mm;
}

void internal_get_qm_atomic_indices(int32_t* qm_atomic_numbers)
{
    for (size_t i_i_qm = 0; i_i_qm < QMMMGlobal::n_qm; i_i_qm++)
    {
        int32_t i_qm = QMMMGlobal::qm_indices[i_i_qm] - 1; // One-index to zero-index
        qm_atomic_numbers[i_i_qm] = tinker::atomid::atomic[i_qm];
    }
}

void internal_get_qm_mass(double* qm_masses)
{
    for (size_t i_i_qm = 0; i_i_qm < QMMMGlobal::n_qm; i_i_qm++)
    {
        int32_t i_qm = QMMMGlobal::qm_indices[i_i_qm] - 1; // One-index to zero-index
        qm_masses[i_i_qm] = tinker::atomid::mass[i_qm];
    }
}

void internal_get_mm_mass(double* mm_masses)
{
    for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
    {
        int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
        mm_masses[i_i_mm] = tinker::atomid::mass[i_mm];
    }
}

void internal_get_qm_xyz(double* qm_coords)
{
    const size_t n_total = tinker::n;
    tinker::real* all_coords = new tinker::real[n_total * 3];
    // tinker::x is defined in include/mdpq.h, set in src/mdpq.cpp::xyz_data().
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 0, tinker::x);
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 1, tinker::y);
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 2, tinker::z);
    tinker::wait_for(tinker::g::q0);

    for (size_t i_i_qm = 0; i_i_qm < QMMMGlobal::n_qm; i_i_qm++)
    {
        int32_t i_qm = QMMMGlobal::qm_indices[i_i_qm] - 1; // One-index to zero-index
        for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            qm_coords[i_i_qm * 3 + i_xyz] = all_coords[i_xyz * n_total + i_qm] / tinker::units::bohr;
    }

    delete[] all_coords;
}

void internal_set_qm_xyz(const double* const qm_coords)
{
    if (tinker::calc::traj & tinker::rc_flag)
        DIE("calc::traj should not be used in QMMM interface.\n")

    for (size_t i_i_qm = 0; i_i_qm < QMMMGlobal::n_qm; i_i_qm++)
    {
        int32_t i_qm = QMMMGlobal::qm_indices[i_i_qm] - 1; // One-index to zero-index

        tinker::atoms::x[i_qm] = qm_coords[i_i_qm * 3 + 0] * tinker::units::bohr;
        tinker::atoms::y[i_qm] = qm_coords[i_i_qm * 3 + 1] * tinker::units::bohr;
        tinker::atoms::z[i_qm] = qm_coords[i_i_qm * 3 + 2] * tinker::units::bohr;
    }

    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::xpos, tinker::atoms::x);
    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::ypos, tinker::atoms::y);
    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::zpos, tinker::atoms::z);
    tinker::copy_pos_to_xyz();
}

void internal_get_mm_xyz(double* mm_coords)
{
    const size_t n_total = tinker::n;
    tinker::real* all_coords = new tinker::real[n_total * 3];
    // tinker::x is defined in include/mdpq.h, set in src/mdpq.cpp::xyz_data().
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 0, tinker::x);
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 1, tinker::y);
    tinker::darray::copyout(tinker::g::q0, n_total, all_coords + n_total * 2, tinker::z);
    tinker::wait_for(tinker::g::q0);

    for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
    {
        int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
        for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            mm_coords[i_i_mm * 3 + i_xyz] = all_coords[i_xyz * n_total + i_mm] / tinker::units::bohr;
    }

    delete[] all_coords;
}

void internal_set_mm_xyz(const double* const mm_coords)
{
    if (tinker::calc::traj & tinker::rc_flag)
        DIE("calc::traj should not be used in QMMM interface.\n")

    for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
    {
        int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
        
        tinker::atoms::x[i_mm] = mm_coords[i_i_mm * 3 + 0] * tinker::units::bohr;
        tinker::atoms::y[i_mm] = mm_coords[i_i_mm * 3 + 1] * tinker::units::bohr;
        tinker::atoms::z[i_mm] = mm_coords[i_i_mm * 3 + 2] * tinker::units::bohr;
    }

    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::xpos, tinker::atoms::x);
    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::ypos, tinker::atoms::y);
    tinker::darray::copyin(tinker::g::q0, tinker::n, tinker::zpos, tinker::atoms::z);
    tinker::copy_pos_to_xyz();
}

void internal_get_mm_charge(double* charges)
{
    int n_charge_source = 0;

    if (tinker::use_potent(tinker::charge_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_charges = new tinker::real[n_total];
        // tinker::pchg is defined in include/mod.charge.h, set in src/elec.cpp::pchg_data().
        tinker::darray::copyout(tinker::g::q0, n_total, all_charges, tinker::pchg);
        tinker::wait_for(tinker::g::q0);

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            charges[i_i_mm] = all_charges[i_mm];
        }

        delete[] all_charges;
        n_charge_source++;
    }

    if (tinker::use_potent(tinker::mpole_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_multipole = new tinker::real[n_total * tinker::mpl_total];
        // tinker::pole and tinker::rpole are defined in include/mod.mpole.h.
        // tinker::pole is set in src/elec.cpp::pole_data(), and tinker::rpole is computed in src/elec.cpp::mpole_init() -> src/acc/rotpole.cpp::rotpole()
        tinker::darray::copyout(tinker::g::q0, n_total, all_multipole, tinker::pole);
        tinker::wait_for(tinker::g::q0);

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            charges[i_i_mm] = all_multipole[i_mm * tinker::mpl_total + 0];
        }

        delete[] all_multipole;
        n_charge_source++;
    }

    if (n_charge_source == 0)
    {
        printf("TC anchor: Warning: the charge is accessed, but neither charge nor multipole is specified in Tinker parameter.\n");
        memset(charges, 0, QMMMGlobal::n_mm * sizeof(double));
    }
    else if (n_charge_source > 1)
        DIE("Both charge and multipole (used automatically with polarize) is specified in Tinker parameter, which will cause inconsistency!\n")
}

int32_t internal_get_mm_static_point_dipole(double* dipoles)
{
    tinker::chkpole();
    tinker::rotpole();

    if (tinker::use_potent(tinker::mpole_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_multipole = new tinker::real[n_total * tinker::mpl_total];
        tinker::darray::copyout(tinker::g::q0, n_total, all_multipole, tinker::rpole);
        tinker::wait_for(tinker::g::q0);

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
                dipoles[i_i_mm * 3 + i_xyz] = all_multipole[i_mm * tinker::mpl_total + (i_xyz + 1)] / tinker::units::bohr;
        }

        delete[] all_multipole;
        return tinker::n;
    }
    else
    {
        printf("TC anchor: Warning: the static dipole is accessed, but mpole is not specified in Tinker parameter.\n");
        memset(dipoles, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
        return 0;
    }
}

void internal_get_mm_polarizability(double* polarizabilities)
{
    if (tinker::use_potent(tinker::polar_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_polarizabilities = new tinker::real[n_total];
        // tinker::polarity is defined in include/mod.polar.h, set in src/epolar.cpp::epolar_data().
        tinker::darray::copyout(tinker::g::q0, n_total, all_polarizabilities, tinker::polarity);
        tinker::wait_for(tinker::g::q0);

        const double angstrom_to_bohr_cube = 1.0 / tinker::units::bohr / tinker::units::bohr / tinker::units::bohr;

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            polarizabilities[i_i_mm] = all_polarizabilities[i_mm] * angstrom_to_bohr_cube;
        }

        delete[] all_polarizabilities;
    }
    else
    {
        printf("TC anchor: Warning: the polarizability is accessed, but polar is not specified in Tinker parameter.\n");
        memset(polarizabilities, 0, QMMMGlobal::n_mm * sizeof(double));
    }
}

double internal_get_energy_nonpolar_mm_contribution()
{
    tinker::energy(tinker::rc_flag);

    printf("TC anchor: removing energy_ep from total energy in get_energy_nonpolar_mm_contribution(), energy_ep = %.10f Kcal/mol\n", tinker::energy_ep);

    return (tinker::esum - tinker::energy_ep) / tinker::units::hartree;
}

void internal_get_gradients_all_atoms_mm_contribution(double* grad)
{
    tinker::energy(tinker::rc_flag);

    std::vector<double> gdx(tinker::n), gdy(tinker::n), gdz(tinker::n);
    tinker::copy_gradient(tinker::calc::grad, gdx.data(), gdy.data(), gdz.data());

    for (size_t i_atom = 0; i_atom < tinker::n; i_atom++)
    {
        grad[i_atom * 3 + 0] = gdx[i_atom];
        grad[i_atom * 3 + 1] = gdy[i_atom];
        grad[i_atom * 3 + 2] = gdz[i_atom];
    }
}

void internal_append_gradient_from_static_dipole_rotation(const double* const mm_torque, double* mm_grad)
{
    if (tinker::use_potent(tinker::mpole_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_torque = new tinker::real[n_total * 3];
        tinker::grad_prec* all_gradient = new tinker::grad_prec[n_total * 3]; // This can be fixed point number!
        memset(all_torque, 0, n_total * 3 * sizeof(tinker::real));
        
        tinker::darray::zero(tinker::g::q0, n_total, tinker::depx, tinker::depy, tinker::depz);

        // The torque has the same unit as energy. Here we multiply a Coulomb constant, convert e^2/Bohr to KCal/mol.
        const double torqueAU_to_kcalPerMol = 1.0 / tinker::units::bohr * tinker::chgpot::electric;

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
                all_torque[i_xyz * n_total + i_mm] = mm_torque[i_i_mm * 3 + i_xyz] * torqueAU_to_kcalPerMol;
        }
        
        tinker::darray::copyin(tinker::g::q0, n_total, tinker::trqx, all_torque + n_total * 0);
        tinker::darray::copyin(tinker::g::q0, n_total, tinker::trqy, all_torque + n_total * 1);
        tinker::darray::copyin(tinker::g::q0, n_total, tinker::trqz, all_torque + n_total * 2);
        tinker::wait_for(tinker::g::q0);
        
        const int vers = tinker::calc::grad;
        tinker::torque(vers, tinker::depx, tinker::depy, tinker::depz);

        tinker::darray::copyout(tinker::g::q0, n_total, all_gradient + n_total * 0, tinker::depx);
        tinker::darray::copyout(tinker::g::q0, n_total, all_gradient + n_total * 1, tinker::depy);
        tinker::darray::copyout(tinker::g::q0, n_total, all_gradient + n_total * 2, tinker::depz);
        tinker::wait_for(tinker::g::q0);

        const double kcalPerMolPerAngstrom_to_hartreePerBohr = tinker::units::bohr / tinker::units::hartree;

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
                mm_grad[i_i_mm * 3 + i_xyz] += tinker::to_flt_host<double>(all_gradient[i_xyz * n_total + i_mm]) * kcalPerMolPerAngstrom_to_hartreePerBohr;
        }

        delete[] all_torque;
        delete[] all_gradient;
    }
    else
    {
        printf("TC anchor: Warning: the static dipole rotation is accessed, but mpole is not specified in Tinker parameter.\n");
    }
}

void internal_get_electric_field_mm_contribution(double* electric_field_direct_mm, double* electric_field_polarization_mm)
{
    if (tinker::use_potent(tinker::polar_term))
    {
        auto* Ed = QMMMGlobal::d_qmmm_electric_field_d;
        auto* Ep = QMMMGlobal::d_qmmm_electric_field_p;

        tinker::dfield(Ed, Ep);
        
        const size_t n_total = tinker::n;
        tinker::real* all_Ed = new tinker::real[n_total * 3];
        tinker::real* all_Ep = new tinker::real[n_total * 3];

        tinker::darray::copyout(tinker::g::q0, n_total, all_Ed, Ed);
        tinker::darray::copyout(tinker::g::q0, n_total, all_Ep, Ep);
        tinker::wait_for(tinker::g::q0);

        const double bohr_to_angstrom_square = tinker::units::bohr * tinker::units::bohr;

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                electric_field_direct_mm[i_i_mm * 3 + i_xyz] = all_Ed[i_mm * 3 + i_xyz] * bohr_to_angstrom_square;
                electric_field_polarization_mm[i_i_mm * 3 + i_xyz] = all_Ep[i_mm * 3 + i_xyz] * bohr_to_angstrom_square;
            }
        }

        delete[] all_Ed;
        delete[] all_Ep;
    }
    else
    {
        printf("TC anchor: Warning: the electric field from MM field is accessed, but polar is not specified in Tinker parameter.\n"
               "           In this case Tinker9 has some global variables not initialized, and it'll cause segmentation fault if you call dfield() function.\n"
               "           So, we just return zero electric fields.\n");
        memset(electric_field_direct_mm, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
        memset(electric_field_polarization_mm, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
    }
}

void internal_evaluate_induced_dipole_from_total_electric_field(
        const double* const electric_field_direct_mm,
        const double* const electric_field_polarization_mm,
        double* induced_dipole_direct,
        double* induced_dipole_polarization
    )
{
    if (tinker::use_potent(tinker::polar_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* workspace_n_mm_times_6 = new tinker::real[n_total * 3 * 2];

        tinker::real* all_Ed = workspace_n_mm_times_6;
        tinker::real* all_Ep = workspace_n_mm_times_6 + n_total * 3;
        // Zero field for QM atoms
        memset(all_Ed, 0, n_total * 3 * sizeof(tinker::real));
        memset(all_Ep, 0, n_total * 3 * sizeof(tinker::real));

        const double angstrom_to_bohr_square = 1.0 / tinker::units::bohr / tinker::units::bohr;

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                all_Ed[i_mm * 3 + i_xyz] = electric_field_direct_mm[i_i_mm * 3 + i_xyz] * angstrom_to_bohr_square;
                all_Ep[i_mm * 3 + i_xyz] = electric_field_polarization_mm[i_i_mm * 3 + i_xyz] * angstrom_to_bohr_square;
            }
        }

        tinker::darray::copyin(tinker::g::q0, n_total, QMMMGlobal::d_qmmm_electric_field_d, all_Ed);
        tinker::darray::copyin(tinker::g::q0, n_total, QMMMGlobal::d_qmmm_electric_field_p, all_Ep);
        tinker::wait_for(tinker::g::q0);

        // Actual induce() function call
        tinker::induce_mutual_pcg1(tinker::uind, tinker::uinp);

        tinker::real* all_mu_d = workspace_n_mm_times_6;
        tinker::real* all_mu_p = workspace_n_mm_times_6 + n_total * 3;

        tinker::darray::copyout(tinker::g::q0, n_total, all_mu_d, tinker::uind);
        tinker::darray::copyout(tinker::g::q0, n_total, all_mu_p, tinker::uinp);
        tinker::wait_for(tinker::g::q0);

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                induced_dipole_direct[i_i_mm * 3 + i_xyz] = all_mu_d[i_mm * 3 + i_xyz] / tinker::units::bohr;
                induced_dipole_polarization[i_i_mm * 3 + i_xyz] = all_mu_p[i_mm * 3 + i_xyz] / tinker::units::bohr;
            }
        }

        delete[] workspace_n_mm_times_6;
    }
    else
    {
        printf("TC anchor: Warning: the induced dipole is computed, but polar is not specified in Tinker parameter.\n");
        memset(induced_dipole_direct, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
        memset(induced_dipole_direct, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
    }
}

void internal_get_mm_induced_dipole(double* mu_d, double* mu_p)
{
    if (tinker::use_potent(tinker::polar_term))
    {
        const size_t n_total = tinker::n;
        tinker::real* all_mu_d = new tinker::real[n_total * 3];
        tinker::real* all_mu_p = new tinker::real[n_total * 3];
        // tinker::uin[d/p] is defined in include/mod.polar.h, allocated in src/elec.cpp::pole_data(), computed in src/acc/induce.cpp::induce_mutual_pcg1_acc().
        tinker::darray::copyout(tinker::g::q0, n_total, all_mu_d, tinker::uind);
        tinker::darray::copyout(tinker::g::q0, n_total, all_mu_p, tinker::uinp);
        tinker::wait_for(tinker::g::q0);

        for (size_t i_i_mm = 0; i_i_mm < QMMMGlobal::n_mm; i_i_mm++)
        {
            int32_t i_mm = QMMMGlobal::mm_indices[i_i_mm] - 1; // One-index to zero-index
            for (size_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                mu_d[i_i_mm * 3 + i_xyz] = all_mu_d[i_mm * 3 + i_xyz] / tinker::units::bohr;
                mu_p[i_i_mm * 3 + i_xyz] = all_mu_p[i_mm * 3 + i_xyz] / tinker::units::bohr;
            }
        }

        delete[] all_mu_d;
        delete[] all_mu_p;
    }
    else
    {
        printf("TC anchor: Warning: the induced dipole is accessed, but polar is not specified in Tinker parameter.\n");
        memset(mu_d, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
        memset(mu_p, 0, QMMMGlobal::n_mm * 3 * sizeof(double));
    }
}
