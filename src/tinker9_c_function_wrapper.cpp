#include <stdio.h>

#include "tinker9_c_function_wrapper.h"
#include "tinker8_fortran_function_wrapper.h"
#include "qmmm_global.h"

#include "tinker_rt.h" // initial(), routine.h::tinker_f_*()
#include "md.h" // mdcalc.h::calc::*, mdegv.h::copy_gradient()
#include "energy.h" // energy()
#include "potent.h" // use_potent(), *_term
#include "nblist.h"

#include <tinker/detail/bndstr.hh>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/urey.hh>
#include <tinker/detail/tors.hh>
#include <tinker/detail/mpole.hh>
#include <tinker/detail/polar.hh>

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

int32_t internal_get_n_mm()
{
    return tinker::n;
}

double internal_get_energy_nonpolar_mm_contribution()
{
    tinker::energy(tinker::rc_flag);
    return tinker::esum;
}

void internal_get_gradients_all_atoms_mm_contribution(double* grad)
{
    std::vector<double> gdx(tinker::n), gdy(tinker::n), gdz(tinker::n);
    tinker::copy_gradient(tinker::calc::grad, gdx.data(), gdy.data(), gdz.data());

    for (size_t i_atom = 0; i_atom < tinker::n; i_atom++)
    {
        grad[i_atom * 3 + 0] = gdx[i_atom];
        grad[i_atom * 3 + 1] = gdy[i_atom];
        grad[i_atom * 3 + 2] = gdz[i_atom];
    }
}
