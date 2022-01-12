#include <stdio.h>

#include "tinker_rt.h" // initial(), routine.h::tinker_f_*()
#include "md.h" // mdcalc.h::calc::*, mdegv.h::copy_gradient()
#include "energy.h" // energy()

#include "nblist.h"

#include "tinker8_fortran_function_wrapper.h"
#include "tinker9_c_function_wrapper.h"

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
    



    printf("Henry: TODO: remove parameters for QM atoms\n"); fflush(stdout);
    
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
    printf("Henry: TODO: make sure this is compatable with Terachem\n"); fflush(stdout);
    tinker::initialize();
    
#if TINKER_CUDART
    printf("Checking GPU compilation: Tinker9 compiled with CUDA\n");
    if (tinker::vlist_version() & tinker::NBL_SPATIAL)
        printf("Checking GPU runtime: vdw_term running with OpenACC\n");
    else
        printf("Checking GPU runtime: vdw_term running with OpenACC\n");
    if (tinker::clist_version() & tinker::NBL_SPATIAL)
        printf("Checking GPU runtime: charge_term running with OpenACC\n");
    else
        printf("Checking GPU runtime: charge_term running with OpenACC\n");
    if (tinker::mlist_version() & tinker::NBL_SPATIAL)
        printf("Checking GPU runtime: mpole_term, polar_term, chgtrn_term, repuls_term running with OpenACC\n");
    else
        printf("Checking GPU runtime: mpole_term, polar_term, chgtrn_term, repuls_term running with OpenACC\n");
    if (tinker::ulist_version() & tinker::NBL_SPATIAL)
        printf("Checking GPU runtime: polar_term running with OpenACC\n");
    else
        printf("Checking GPU runtime: polar_term running with OpenACC\n");
#else
    printf("Checking GPU compilation: Tinker9 compiled with only OpenACC\n");
    printf("Checking GPU runtime: All running with OpenACC\n");
#endif
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
