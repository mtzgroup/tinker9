#include <stdio.h>

#include "tinker_rt.h" // initial(), routine.h::tinker_f_*()

#include "tinker8_fortran_function_wrapper.h"
#include "tinkerbox.h"

void initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile)
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
    mechanic2();
    



    printf("Henry: TODO: remove parameters for QM atoms\n"); fflush(stdout);
    
    /**
     * Copied from src/testgrad_x.cpp::x_testgrad()
     * rc_flag defined in include/mdpq.h, calc namespace defined in include/mdcalc.h
     */
    int flags = (calc::xyz + calc::mass);
    flags += (calc::energy + calc::grad);
#if TINKER_TESTGRAD_VIRIAL
    flags += calc::virial;
#endif
    rc_flag = flags;
    
    /**
     * Implemented in src/rc_man.cpp
     * Initialize Tinker device resource manager.
     * Specifically, it calls src/cudart/gpu_card.cpp::gpu_card_data() to query GPU information,
     * then allocate GPU memory in subsequent rc_man constructor calls.
     */
    printf("Henry: TODO: make sure this is compatable with Terachem\n"); fflush(stdout);
    initialize();
}

    // energy(rc_flag);

    // std::vector<double> gdx(n), gdy(n), gdz(n);
    // copy_gradient(calc::grad, gdx.data(), gdy.data(), gdz.data());