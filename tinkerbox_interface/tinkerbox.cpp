#include <stdio.h>

#include "tinker_rt.h" // initial()

#include "tinker8_fortran_function_wrapper.h"
#include "tinkerbox.h"

void initialize_tinker(const int32_t* const qm_indices, const int32_t n_qm, const char* const xyzfile)
{
    /** Implemented in src/initial.cpp
     *  This function initialize a bunch of global variables to zero-start.
     *  It also prints tinker header message (TINKER9_PROMO_STRING).
     **/
    tinker::initial();

    const int32_t xyzfile_string_length = strlen(xyzfile);
    tinkerbox_getxyz_(xyzfile, &xyzfile_string_length);

    printf("Done!\n"); fflush(stdout);
}
