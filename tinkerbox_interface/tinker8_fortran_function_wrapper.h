#ifndef TINKER8_FORTRAN_FUNCTION_WRAPPER_
#define TINKER8_FORTRAN_FUNCTION_WRAPPER_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void tinkerbox_getxyz_(const char* const xyzfile_c, const int32_t* const xyzfile_string_length);

#ifdef __cplusplus
}
#endif

#endif // #ifndef TINKER8_FORTRAN_FUNCTION_WRAPPER_
