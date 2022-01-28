#ifndef QMMM_GLOBAL_H_
#define QMMM_GLOBAL_H_

#include "macro.h" // for real
#include <stdint.h>

namespace QMMMGlobal
{
    extern int32_t n_qm;
    extern int32_t n_mm;
    extern int32_t* qm_indices;
    extern int32_t* mm_indices;
    
    TINKER_EXTERN tinker::real (*d_qmmm_electric_field_d)[3];
    TINKER_EXTERN tinker::real (*d_qmmm_electric_field_p)[3];
}

#endif
