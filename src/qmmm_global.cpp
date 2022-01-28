#include "qmmm_global.h"

namespace QMMMGlobal
{
    int32_t n_qm;
    int32_t n_mm;
    int32_t* qm_indices;
    int32_t* mm_indices;

    tinker::real (*d_qmmm_electric_field_d)[3];
    tinker::real (*d_qmmm_electric_field_p)[3];
}
