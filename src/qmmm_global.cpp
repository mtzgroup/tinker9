#include "qmmm_global.h"

namespace QMMMGlobal
{
    int32_t n_qm;
    int32_t n_mm;
    int32_t* qm_indices;
    int32_t* mm_indices;
    bool if_replace_electric_field_for_compute_induced_dipole;
    bool if_replace_induced_dipole_for_compute_polar_gradient;

    tinker::real (*d_qmmm_electric_field_d)[3];
    tinker::real (*d_qmmm_electric_field_p)[3];

    tinker::real (*d_qmmm_uind_temporary)[3];
    tinker::real (*d_qmmm_uinp_temporary)[3];
}
