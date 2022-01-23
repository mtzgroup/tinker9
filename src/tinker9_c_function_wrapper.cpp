#include <stdio.h>

#include "tinker9_c_function_wrapper.h"
#include "tinker8_fortran_function_wrapper.h"
#include "qmmm_global.h"

#include "tinker_rt.h" // initial(), routine.h::tinker_f_*()
#include "md.h" // mdcalc.h::calc::*, mdegv.h::copy_gradient()
#include "energy.h" // energy()

#include "nblist.h"

#include <tinker/detail/bndstr.hh>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/urey.hh>
#include <tinker/detail/tors.hh>

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
    
    // ! VDW, there are 2 types of them

    // ! with the next loop, we are checking that the last column and last row

    // ! In order to zero the vdW interactions, I will exploit the fact that
    // ! radmin and epsilon matrices are zeroed and no one is ever using the last row and the
    // ! last column. These matrices are created of fixed size of 1000,1000, dimension
    // ! given into file MOD_size.f.

    // ! count how many unique vdW types are in just the qm region.
    // ! Also assigns a new vdw atom type to each new atom
    // counter = 0
    // do i=1, n_qm
    //   this_number = jvdw(qm_ind(i))
    //   is_unique = .True.
    //   do j=1, i-1
    //     if (this_number .eq. jvdw(qm_ind(j))) then
    //       mapping_qm(i) = mapping_qm(j)
    //       is_unique = .False.
    //     end if
    //   end do
    //   if (is_unique) then
    //     mapping_qm(i) = maxclass - counter
    //     counter = counter + 1
    //   end if
    // end do
    // print *, "There are ", counter, " distinct atom types in the QM region"
    // print *, "after", mapping_qm


    // ! Check if epsilon matrix is symmetric: it always is
    // ! print *, "Checking if epsilon matrix is symmetric (this can be removed later)"
    // ! do i=1, maxclass
    // !   do j=1, i-1
    // !     if (epsilon(j,i) .ne. epsilon(i,j)) then
    // !       print *, "FATAL ERROR"
    // !       call exit(1)
    // !     end if
    // !   end do
    // ! end do


    // ! Checks both that the rows and columns at the very end of the epsilon
    // ! matrix that we are assigning to the qm region
    // ! are all 0, meaning all `maxclass`-counter atom types have not been used.
    // ! Additional assumption: radmin matrix has the exactly same structure as epsilon matrix.
    // sumcheck_epsilon=0
    // do j=(maxclass-counter)+1, maxclass
    //   do i=1, j
    //     sumcheck_epsilon = sumcheck_epsilon+epsilon(j,i)+epsilon(i,j)
    //   end do
    // end do
    // if (sumcheck_epsilon .ne. 0.0) then
    //   print *, "FATAL ERROR: Need to set `maxclass` higher in MOD_sizes.f"
    //   call exit(1)
    // end if

    // ! the size of n_qm is small and this is happening at initialization. So we
    // ! will use this simple loop and overwrite many things
    // ! in this version we are just looping until (maxclass-counter) because we
    // ! want zeros in the bottom right corner of the epsilon matrix
    // do i=1, n_qm
    //   do j=1, (maxclass-counter)
    //     radmin4(mapping_qm(i),j) = radmin4(jvdw(qm_ind(i)),j)
    //     radmin4(j,mapping_qm(i)) = radmin4(j,jvdw(qm_ind(i)))
    //     radmin(mapping_qm(i),j) = radmin(jvdw(qm_ind(i)),j)
    //     radmin(j,mapping_qm(i)) = radmin(j,jvdw(qm_ind(i)))
    //     epsilon4(mapping_qm(i),j) = epsilon4(jvdw(qm_ind(i)),j)
    //     epsilon4(j,mapping_qm(i)) = epsilon4(j,jvdw(qm_ind(i)))
    //     epsilon(mapping_qm(i),j) = epsilon(jvdw(qm_ind(i)),j)
    //     epsilon(j,mapping_qm(i)) = epsilon(j,jvdw(qm_ind(i)))
    //   end do
    //   ! van der waals radius
    //   rad(mapping_qm(i)) = rad(jvdw(qm_ind(i)))
    //   rad4(mapping_qm(i)) = rad4(jvdw(qm_ind(i)))

    //   atmcls(mapping_qm(i)) = atmcls(jvdw(qm_ind(i))) ! TODO: understand this
    //   atmnum(mapping_qm(i)) = atmnum(jvdw(qm_ind(i))) ! atomic number
    //   ligand(mapping_qm(i)) = ligand(jvdw(qm_ind(i))) ! num ligands
    //   symbol(mapping_qm(i)) = symbol(jvdw(qm_ind(i))) ! atomic symbol, 3 chars
    //   describe(mapping_qm(i)) = describe(jvdw(qm_ind(i))) ! string description

    //   ! TODO: Right now, we are copying over a lot of atom type information to
    //   ! the new qm atom types indiscriminately, even though we may only need
    //   ! the vdw parameters copied. We should investigate whether it is better to
    //   ! copy the rest of these or not by default

    //   ! parameters related to atomic multipoles
    //   ! TODO: maybe don't copy this over.
    //   ! TODO: This is how we can write the effective multipoles of the QM region
    //   ! after SCF
    //   do j=1,3
    //     sibfacp(j,mapping_qm(i)) = sibfacp(j,jvdw(qm_ind(i)))
    //   end do

    //   ! Polarization group connected atoms.
    //   ! TODO: maybe need to update in a more nuanced way to account for 
    //   do j=1,maxvalue
    //     pgrp(j, mapping_qm(i)) = pgrp(j, jvdw(qm_ind(i)))
    //   end do
    //   polr(mapping_qm(i)) = polr(jvdw(qm_ind(i))) ! dipole polarizability
    //   athl(mapping_qm(i)) = athl(jvdw(qm_ind(i))) ! Thole damping value

    //   ! aromatic parameters
    //   do j=1,6
    //     mmffarom(mapping_qm(i),j) = mmffarom(jvdw(qm_ind(i)),j) ! Aromatic param
    //     mmffaromc(mapping_qm(i),j) = mmffaromc(jvdw(qm_ind(i)),j) ! cationic param
    //     mmffaroma(mapping_qm(i),j) = mmffaroma(jvdw(qm_ind(i)),j) ! anionic param
    //   end do

    //   ! After copying everything related to the atom type over to new atom types
    //   ! can finally set the atom type itself and forget old atom types
    //   jvdw(qm_ind(i)) = mapping_qm(i)
    // end do

    // ! Set qm charges to 0 to remove charge-charge energy
    // if (use_charge) then
    //   chglist(qm_ind(1:n_qm)) = 0
    //   pchg(qm_ind(1:n_qm)) = 0.0
    // end if
    // if (use_mpole) then
    //   pollist(qm_ind(1:n_qm)) = 0
    //   pole(:, qm_ind(1:n_qm)) = 0.0
    // end if


    // ! Set the polarizability and Thole parameter to 0
    // polarity(qm_ind(1:n_qm)) = 0.0
    // thole(qm_ind(1:n_qm)) = 0.0

    /**
     * Stores these variables in tinker global memory.
     * They're used in src/rc_man.cpp::initialize() -> device_data() -> src/energy.cpp::energy_data() -> src/evdw.cpp::evdw_data(),
     * so must be set before initalize() function call.
     */
    QMMMGlobal::n_qm = n_qm;
    QMMMGlobal::qm_indices = new int32_t[n_qm];
    memcpy(QMMMGlobal::qm_indices, qm_indices, n_qm * sizeof(int32_t));

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
