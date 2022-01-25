#include "echglj.h"
#include "box.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "qmmm_global.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include "tool/io_text.h"
#include <array>
#include <map>
#include <set>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/kvdws.hh>
#include <tinker/detail/params.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/vdwpot.hh>


namespace tinker {
void echglj_data(rc_op op)
{
   if (!(clist_version() & NBL_SPATIAL))
      return;


   if (!use_potent(charge_term) || !use_potent(vdw_term))
      return;


   if (op & rc_dealloc) {
      ncvexclude = 0;
      vdwpr_in_use = false;
      darray::deallocate(cvexclude, cvexclude_scale);
      darray::deallocate(atom_rad, atom_eps);
      darray::deallocate(mut_coalesced, chg_coalesced, radeps_coalesced);
   }


   if (op & rc_alloc) {
      // check "VDWPR" keyword
      if (!vdwpr_in_use) {
         auto parse_vdwpr = [](std::string line, int& i, int& k, double& rad,
                               double& eps) {
            try {
               auto vs = Text::split(line);
               std::string ke = vs.at(0);
               Text::upcase(ke);
               if (ke == "VDWPR") {
                  i = std::stoi(vs.at(1));
                  k = std::stoi(vs.at(2));
                  rad = std::stod(vs.at(3));
                  eps = std::stod(vs.at(4));
                  return true;
               }
               return false;
            } catch (...) {
               return false;
            }
         };
         const int* src = atoms::type;
         if (vdwindex == evdw_t::atom_class)
            src = atomid::class_;
         std::set<int> all_tid(src, src + n);
         auto end = all_tid.end();
         int i, k;
         double rad, eps;
         std::string record;
         // first prm
         for (int ii = 0; ii < params::nprm && !vdwpr_in_use; ++ii) {
            fstr_view fsv = params::prmline[ii];
            record = fsv.trim();
            bool okay = parse_vdwpr(record, i, k, rad, eps);
            if (okay) {
               auto iit = all_tid.find(i);
               auto kit = all_tid.find(k);
               if (iit != end && kit != end) {
                  vdwpr_in_use = true;
               }
            }
         }
         // then key
         for (int ii = 0; ii < keys::nkey && !vdwpr_in_use; ++ii) {
            fstr_view fsv = keys::keyline[ii];
            record = fsv.trim();
            bool okay = parse_vdwpr(record, i, k, rad, eps);
            if (okay) {
               auto iit = all_tid.find(i);
               auto kit = all_tid.find(k);
               if (iit != end && kit != end) {
                  vdwpr_in_use = true;
               }
            }
         }
      }


      struct cv
      {
         real c, v;
      };
      auto insert_cv = [](std::map<std::pair<int, int>, cv>& a, int i, int k,
                          real val, char ch) {
         std::pair<int, int> key;
         key.first = i;
         key.second = k;
         auto it = a.find(key);
         if (it == a.end()) {
            cv x;
            x.c = 1;
            x.v = 1;
            if (ch == 'c')
               x.c = val;
            else if (ch == 'v')
               x.v = val;
            a[key] = x;
         } else {
            if (ch == 'c')
               it->second.c = val;
            else if (ch == 'v')
               it->second.v = val;
         }
      };


      // see also attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      const real c2scale = chgpot::c2scale;
      const real c3scale = chgpot::c3scale;
      const real c4scale = chgpot::c4scale;
      const real c5scale = chgpot::c5scale;
      const real v2scale = vdwpot::v2scale;
      const real v3scale = vdwpot::v3scale;
      const real v4scale = vdwpot::v4scale;
      const real v5scale = vdwpot::v5scale;


      std::map<std::pair<int, int>, cv> ik_scale;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;


         // c
         if (c2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c2scale, 'c');
               }
            }
         }
         if (c3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c3scale, 'c');
               }
            }
         }
         if (c4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c4scale, 'c');
               }
            }
         }
         if (c5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c5scale, 'c');
               }
            }
         }


         // v
         if (v2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v2scale, 'v');
               }
            }
         }
         if (v3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v3scale, 'v');
               }
            }
         }
         if (v4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v4scale, 'v');
               }
            }
         }
         if (v5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v5scale, 'v');
               }
            }
         }
      }

      if (QMMMGlobal::n_qm > 0)
      {
         printf("TC anchor: Removing VDW interaction between QM atoms.\n");
         printf("\n\nWarning!\nWarning!\nWarning!\nHenry 20220123 don't know how to run this part of the code, so this piece of code is not tested!\n\n\n\n");
      }
      // Henry 20220121: See comments in src/evdw.cpp. These are the 2 files that set up the exclusion list.
      for (size_t i_i_qm = 0; i_i_qm < QMMMGlobal::n_qm; i_i_qm++)
      {
         int32_t i_qm = QMMMGlobal::qm_indices[i_i_qm] - 1; // From one-indexed to zero-indexed

         for (size_t j_i_qm = 0; j_i_qm < QMMMGlobal::n_qm; j_i_qm++)
         {
            int32_t j_qm = QMMMGlobal::qm_indices[j_i_qm] - 1; // From one-indexed to zero-indexed

            if (j_qm > i_qm)
               insert_cv(ik_scale, i_qm, j_qm, 0.0, 'v'); // It handles appending-or-replacing.
         }
      }

      std::vector<int> ik_vec;
      std::vector<real> scal_vec;
      for (auto& it : ik_scale) {
         ik_vec.push_back(it.first.first);
         ik_vec.push_back(it.first.second);
         scal_vec.push_back(it.second.c);
         scal_vec.push_back(it.second.v);
      }
      ncvexclude = ik_scale.size();
      darray::allocate(ncvexclude, &cvexclude, &cvexclude_scale);
      darray::copyin(g::q0, ncvexclude, cvexclude, ik_vec.data());
      darray::copyin(g::q0, ncvexclude, cvexclude_scale, scal_vec.data());
      wait_for(g::q0);
      darray::allocate(n, &atom_rad, &atom_eps);
      darray::allocate(n, &mut_coalesced);
      darray::allocate(n, &chg_coalesced);
      darray::allocate(2 * n, &radeps_coalesced);
   }


   if (op & rc_init) {
      std::vector<real> vrad(n), veps(n);
      if (vdwindex == evdw_t::atom_type) {
         for (int i = 0; i < n; ++i) {
            int jj = atoms::type[i] - 1;
            vrad[i] = kvdws::rad[jj];
            veps[i] = kvdws::eps[jj];
         }
      } else if (vdwindex == evdw_t::atom_class) {
         for (int i = 0; i < n; ++i) {
            int jj = atomid::class_[i] - 1;
            vrad[i] = kvdws::rad[jj];
            veps[i] = kvdws::eps[jj];
         }
      }
      darray::copyin(g::q0, n, atom_rad, vrad.data());
      darray::copyin(g::q0, n, atom_eps, veps.data());
      wait_for(g::q0);
   }


#if TINKER_CUDART
   echglj_data_cu(op);
#endif
}


void echglj(int vers)
{
#if TINKER_CUDART
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;


   host_zero(energy_ec, virial_ec);
   host_zero(energy_ev, virial_ev);


   assert(vdwtyp == evdw_t::lj);
   assert(radrule == evdw_t::arithmetic);
   assert(epsrule == evdw_t::geometric);
   if (use_ewald()) {
      echarge_ewald_recip_self(vers);
      echglj_rad_arith_eps_geom_ewald_real_cu(vers);
   } else {
      echglj_rad_arith_eps_geom_nonewald_cu(vers);
   }


   if (do_e) {
      if (elrc_vol != 0) {
         energy_prec corr = elrc_vol / volbox();
         energy_ev += corr;
         energy_vdw += corr;
      }
   }
   if (do_v) {
      if (vlrc_vol != 0) {
         virial_prec term = vlrc_vol / volbox();
         virial_ev[0] += term; // xx
         virial_ev[4] += term; // yy
         virial_ev[8] += term; // zz
         virial_vdw[0] += term;
         virial_vdw[4] += term;
         virial_vdw[8] += term;
      }
   }
#endif
}


int* mut_coalesced;
real* chg_coalesced;
real* radeps_coalesced;
grad_prec* gx_coalesced;
grad_prec* gy_coalesced;
grad_prec* gz_coalesced;
}
