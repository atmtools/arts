/* copyright (C) 2003-2012 Cory Davis <cory.davis@metservice.com>
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
 * @file   montecarlo.cc
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-19
 *
 * @brief  functions used by MCGeneral
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "montecarlo.h"

#include <workspace.h>

#include <cfloat>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "debug.h"
#include "rte.h"
#include "special_interp.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

// Some root-finding helper functions (for MCRadar) that don't need
// visibility outside this source file
//
//

/** ext_I.
 *
 * Calculate the extinction of I for a propagating "photon"
 *
 * @param[in] I  1st Stokes element *** FIXMEDOC ***
 * @param[in] Q  2nd Stokes element *** FIXMEDOC ***
 * @param[in] KI  Extinction matrix element 0,0 *** FIXMEDOC ***
 * @param[in] KQ  Extinction matrix element 0,1 *** FIXMEDOC ***
 * @param[in] s   Pathlength *** FIXMEDOC ***
 *
 * @author Ian Adams
 * @date  2016-06-15
 */
Numeric ext_I(const Numeric& I,
              const Numeric& Q,
              const Numeric& kI,
              const Numeric& kQ,
              const Numeric& s) {
  Numeric fs;

  fs = exp(-kI * s) * (I * cosh(kQ * s) + Q * sinh(kQ * s));

  return fs;
}

/** brent_zero.
 *
 *
 * This function seeks the root of a function F(X) in an interval [A,B].
 *
 * Discussion:
 * The interval [A,B] must be a change of sign interval for F. That is,
 * F(A) and F(B) must be of opposite signs. Then assuming that F is
 * continuous implies the existence of at least one value C between A and
 * B for which F(C) = 0.
 *
 * The location of the zero is determined to within an accuracy of
 * 6 * MACHEPS * abs ( C ) + 2 * T.
 *
 * Thanks to Thomas Secretin for pointing out a transcription error in the
 * setting of the value of P, 11 February 2013.
 *
 * Modifications by Ian S. Adams, U.S. Naval Research Laboratory to conform
 * to ARTS and to hardcode function for root finding while passing in
 * mulitple args for function.
 *
 * Licensing:
 *
 * This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *
 * 11 February 2013, J. Burkardt
 * 15 July 2016    , I. Adams
 *
 * Author:
 *
 * Original FORTRAN77 version by Richard Brent.
 * C++ version by John Burkardt.
 *
 * Reference:
 *
 * Richard Brent,
 * Algorithms for Minimization Without Derivatives,
 * Dover, 2002,
 * ISBN: 0-486-41998-3,
 * LC: QA402.5.B74.
 *
 * Parameters:
 *
 * Input, double A, B, the endpoints of the change of sign interval.
 * Input, double T, a positive error tolerance.
 * Output, double ZERO, the estimated value of a zero of the function F.
 *
 * @param[out] sb  The estimated value of a zero of the function F.
 * @param[in]  a   The lower endpoint of the change of sign interval
 * @param[in]  b   The upper endpoint of the change of sign interval
 * @param[in]  t   A positive error tolerance.
 * @param[in]  rn  A Random number.
 * @param[in]  I   1st Stokes element *** FIXMEDOC ***
 * @param[in]  Q   2nd Stokes element *** FIXMEDOC ***
 * @param[in]  KI  Extinction matrix element 0,0 *** FIXMEDOC ***
 * @param[in]  KQ  Extinction matrix element 0,1 *** FIXMEDOC ***
 *
 * @author J. Burkhardt
 * @date 20??-??-??
 */
void brent_zero(Numeric& sb,
                const Numeric& a,
                const Numeric& b,
                const Numeric& t,
                const Numeric& rn,
                const Numeric& I,
                const Numeric& Q,
                const Numeric& KI,
                const Numeric& KQ) {
  Numeric c;
  Numeric d;
  Numeric e;
  Numeric fa;
  Numeric fb;
  Numeric fc;
  Numeric m;
  Numeric macheps;
  Numeric p;
  Numeric q;
  Numeric r;
  Numeric s;
  Numeric sa;
  Numeric tol;
  //
  //  Make local copies of A and B.
  //
  sa = a;
  sb = b;
  fa = ext_I(I, Q, KI, KQ, sa);  // - rn;
  fa -= rn;
  fb = ext_I(I, Q, KI, KQ, sb);  // - rn;
  fb -= rn;

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;

  macheps = DBL_EPSILON;

  for (;;) {
    if (abs(fc) < abs(fb)) {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * macheps * abs(sb) + t;
    m = 0.5 * (c - sb);

    if (abs(m) <= tol || fb == 0.0) {
      break;
    }

    if (abs(e) < tol || abs(fa) <= abs(fb)) {
      e = m;
      d = e;
    } else {
      s = fb / fa;

      if (sa == c) {
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if (0.0 < p) {
        q = -q;
      } else {
        p = -p;
      }

      s = e;
      e = d;

      if (2.0 * p < 3.0 * m * q - abs(tol * q) && p < abs(0.5 * s * q)) {
        d = p / q;
      } else {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if (tol < abs(d)) {
      sb = sb + d;
    } else if (0.0 < m) {
      sb = sb + tol;
    } else {
      sb = sb - tol;
    }

    fb = ext_I(I, Q, KI, KQ, sb);
    fb -= rn;

    if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
}

void clear_rt_vars_at_gp(const Workspace& ws,
                         MatrixView ext_mat_mono,
                         VectorView abs_vec_mono,
                         Numeric& temperature,
                         const Agenda& propmat_clearsky_agenda,
                         const Numeric& f_mono,
                         const GridPos& gp_p,
                         const GridPos& gp_lat,
                         const GridPos& gp_lon,
                         ConstVectorView p_grid,
                         ConstTensor3View t_field,
                         ConstTensor4View vmr_field) {
  const Index ns = vmr_field.nbooks();
  Vector t_vec(1);  //vectors are required by interp_atmfield_gp2itw etc.
  Vector p_vec(1);  //may not be efficient with unecessary vectors
  Matrix itw_p(1, 2);
  ArrayOfGridPos ao_gp_p(1), ao_gp_lat(1), ao_gp_lon(1);
  Matrix vmr_mat(ns, 1), itw_field;

  //local versions of workspace variables
  StokvecVector local_abs_vec;
  StokvecVector local_nlte_source_dummy;
  PropmatVector local_ext_mat;
  PropmatVector local_propmat_clearsky;
  PropmatMatrix local_partial_dummy;
  StokvecMatrix local_dnlte_source_dx_dummy;
  ao_gp_p[0] = gp_p;
  ao_gp_lat[0] = gp_lat;
  ao_gp_lon[0] = gp_lon;

  // Determine the pressure
  interpweights(itw_p, ao_gp_p);
  itw2p(p_vec, p_grid, ao_gp_p, itw_p);

  // Determine the atmospheric temperature and species VMR
  //
  interp_atmfield_gp2itw(itw_field, ao_gp_p, ao_gp_lat, ao_gp_lon);
  //
  interp_atmfield_by_itw(
      t_vec, t_field, ao_gp_p, ao_gp_lat, ao_gp_lon, itw_field);
  //
  for (Index is = 0; is < ns; is++) {
    interp_atmfield_by_itw(vmr_mat(is, joker),
                           vmr_field(is, joker, joker, joker),
                           ao_gp_p,
                           ao_gp_lat,
                           ao_gp_lon,
                           itw_field);
  }

  temperature = t_vec[0];

  //calcualte absorption coefficient
  propmat_clearsky_agendaExecute(ws,
                                 local_propmat_clearsky,
                                 local_nlte_source_dummy,
                                 local_partial_dummy,
                                 local_dnlte_source_dx_dummy,
                                 {},
                                 {},
                                 Vector(1, f_mono),
                                 {},
                                 AtmPoint{},  // FIXME: DUMMY VALUE
                                 propmat_clearsky_agenda);

  opt_prop_sum_propmat_clearsky(
      local_ext_mat, local_abs_vec, local_propmat_clearsky);

  ext_mat_mono = to_matrix(local_ext_mat[0]);
  abs_vec_mono = to_vector(local_abs_vec[0]);
}

void cloudy_rt_vars_at_gp(const Workspace& ws,
                          MatrixView ext_mat_mono,
                          VectorView abs_vec_mono,
                          VectorView pnd_vec,
                          Numeric& temperature,
                          const Agenda& propmat_clearsky_agenda,
                          const Index f_index,
                          const Vector& f_grid,
                          const GridPos& gp_p,
                          const GridPos& gp_lat,
                          const GridPos& gp_lon,
                          ConstVectorView p_grid_cloud,
                          ConstTensor3View t_field_cloud,
                          ConstTensor4View vmr_field_cloud,
                          const Tensor4& pnd_field,
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const ArrayOfIndex& cloudbox_limits,
                          const Vector& rte_los)

{
  const Index ns = vmr_field_cloud.nbooks();
  const Index N_se = pnd_field.nbooks();
  Matrix pnd_ppath(N_se, 1);
  Vector t_ppath(1);
  Vector p_ppath(1);  //may not be efficient with unecessary vectors
  Matrix itw_p(1, 2);
  ArrayOfGridPos ao_gp_p(1), ao_gp_lat(1), ao_gp_lon(1);
  Matrix vmr_ppath(ns, 1), itw_field;

  //local versions of workspace variables
  PropmatMatrix
      local_partial_dummy;  // This is right since there should be only clearsky partials
  StokvecMatrix local_dnlte_source_dx_dummy;  // This is right since there should be only clearsky partials
  PropmatVector local_propmat_clearsky;
  StokvecVector local_nlte_source_dummy;  //FIXME: Do this right?
  StokvecVector local_abs_vec;
  PropmatVector local_ext_mat;

  ao_gp_p[0] = gp_p;
  ao_gp_lat[0] = gp_lat;
  ao_gp_lon[0] = gp_lon;

  cloud_atm_vars_by_gp(p_ppath,
                       t_ppath,
                       vmr_ppath,
                       pnd_ppath,
                       ao_gp_p,
                       ao_gp_lat,
                       ao_gp_lon,
                       cloudbox_limits,
                       p_grid_cloud,
                       t_field_cloud,
                       vmr_field_cloud,
                       pnd_field);
  pnd_vec = pnd_ppath(joker, 0);
  temperature = t_ppath[0];

  //rtp_vmr    = vmr_ppath(joker,0);
  propmat_clearsky_agendaExecute(ws,
                                 local_propmat_clearsky,
                                 local_nlte_source_dummy,
                                 local_partial_dummy,
                                 local_dnlte_source_dx_dummy,
                                 {},
                                 {},
                                 Vector{f_grid[Range(f_index, 1)]},
                                 {},
                                 AtmPoint{},  // FIXME: DUMMY VALUE
                                 propmat_clearsky_agenda);

  opt_prop_sum_propmat_clearsky(
      local_ext_mat, local_abs_vec, local_propmat_clearsky);

  ext_mat_mono = to_matrix(local_ext_mat[0]);
  abs_vec_mono = to_vector(local_abs_vec[0]);

  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk;
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  Vector sca_dir;
  mirror_los(sca_dir, rte_los);
  Matrix dir_array(1, 2, 0.);
  dir_array(0, joker) = sca_dir;
  //
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      t_ppath,
                      dir_array,
                      f_index);
  //
  opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                        abs_vec_ssbulk,
                        ptype_ssbulk,
                        ext_mat_Nse,
                        abs_vec_Nse,
                        ptypes_Nse,
                        pnd_ppath,
                        t_ok);
  opt_prop_Bulk(ext_mat_bulk,
                abs_vec_bulk,
                ptype_bulk,
                ext_mat_ssbulk,
                abs_vec_ssbulk,
                ptype_ssbulk);

  ext_mat_mono += ext_mat_bulk(0, 0, 0, joker, joker);
  abs_vec_mono += abs_vec_bulk(0, 0, 0, joker);
}

void cloud_atm_vars_by_gp(VectorView pressure,
                          VectorView temperature,
                          MatrixView vmr,
                          MatrixView pnd,
                          const ArrayOfGridPos& gp_p,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstVectorView p_grid_cloud,
                          ConstTensor3View t_field_cloud,
                          ConstTensor4View vmr_field_cloud,
                          ConstTensor4View pnd_field) {
  Index np = gp_p.size();
  ARTS_ASSERT(pressure.size() == np);
  Index ns = vmr_field_cloud.nbooks();
  Index N_se = pnd_field.nbooks();
  ArrayOfGridPos gp_p_cloud = gp_p;
  ArrayOfGridPos gp_lat_cloud = gp_lat;
  ArrayOfGridPos gp_lon_cloud = gp_lon;

  for (Index i = 0; i < np; i++) {
    // Calculate grid positions inside the cloud.
    gp_p_cloud[i].idx -= cloudbox_limits[0];
    gp_lat_cloud[i].idx -= cloudbox_limits[2];
    gp_lon_cloud[i].idx -= cloudbox_limits[4];
  }

  const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
  const Index n2 = cloudbox_limits[3] - cloudbox_limits[2];
  const Index n3 = cloudbox_limits[5] - cloudbox_limits[4];
  gridpos_upperend_check(gp_p_cloud[0], n1);
  gridpos_upperend_check(gp_p_cloud[np - 1], n1);
  gridpos_upperend_check(gp_lat_cloud[0], n2);
  gridpos_upperend_check(gp_lat_cloud[np - 1], n2);
  gridpos_upperend_check(gp_lon_cloud[0], n3);
  gridpos_upperend_check(gp_lon_cloud[np - 1], n3);

  // Determine the pressure at each propagation path point
  Matrix itw_p(np, 2);
  //
  //interpweights( itw_p, ppath.gp_p );
  interpweights(itw_p, gp_p_cloud);
  itw2p(pressure, p_grid_cloud, gp_p_cloud, itw_p);

  // Determine the atmospheric temperature and species VMR at
  // each propagation path point
  Matrix itw_field;
  //
  interp_atmfield_gp2itw(
      itw_field, gp_p_cloud, gp_lat_cloud, gp_lon_cloud);
  //
  interp_atmfield_by_itw(temperature,
                         t_field_cloud,
                         gp_p_cloud,
                         gp_lat_cloud,
                         gp_lon_cloud,
                         itw_field);
  //
  for (Index is = 0; is < ns; is++) {
    interp_atmfield_by_itw(vmr(is, joker),
                           vmr_field_cloud(is, joker, joker, joker),
                           gp_p_cloud,
                           gp_lat_cloud,
                           gp_lon_cloud,
                           itw_field);
  }

  //Determine the particle number density for every scattering element at
  // each propagation path point
  for (Index i_se = 0; i_se < N_se; i_se++) {
    // if grid positions still outside the range the propagation path step
    // must be outside the cloudbox and pnd is set to zero
    interp_atmfield_by_itw(pnd(i_se, joker),
                           pnd_field(i_se, joker, joker, joker),
                           gp_p_cloud,
                           gp_lat_cloud,
                           gp_lon_cloud,
                           itw_field);
  }
}

void ext_mat_case(Index& icase,
                  ConstMatrixView ext_mat) {
  if (icase == 0) {
    icase = 1;  // Start guess is diagonal

    {
      // Check symmetries and analyse structure of exp_mat:
      ARTS_ASSERT(ext_mat(1, 1) == ext_mat(0, 0));
      ARTS_ASSERT(ext_mat(1, 0) == ext_mat(0, 1));

      if (ext_mat(1, 0) != 0) {
        icase = 2;
      }

      {
        ARTS_ASSERT(ext_mat(2, 2) == ext_mat(0, 0));
        ARTS_ASSERT(ext_mat(2, 1) == -ext_mat(1, 2));
        ARTS_ASSERT(ext_mat(2, 0) == ext_mat(0, 2));

        if (ext_mat(2, 0) != 0 || ext_mat(2, 1) != 0) {
          icase = 3;
        }

        {
          ARTS_ASSERT(ext_mat(3, 3) == ext_mat(0, 0));
          ARTS_ASSERT(ext_mat(3, 2) == -ext_mat(2, 3));
          ARTS_ASSERT(ext_mat(3, 1) == -ext_mat(1, 3));
          ARTS_ASSERT(ext_mat(3, 0) == ext_mat(0, 3));

          if (icase < 3)  // if icase==3, already at most complex case
          {
            if (ext_mat(3, 0) != 0 || ext_mat(3, 1) != 0) {
              icase = 3;
            } else if (ext_mat(3, 2) != 0) {
              icase = 2;
            }
          }
        }
      }
    }
  }
}

/** Converts an extinction matrix to a transmission matrix

    The function performs the calculations differently depending on the
    conditions, to improve the speed. There are three cases: <br>
       1. Scalar RT and/or the matrix ext_mat_av is diagonal. <br>
       2. Special expression for "azimuthally_random" case. <br>
       3. The total general case.

    If the structure of *ext_mat* is known, *icase* can be set to "case index"
    (1, 2 or 3) and some time is saved. This includes that no asserts are
    performed on *ext_mat*.

    Otherwise, *icase* must be set to 0. *ext_mat* is then analysed and *icase*
    is set by the function and is returned.

    trans_mat must be sized before calling the function.

    @param[out]   trans_mat      Transmission matrix of slab.
    @param[out]   icase          Index giving ext_mat case.
    @param[in]    ext_mat        Averaged extinction matrix.
    @param[in]    lstep          The length of the RTE step.

    @author Patrick Eriksson (based on earlier version started by Claudia)
    @date   2013-05-17 
 */
void ext2trans(MatrixView trans_mat,
               Index& icase,
               ConstMatrixView ext_mat,
               const Numeric& lstep) {
  ARTS_ASSERT(ext_mat.nrows() == 4);
  ARTS_ASSERT(trans_mat.nrows() == 4 && trans_mat.ncols() == 4);

  // Theoretically ext_mat(0,0) >= 0, but to demand this can cause problems for
  // iterative retrievals, and the assert is skipped. Negative should be a
  // result of negative vmr, and this issue is checked in basics_checkedCalc.
  //ARTS_ASSERT( ext_mat(0,0) >= 0 );

  ARTS_ASSERT(icase >= 0 && icase <= 3);
  ARTS_ASSERT(!is_singular(ext_mat));
  ARTS_ASSERT(lstep >= 0);

  // Analyse ext_mat?
  ext_mat_case(icase, ext_mat);

  // Calculation options:
  if (icase == 1) {
    trans_mat = 0;
    trans_mat(0, 0) = exp(-ext_mat(0, 0) * lstep);
    for (Index i = 1; i < 4; i++) {
      trans_mat(i, i) = trans_mat(0, 0);
    }
  }
  
  else {
    Matrix ext_mat_ds{ext_mat};
    ext_mat_ds *= -lstep;
    //
    constexpr Index q = 10;  // index for the precision of the matrix exp function
    matrix_exp(trans_mat, ext_mat_ds, q);
  }
}

bool is_anyptype_nonTotRan(
    const ArrayOfArrayOfSingleScatteringData& scat_data) {
  bool is_anyptype_nonTotRan = false;
  for (Size i_ss = 0;
       is_anyptype_nonTotRan == false && i_ss < scat_data.size();
       i_ss++) {
    for (Size i_se = 0;
         is_anyptype_nonTotRan == false && i_se < scat_data[i_ss].size();
         i_se++) {
      if (scat_data[i_ss][i_se].ptype > PTYPE_TOTAL_RND) {
        is_anyptype_nonTotRan = true;
      }
    }
  }
  return is_anyptype_nonTotRan;
}

void Sample_los(VectorView new_rte_los,
                Numeric& g_los_csc_theta,
                MatrixView Z,
                RandomNumberGenerator<>& rng,
                ConstVectorView rte_los,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Index f_index,
                ConstVectorView pnd_vec,
                ConstVectorView Z11maxvector,
                const Numeric Csca,
                const Numeric rtp_temperature,
                const Index t_interp_order) {
  Numeric Z11max = 0;
  bool tryagain = true;

  Vector sca_dir;
  mirror_los(sca_dir, rte_los);

  // Rejection method http://en.wikipedia.org/wiki/Rejection_sampling
  Index np = pnd_vec.size();
  ARTS_ASSERT(TotalNumberOfElements(scat_data) == np);
  for (Index i = 0; i < np; i++) {
    Z11max += Z11maxvector[i] * pnd_vec[i];
  }

  ///////////////////////////////////////////////////////////////////////
  // allocating variables needed for pha_mat extraction - this seems a little
  // disadvantageous. If we move the whole function back into the calling one
  // (it's used only at one place anyways), we can do the allocation there once
  // and avoid that it is done everytime this function is called within a loop.
  ArrayOfArrayOfTensor6 pha_mat_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 pha_mat_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 pha_mat_bulk;
  Index ptype_bulk;
  Matrix pdir(1, 2), idir(1, 2);
  Vector t(1, rtp_temperature);
  Matrix pnds(np, 1);
  pnds(joker, 0) = pnd_vec;

  while (tryagain) {
    new_rte_los[0] = Conversion::acosd(rng.get(-1.0, 1.0)());
    new_rte_los[1] = rng.get(-180.0, 180.0)();

    //Calculate Phase matrix////////////////////////////////
    Vector inc_dir;
    mirror_los(inc_dir, new_rte_los);

    //pha_mat_singleCalc( Z, sca_dir[0], sca_dir[1], inc_dir[0], inc_dir[1],
    //                    scat_data_mono, pnd_vec, rtp_temperature,
    //                     );

    pdir(0, joker) = sca_dir;
    idir(0, joker) = inc_dir;
    pha_mat_NScatElems(pha_mat_Nse,
                       ptypes_Nse,
                       t_ok,
                       scat_data,
                       t,
                       pdir,
                       idir,
                       f_index,
                       t_interp_order);
    pha_mat_ScatSpecBulk(
        pha_mat_ssbulk, ptype_ssbulk, pha_mat_Nse, ptypes_Nse, pnds, t_ok);
    pha_mat_Bulk(pha_mat_bulk, ptype_bulk, pha_mat_ssbulk, ptype_ssbulk);
    Z = pha_mat_bulk(0, 0, 0, 0, joker, joker);

    if (rng.get(0.0, 1.0)() <= Z(0, 0) / Z11max)  //then new los is accepted
    {
      tryagain = false;
    }
  }
  g_los_csc_theta = Z(0, 0) / Csca;
}

void Sample_los_uniform(VectorView new_rte_los, RandomNumberGenerator<>& rng) {
  new_rte_los[1] = rng.get<>(-180., 180.)();
  new_rte_los[0] = Conversion::acosd(rng.get<>(-1., 1.)());
}
