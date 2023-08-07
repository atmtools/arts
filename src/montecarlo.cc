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

#include <cfloat>
#include <sstream>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include <workspace.h>
#include "debug.h"
#include "geodetic.h"
#include "mc_interp.h"
#include "montecarlo.h"

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

  const Vector rtp_mag_dummy(3, 0);
  const Vector ppath_los_dummy;

  //calcualte absorption coefficient
  propmat_clearsky_agendaExecute(ws,
                                 local_propmat_clearsky,
                                 local_nlte_source_dummy,
                                 local_partial_dummy,
                                 local_dnlte_source_dx_dummy,
                                 ArrayOfRetrievalQuantity(0),
                                 {},
                                 Vector(1, f_mono),
                                 ppath_los_dummy,
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

  const Vector rtp_mag_dummy(3, 0);
  const Vector ppath_los_dummy;

  //rtp_vmr    = vmr_ppath(joker,0);
  propmat_clearsky_agendaExecute(ws,
                                 local_propmat_clearsky,
                                 local_nlte_source_dummy,
                                 local_partial_dummy,
                                 local_dnlte_source_dx_dummy,
                                 ArrayOfRetrievalQuantity(0),
                                 {},
                                 Vector{f_grid[Range(f_index, 1)]},
                                 ppath_los_dummy,
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
  Index np = gp_p.nelem();
  ARTS_ASSERT(pressure.nelem() == np);
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

void get_ppath_transmat(const Workspace& ws,
                        MatrixView trans_mat,
                        const Ppath& ppath,
                        const Agenda& propmat_clearsky_agenda,
                        const Index f_index,
                        const Vector& f_grid,
                        const AtmField& atm_field,
                        const ArrayOfIndex& cloudbox_limits,
                        const Tensor4& pnd_field,
                        const ArrayOfArrayOfSingleScatteringData& scat_data) {
  bool inside_cloud;
  const Index np = ppath.np;
  ArrayOfMatrix ext_matArray(2);
  ArrayOfMatrix trans_matArray(2);
  Index N_se = pnd_field.nbooks();  //Number of scattering elements
  Vector pnd_vec(N_se);
  Vector abs_vec_mono(4);
  Matrix ext_mat(4, 4);
  Matrix ext_mat_mono(4, 4);
  Matrix incT(4, 4, 0.0);
  Numeric temperature;
  Numeric dl = -999;

  id_mat(trans_mat);

  if (np > 1) {
    // range defining cloudbox
    Range p_range(cloudbox_limits[0],
                  cloudbox_limits[1] - cloudbox_limits[0] + 1);
    Range lat_range(cloudbox_limits[2],
                    cloudbox_limits[3] - cloudbox_limits[2] + 1);
    Range lon_range(cloudbox_limits[4],
                    cloudbox_limits[5] - cloudbox_limits[4] + 1);

    inside_cloud = is_gp_inside_cloudbox(ppath.gp_p[np - 1],
                                         ppath.gp_lat[np - 1],
                                         ppath.gp_lon[np - 1],
                                         cloudbox_limits,
                                         0);
    if (inside_cloud) {
      /*  FIXME: CLOUD STUFF
      cloudy_rt_vars_at_gp(ws,
                           ext_mat_mono,
                           abs_vec_mono,
                           pnd_vec,
                           temperature,
                           propmat_clearsky_agenda,
                           f_index,
                           f_grid,
                           ppath.gp_p[np - 1],
                           ppath.gp_lat[np - 1],
                           ppath.gp_lon[np - 1],
                           p_grid[p_range],
                           t_field(p_range, lat_range, lon_range),
                           vmr_field(joker, p_range, lat_range, lon_range),
                           pnd_field,
                           scat_data,
                           cloudbox_limits,
                           Vector{ppath.los(np - 1, joker)});
                           */
    } else {
      /*  FIXME: CLOUD STUFF
      clear_rt_vars_at_gp(ws,
                          ext_mat_mono,
                          abs_vec_mono,
                          temperature,
                          propmat_clearsky_agenda,
                          f_grid[f_index],
                          ppath.gp_p[np - 1],
                          ppath.gp_lat[np - 1],
                          ppath.gp_lon[np - 1],
                          p_grid,
                          t_field,
                          vmr_field);
                          */
      pnd_vec = 0.0;
    }

    trans_matArray[1] = trans_mat;
    ext_matArray[1] = ext_mat_mono;

    // Index in ppath of end point considered presently
    for (Index ip = np - 2; ip >= 0; ip--) {
      dl = ppath.lstep[ip];

      ext_matArray[0] = ext_matArray[1];
      trans_matArray[0] = trans_matArray[1];

      inside_cloud = is_gp_inside_cloudbox(ppath.gp_p[ip],
                                           ppath.gp_lat[ip],
                                           ppath.gp_lon[ip],
                                           cloudbox_limits,
                                           0);
      if (inside_cloud) {
      /*  FIXME: CLOUD STUFF
        cloudy_rt_vars_at_gp(ws,
                             ext_mat_mono,
                             abs_vec_mono,
                             pnd_vec,
                             temperature,
                             propmat_clearsky_agenda,
                             f_index,
                             f_grid,
                             ppath.gp_p[ip],
                             ppath.gp_lat[ip],
                             ppath.gp_lon[ip],
                             p_grid[p_range],
                             t_field(p_range, lat_range, lon_range),
                             vmr_field(joker, p_range, lat_range, lon_range),
                             pnd_field,
                             scat_data,
                             cloudbox_limits,
                             Vector{ppath.los(ip, joker)});
                             */
      } else {
      /*  FIXME: CLOUD STUFF
        clear_rt_vars_at_gp(ws,
                            ext_mat_mono,
                            abs_vec_mono,
                            temperature,
                            propmat_clearsky_agenda,
                            f_grid[f_index],
                            ppath.gp_p[ip],
                            ppath.gp_lat[ip],
                            ppath.gp_lon[ip],
                            p_grid,
                            t_field,
                            vmr_field);
                            */
        pnd_vec = 0.0;
      }

      ext_matArray[1] = ext_mat_mono;
      ext_mat = ext_matArray[0];
      ext_mat += ext_matArray[1];  // Factor 2 fixed by using dl/2
      //
      {
        Index extmat_case = 0;
        ext2trans(incT, extmat_case, ext_mat, dl / 2);
      }
      //
      mult(trans_mat, incT, trans_matArray[1]);
      trans_matArray[1] = trans_mat;

    }  // for( ip... )
  }    // if( np > 1 )
}

bool is_anyptype_nonTotRan(
    const ArrayOfArrayOfSingleScatteringData& scat_data) {
  bool is_anyptype_nonTotRan = false;
  for (Index i_ss = 0;
       is_anyptype_nonTotRan == false && i_ss < scat_data.nelem();
       i_ss++) {
    for (Index i_se = 0;
         is_anyptype_nonTotRan == false && i_se < scat_data[i_ss].nelem();
         i_se++) {
      if (scat_data[i_ss][i_se].ptype > PTYPE_TOTAL_RND) {
        is_anyptype_nonTotRan = true;
      }
    }
  }
  return is_anyptype_nonTotRan;
}

void mcPathTraceGeneral(const Workspace& ws,
                        MatrixView evol_op,
                        Vector& abs_vec_mono,
                        Numeric& temperature,
                        MatrixView ext_mat_mono,
                        RandomNumberGenerator<>& rng,
                        Vector& rte_pos,
                        Vector& rte_los,
                        Vector& pnd_vec,
                        Numeric& g,
                        Ppath& ppath_step,
                        Index& termination_flag,
                        bool& inside_cloud,
                        const Agenda& ppath_step_agenda,
                        const Numeric& ppath_lmax,
                        const Numeric& ppath_lraytrace,
                        const Numeric& taustep_limit,
                        const Agenda& propmat_clearsky_agenda,
                        const Index f_index,
                        const Vector& f_grid,
                        const SurfaceField& surface_field,
                        const AtmField& atm_field,
                        const ArrayOfIndex& cloudbox_limits,
                        const Tensor4& pnd_field,
                        const ArrayOfArrayOfSingleScatteringData& scat_data) {
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2);
  ArrayOfVector pnd_vecArray(2);
  Matrix ext_mat(4, 4);
  Matrix incT(4, 4, 0.0);
  Vector tArray(2);
  Matrix T(4, 4);
  Numeric k;
  Numeric ds, dl = -999;
  Index istep = 0;  // Counter for number of steps
  Matrix old_evol_op(4, 4);

  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1] = evol_op;

  // range defining cloudbox
  Range p_range(cloudbox_limits[0],
                cloudbox_limits[1] - cloudbox_limits[0] + 1);
  Range lat_range(cloudbox_limits[2],
                  cloudbox_limits[3] - cloudbox_limits[2] + 1);
  Range lon_range(cloudbox_limits[4],
                  cloudbox_limits[5] - cloudbox_limits[4] + 1);

  //initialise Ppath with ppath_start_stepping
  ARTS_USER_ERROR("NOT PORTED ppath_start_stepping TO USE ATMFIELD")
  /*
  ppath_start_stepping(ppath_step,
                       3,
                       p_grid,
                       lat_grid,
                       lon_grid,
                       z_field,
                       refellipsoid,
                       z_surface,
                       0,
                       cloudbox_limits,
                       false,
                       rte_pos,
                       rte_los,
                       );
                       */

  // Check if we have already has radiative background
  // if (ppath_what_background(ppath_step)) {
  //   termination_flag = ppath_what_background(ppath_step);
  //   g = 1;
  //   return;
  // }
  ARTS_USER_ERROR("ERROR")

  // Index in ppath_step of end point considered presently
  Index ip = 0;

  // Is point ip inside the cloudbox?
  inside_cloud = is_gp_inside_cloudbox(ppath_step.gp_p[ip],
                                       ppath_step.gp_lat[ip],
                                       ppath_step.gp_lon[ip],
                                       cloudbox_limits,
                                       0);

  // Determine radiative properties at point
  if (inside_cloud) {
    ARTS_USER_ERROR("NOT PORTED TO USE ATMFIELD BECAUSE I CANNOT DEAL WITH CLOUDS")
    /*
    cloudy_rt_vars_at_gp(ws,
                         ext_mat_mono,
                         abs_vec_mono,
                         pnd_vec,
                         temperature,
                         propmat_clearsky_agenda,
                         f_index,
                         f_grid,
                         ppath_step.gp_p[0],
                         ppath_step.gp_lat[0],
                         ppath_step.gp_lon[0],
                         p_grid[p_range],
                         t_field(p_range, lat_range, lon_range),
                         vmr_field(joker, p_range, lat_range, lon_range),
                         pnd_field,
                         scat_data,
                         cloudbox_limits,
                         Vector{ppath_step.los(0, joker)});
                         */
  } else {
    ARTS_USER_ERROR("NOT PORTED TO USE ATMFIELD BECAUSE I CANNOT DEAL WITH THE INTERPOLATION")
    /*
    clear_rt_vars_at_gp(ws,
                        ext_mat_mono,
                        abs_vec_mono,
                        temperature,
                        propmat_clearsky_agenda,
                        f_grid[f_index],
                        ppath_step.gp_p[0],
                        ppath_step.gp_lat[0],
                        ppath_step.gp_lon[0],
                        p_grid,
                        t_field,
                        vmr_field);
                        */
    pnd_vec = 0.0;
  }

  // Move the data to end point containers
  ext_matArray[1] = ext_mat_mono;
  abs_vecArray[1] = abs_vec_mono;
  tArray[1] = temperature;
  pnd_vecArray[1] = pnd_vec;

  //draw random number to determine end point
  Numeric r = rng.get(0.0, 1.0)();

  termination_flag = 0;

  while ((evol_op(0, 0) > r) && (!termination_flag)) {
    istep++;

    ARTS_USER_ERROR_IF (istep > 100000,
          "100000 path points have been reached. "
          "Is this an infinite loop?");

    evol_opArray[0] = evol_opArray[1];
    ext_matArray[0] = ext_matArray[1];
    abs_vecArray[0] = abs_vecArray[1];
    tArray[0] = tArray[1];
    pnd_vecArray[0] = pnd_vecArray[1];

    // The algorith to meet taustep_lim:
    // When first calculating a new ppath_step, it assumed that the present
    // ppath position holds the highest extinction. If the extinction at the
    // next position is higher, the criterion is checked and a new ppath_step
    // calculation is triggered if found necessary.
    // This should work in most cases, but is not 100% safe. Consider a case
    // with ppath_lmax = -1 and the extinction is zero at all grid box
    // corners except one. The two "test points" can then both get an
    // extinction of zero, while in fact is non-zero through the grid box and
    // the optical depth is underestimated. But this was handled equally bad
    // before taustep_limit was introduced (2016-10-10, PE)
    bool oktaustep = false;
    Index ppath_try = 1;
    const Index lmax_limit = 10;

    while (!oktaustep) {
      // Shall new ppath_step be calculated?
      if (ip == ppath_step.np - 1) {
        Numeric lmax = taustep_limit / ext_mat_mono(0, 0);
        if (ppath_lmax > 0) {
          lmax = min(ppath_lmax, lmax);
        }
        if (lmax < lmax_limit) {
          lmax = lmax_limit;
        }
        //cout << ppath_try << ", lmax = " << lmax << endl;
        //Print( ppath_step, 0 );

        ppath_step_agendaExecute(ws,
                                 ppath_step,
                                 lmax,
                                 ppath_lraytrace,
                                 Vector{f_grid[Range(f_index, 1)]},
                                 ppath_step_agenda);
        ip = 1;

        inside_cloud = is_gp_inside_cloudbox(ppath_step.gp_p[ip],
                                             ppath_step.gp_lat[ip],
                                             ppath_step.gp_lon[ip],
                                             cloudbox_limits,
                                             0);
      } else {
        ip++;
      }

      if (inside_cloud) {
    ARTS_USER_ERROR("NOT PORTED")
    /*
        cloudy_rt_vars_at_gp(ws,
                             ext_mat_mono,
                             abs_vec_mono,
                             pnd_vec,
                             temperature,
                             propmat_clearsky_agenda,
                             f_index,
                             f_grid,
                             ppath_step.gp_p[ip],
                             ppath_step.gp_lat[ip],
                             ppath_step.gp_lon[ip],
                             p_grid[p_range],
                             t_field(p_range, lat_range, lon_range),
                             vmr_field(joker, p_range, lat_range, lon_range),
                             pnd_field,
                             scat_data,
                             cloudbox_limits,
                             Vector{ppath_step.los(ip, joker)});
                             */
      } else {
    ARTS_USER_ERROR("NOT PORTED")
    /*
        clear_rt_vars_at_gp(ws,
                            ext_mat_mono,
                            abs_vec_mono,
                            temperature,
                            propmat_clearsky_agenda,
                            f_grid[f_index],
                            ppath_step.gp_p[ip],
                            ppath_step.gp_lat[ip],
                            ppath_step.gp_lon[ip],
                            p_grid,
                            t_field,
                            vmr_field);
                            */
        pnd_vec = 0.0;
      }

      dl = ppath_step.lstep[ip - 1];

      // Check if taustep_limit criterion is met. OK if:
      // 1. Ppath step already recalculated
      // 2. New ext_mat <= old one
      // 3. New ext_mat bigger, but tau of step still below limit
      if (ppath_try > 1 || ext_mat_mono(0, 0) <= ext_matArray[0](0, 0) ||
          (ext_mat_mono(0, 0) + ext_matArray[0](0, 0)) * dl / 2 <=
              taustep_limit) {
        oktaustep = true;
      } else {
        // We trigger a recalculation of ppath_step, from the previous
        // point
        ppath_step.np = ip;
        ip--;
        ppath_try = 2;
        // If a background found in first try this has to be reset:
        //ppath_set_background(ppath_step, 0);
        ARTS_USER_ERROR("ERROR")
      }
    }  // while !oktuastep

    ext_matArray[1] = ext_mat_mono;
    abs_vecArray[1] = abs_vec_mono;
    tArray[1] = temperature;
    pnd_vecArray[1] = pnd_vec;
    ext_mat = ext_matArray[1];
    ext_mat += ext_matArray[0];  // Factor 2 fixed by using dl/2
    //
    {
      Index extmat_case = 0;
      ext2trans(incT, extmat_case, ext_mat, dl / 2);
    }
    //
    mult(evol_op, evol_opArray[0], incT);
    evol_opArray[1] = evol_op;

    if (evol_op(0, 0) > r) {
      // Check whether hit ground or space.
      // path_step_agenda just detects surface intersections, and
      // if TOA is reached requires a special check.
      // But we are already ready if evol_op<=r
      if (ip == ppath_step.np - 1) {
        // if (ppath_what_background(ppath_step)) {
        //   termination_flag = 2;
        // }  //we have hit the surface
        ARTS_USER_ERROR("ERROR")
    ARTS_USER_ERROR("NOT PORTED")
    /*
        else if (fractional_gp(ppath_step.gp_p[ip]) >=
                 (Numeric)(p_grid.nelem() - 1) - 1e-3) {
          termination_flag = 1;
        } */ //we are at TOA
      }
    }
  }  // while

  if (termination_flag) {  //we must have reached the cloudbox boundary
    rte_pos = ppath_step.pos(ip, joker);
    rte_los = ppath_step.los(ip, joker);
    g = evol_op(0, 0);
  } else {
    //find position...and evol_op..and everything else required at the new
    //scattering/emission point
    // GH 2011-09-14:
    //   log(incT(0,0)) = log(exp(opt_depth_mat(0, 0))) = opt_depth_mat(0, 0)
    //   Avoid loss of precision, use opt_depth_mat directly
    //k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
    // PE 2013-05-17, Now the above comes directly from ext_mat:
    k = ext_mat(0, 0) / 2;  // Factor 2 as sum of end point values
    ds = log(evol_opArray[0](0, 0) / r) / k;
    g = k * r;
    Vector x(2, 0.0);
    //interpolate atmospheric variables required later
    ArrayOfGridPos gp(1);
    x[1] = dl;
    Vector itw(2);

    gridpos(gp, x, ConstVectorView{ds});
    ARTS_ASSERT(gp[0].idx == 0);
    interpweights(itw, gp[0]);
    interp(ext_mat_mono, itw, ext_matArray, gp[0]);
    ext_mat = ext_mat_mono;
    ext_mat += ext_matArray[gp[0].idx];
    //
    {
      Index extmat_case = 0;
      ext2trans(incT, extmat_case, ext_mat, ds / 2);
    }
    //
    mult(evol_op, evol_opArray[gp[0].idx], incT);
    interp(abs_vec_mono, itw, abs_vecArray, gp[0]);
    temperature = interp(itw, tArray, gp[0]);
    interp(pnd_vec, itw, pnd_vecArray, gp[0]);
    for (Index i = 0; i < 2; i++) {
      rte_pos[i] = interp(itw, ppath_step.pos(Range(ip - 1, 2), i), gp[0]);
      rte_los[i] = interp(itw, ppath_step.los(Range(ip - 1, 2), i), gp[0]);
    }
    rte_pos[2] = interp(itw, ppath_step.pos(Range(ip - 1, 2), 2), gp[0]);
  }

  ARTS_ASSERT(isfinite(g));

  // A dirty trick to avoid copying ppath_step
  const Index np = ip + 1;
  ppath_step.np = np;
}

void mcPathTraceRadar(const Workspace& ws,
                      MatrixView evol_op,
                      Vector& abs_vec_mono,
                      Numeric& temperature,
                      MatrixView ext_mat_mono,
                      RandomNumberGenerator<>& rng,
                      Vector& rte_pos,
                      Vector& rte_los,
                      Vector& pnd_vec,
                      Numeric& stot,
                      Numeric& ttot,
                      Ppath& ppath_step,
                      Index& termination_flag,
                      bool& inside_cloud,
                      const Agenda& ppath_step_agenda,
                      const Numeric& ppath_lmax,
                      const Numeric& ppath_lraytrace,
                      const Agenda& propmat_clearsky_agenda,
                      const bool& anyptype_nonTotRan,
                      const Index f_index,
                      const Vector& f_grid,
                      const Vector& Iprop,
                      const SurfaceField& surface_field,
                      const AtmField& atm_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Tensor4& pnd_field,
                      const ArrayOfArrayOfSingleScatteringData& scat_data) {
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2);
  ArrayOfVector pnd_vecArray(2);
  Matrix ext_mat(4, 4);
  Matrix incT(4, 4, 0.0);
  Vector tArray(2);
  Matrix T(4, 4);
  Numeric kI, kQ;
  Numeric ds, dt = -999, dl = -999;
  Index istep = 0;  // Counter for number of steps
  Matrix old_evol_op(4, 4);
  Vector local_rte_los(2);

  // Total path length starts at zero
  stot = 0.0;
  ttot = 0.0;

  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1] = evol_op;

  // range defining cloudbox
  Range p_range(cloudbox_limits[0],
                cloudbox_limits[1] - cloudbox_limits[0] + 1);
  Range lat_range(cloudbox_limits[2],
                  cloudbox_limits[3] - cloudbox_limits[2] + 1);
  Range lon_range(cloudbox_limits[4],
                  cloudbox_limits[5] - cloudbox_limits[4] + 1);

  //initialise Ppath with ppath_start_stepping
  /*  FIXME: OLD CODE MUST BE UPDATED
  ppath_start_stepping(ppath_step,
                       3,
                       p_grid,
                       lat_grid,
                       lon_grid,
                       z_field,
                       refellipsoid,
                       z_surface,
                       0,
                       cloudbox_limits,
                       false,
                       rte_pos,
                       rte_los,
                       );
                       */

  if (ppath_step.np == 0) {
    termination_flag = 1;
    return;
  }

  // Index in ppath_step of end point considered presently
  Index ip = 0;

  // Is point ip inside the cloudbox?
  inside_cloud = is_gp_inside_cloudbox(ppath_step.gp_p[ip],
                                       ppath_step.gp_lat[ip],
                                       ppath_step.gp_lon[ip],
                                       cloudbox_limits,
                                       0);

  // Determine radiative properties at point
  if (inside_cloud) {
    local_rte_los[0] = 180 - ppath_step.los(0, 0);
    local_rte_los[1] = ppath_step.los(0, 1) - 180;
      /*  FIXME: CLOUD STUFF
    cloudy_rt_vars_at_gp(ws,
                         ext_mat_mono,
                         abs_vec_mono,
                         pnd_vec,
                         temperature,
                         propmat_clearsky_agenda,
                         f_index,
                         f_grid,
                         ppath_step.gp_p[0],
                         ppath_step.gp_lat[0],
                         ppath_step.gp_lon[0],
                         p_grid[p_range],
                         t_field(p_range, lat_range, lon_range),
                         vmr_field(joker, p_range, lat_range, lon_range),
                         pnd_field,
                         scat_data,
                         cloudbox_limits,
                         local_rte_los);
                         */
  } else {
      /*  FIXME: CLOUD STUFF
    clear_rt_vars_at_gp(ws,
                        ext_mat_mono,
                        abs_vec_mono,
                        temperature,
                        propmat_clearsky_agenda,
                        f_grid[f_index],
                        ppath_step.gp_p[0],
                        ppath_step.gp_lat[0],
                        ppath_step.gp_lon[0],
                        p_grid,
                        t_field,
                        vmr_field);
                        */
    pnd_vec = 0.0;
  }

  // Move the data to end point containers
  ext_matArray[1] = ext_mat_mono;
  abs_vecArray[1] = abs_vec_mono;
  tArray[1] = temperature;
  pnd_vecArray[1] = pnd_vec;

  //draw random number to determine end point
  Numeric r = rng.get(0.0, 1.0)();

  termination_flag = 0;

  stot = ppath_step.end_lstep;
  ttot = ppath_step.end_lstep / SPEED_OF_LIGHT;

  Numeric evop0 = 1;
  while ((evop0 > r) && (!termination_flag)) {
    istep++;

    ARTS_USER_ERROR_IF (istep > 25000,
          "25000 path points have been reached. "
          "Is this an infinite loop?");

    evol_opArray[0] = evol_opArray[1];
    ext_matArray[0] = ext_matArray[1];
    abs_vecArray[0] = abs_vecArray[1];
    tArray[0] = tArray[1];
    pnd_vecArray[0] = pnd_vecArray[1];

    // Shall new ppath_step be calculated?
    if (ip >= ppath_step.np - 1) {
      ip = 1;
      ppath_step_agendaExecute(ws,
                               ppath_step,
                               ppath_lmax,
                               ppath_lraytrace,
                               Vector{f_grid[Range(f_index, 1)]},
                               ppath_step_agenda);

      if (ppath_step.np <= 1) {
        termination_flag = 1;
        break;
      }
      inside_cloud = is_gp_inside_cloudbox(ppath_step.gp_p[ip],
                                           ppath_step.gp_lat[ip],
                                           ppath_step.gp_lon[ip],
                                           cloudbox_limits,
                                           0);
    } else {
      ip++;
    }

    dl = ppath_step.lstep[ip - 1];
    dt = dl * 0.5 * (ppath_step.ngroup[ip - 1] + ppath_step.ngroup[ip]) /
         SPEED_OF_LIGHT;
    stot += dl;
    ttot += dt;
    if (inside_cloud) {
      local_rte_los[0] = 180 - ppath_step.los(ip, 0);
      local_rte_los[1] = ppath_step.los(ip, 1) - 180;
      /*  FIXME: CLOUD STUFF
      cloudy_rt_vars_at_gp(ws,
                           ext_mat_mono,
                           abs_vec_mono,
                           pnd_vec,
                           temperature,
                           propmat_clearsky_agenda,
                           f_index,
                           f_grid,
                           ppath_step.gp_p[ip],
                           ppath_step.gp_lat[ip],
                           ppath_step.gp_lon[ip],
                           p_grid[p_range],
                           t_field(p_range, lat_range, lon_range),
                           vmr_field(joker, p_range, lat_range, lon_range),
                           pnd_field,
                           scat_data,
                           cloudbox_limits,
                           local_rte_los);
                           */
    } else {
      /*  FIXME: CLOUD STUFF
      clear_rt_vars_at_gp(ws,
                          ext_mat_mono,
                          abs_vec_mono,
                          temperature,
                          propmat_clearsky_agenda,
                          f_grid[f_index],
                          ppath_step.gp_p[ip],
                          ppath_step.gp_lat[ip],
                          ppath_step.gp_lon[ip],
                          p_grid,
                          t_field,
                          vmr_field);
                          */
      pnd_vec = 0.0;
    }

    ext_matArray[1] = ext_mat_mono;
    abs_vecArray[1] = abs_vec_mono;
    tArray[1] = temperature;
    pnd_vecArray[1] = pnd_vec;
    ext_mat = ext_matArray[1];
    ext_mat += ext_matArray[0];  // Factor 2 fixed by using dl/2
    //
    {
      Index extmat_case = 0;
      ext2trans(incT, extmat_case, ext_mat, dl / 2);
    }
    //

    mult(evol_op, incT, evol_opArray[0]);
    evol_opArray[1] = evol_op;
    evop0 = evol_op(0, 0);

    // Handle cross-talk for ptype==30
    if (anyptype_nonTotRan) {
      const Numeric Q1 = evol_op(0, 1) * Iprop[1] / Iprop[0];
      evop0 += Q1;
    }
    if (evop0 > r) {
      // Check whether hit ground or space.
      // path_step_agenda just detects surface intersections, and
      // if TOA is reached requires a special check.
      // But we are already ready if evol_op<=r
      if (ip >= ppath_step.np - 1) {
        // if (ppath_what_background(ppath_step)) {
        //   termination_flag = 2;
        // }  //we have hit the surface
        // else if (fractional_gp(ppath_step.gp_p[ip]) >=
        //          (Numeric)(atm_field.regularized_shape()[0] - 1) - 1e-3) {
        //   termination_flag = 1;
        // }  //we are at TOA
        ARTS_USER_ERROR("ERROR")
      }
    }
  }  // while

  if (termination_flag != 0) {  //we must have reached the cloudbox boundary
    if (ip < ppath_step.np) {
      // This is not used if termination flag set,
      // so not an issue of ip is too large.
      // Need to think about this.
      rte_pos = ppath_step.pos(ip, joker);
      rte_los = ppath_step.los(ip, joker);
    }
  } else {
    //find position...and evol_op..and everything else required at the new
    //scattering/emission point
    const Numeric tol = 0.1;  // Tolerance of 10 cm
    stot -= dl;  // Take out last step because we are between stepping points
    ttot -= dt;  // Take out last step because we are between stepping points
    kI = ext_mat(0, 0) / 2;  // Factor 2 as sum of end point values
    kQ = ext_mat(0, 1) / 2;  // Factor 2 as sum of end point values
    if (anyptype_nonTotRan) {
      const Numeric I1 = evol_opArray[0](0, 0);
      const Numeric Q1 = evol_opArray[0](0, 1) * Iprop[1] / Iprop[0];

      // Need to use root finding to solve for ds
      brent_zero(
          ds, (Numeric)0.0, ppath_step.lstep[ip - 1], tol, r, I1, Q1, kI, kQ);
    } else {
      // Simple inversion when no cross-talk between I and Q
      ds = log(evol_opArray[0](0, 0) / r) / kI;
    }
    stot += ds;
    ttot += ds * dt / dl;
    Vector x(2, 0.0);

    //interpolate atmospheric variables required later
    ArrayOfGridPos gp(1);
    x[1] = dl;
    Vector itw(2);
    gridpos(gp, x, ConstVectorView{ds});
    ARTS_ASSERT(gp[0].idx == 0);
    interpweights(itw, gp[0]);
    interp(ext_mat_mono, itw, ext_matArray, gp[0]);
    ext_mat = ext_mat_mono;
    ext_mat += ext_matArray[gp[0].idx];
    //
    {
      Index extmat_case = 0;
      ext2trans(incT, extmat_case, ext_mat, ds / 2);
    }
    //
    mult(evol_op, incT, evol_opArray[gp[0].idx]);
    interp(abs_vec_mono, itw, abs_vecArray, gp[0]);
    temperature = interp(itw, tArray, gp[0]);
    interp(pnd_vec, itw, pnd_vecArray, gp[0]);
    for (Index i = 0; i < 2; i++) {
      rte_pos[i] = interp(itw, ppath_step.pos(Range(ip - 1, 2), i), gp[0]);
      rte_los[i] = interp(itw, ppath_step.los(Range(ip - 1, 2), i), gp[0]);
    }
    rte_pos[2] = interp(itw, ppath_step.pos(Range(ip - 1, 2), 2), gp[0]);
  }

  // A dirty trick to avoid copying ppath_step
  const Index np = ip + 1;
  ppath_step.np = np;
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
  Index np = pnd_vec.nelem();
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
