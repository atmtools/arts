/* Copyright (C) 2003-2012 Claudia Emde <claudia.emde@dlr.de>

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

/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   optproperties.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu Mar  6 11:29:59 2003

  \brief  This file contains definitions and functions related to the
          optical properties of particles.


 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "optproperties.h"
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "check_input.h"
#include "constants.h"
#include "interpolation.h"
#include "interpolation_lagrange.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "xml_io.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

#define F11 pha_mat_int[0]
#define F12 pha_mat_int[1]
#define F22 pha_mat_int[2]
#define F33 pha_mat_int[3]
#define F34 pha_mat_int[4]
#define F44 pha_mat_int[5]

//! one-line descript
/*!
  Descript/Doc

  \param[out] name  desc.
  \param[in]  name  desc.

  \author Jana Mendrok
  \date   2018-01-15
*/
/*
void methodname(//Output
                type& name,
                //Input
                const type& name)
{
}
*/

//! Scattering species bulk extinction and absorption.
/*!
  Derives bulk properties from per-scat-species bulk properties.

  Ptype is defined by the most complex ptype of the individual scattering
  elements.

  \param[out] ext_mat     Bulk extinction matrix (over freq, temp/location,
                            propagation direction).
  \param[out] abs_vec     Bulk absorption vector (over freq, temp/location,
                            propagation direction).
  \param[out] ptype       Bulk ptype.
  \param[in]  ext_mat_ss  Bulk extinction matrix per scattering species (over
                            freq, temp/location, propagation direction).
  \param[in]  abs_vec_ss  Bulk absorption vector per scattering species (over
                            freq, temp/location, propagation direction).
  \param[in]  ptypes_ss   Scattering species ptypes.

  \author Jana Mendrok
  \date   2018-01-16
*/
void opt_prop_Bulk(    //Output
    Tensor5& ext_mat,  // (nf,nT,ndir,nst,nst)
    Tensor4& abs_vec,  // (nf,nT,ndir,nst)
    Index& ptype,
    //Input
    const ArrayOfTensor5& ext_mat_ss,  // [nss](nf,nT,ndir,nst,nst)
    const ArrayOfTensor4& abs_vec_ss,  // [nss](nf,nT,ndir,nst)
    const ArrayOfIndex& ptypes_ss) {
  ARTS_ASSERT(ext_mat_ss.nelem() == abs_vec_ss.nelem());

  ext_mat = ext_mat_ss[0];
  abs_vec = abs_vec_ss[0];

  for (Index i_ss = 1; i_ss < ext_mat_ss.nelem(); i_ss++) {
    ext_mat += ext_mat_ss[i_ss];
    abs_vec += abs_vec_ss[i_ss];
  }
  ptype = max(ptypes_ss);
}

//! Scattering species bulk extinction and absorption.
/*!
  Derives bulk properties separately per scattering species from (input)
  per-scattering-elements extinction and absorption given over frequency,
  temperature and propagation direction.

  Temperature dimension in per-scattering-element properties is assumed to
  correspond to column dimension of pnds (leading dimension, i.e. row,
  corresponds to scattering elements).

  Ptype of each scattering species is defined by the most complex ptype of the
  individual scattering elements within this species.

  \param[out] ext_mat     Bulk extinction matrix per scattering species (over
                            freq, temp/location, propagation direction).
  \param[out] abs_vec     Bulk absorption vector per scattering species (over
                            freq, temp/location, propagation direction).
  \param[out] ptype       Scattering species ptypes.
  \param[in]  ext_mat_se  Extinction matrix per scattering element.
  \param[in]  abs_vec_se  Absorption vector per scattering element.
  \param[in]  ptypes_se   Scattering element types.
  \param[in]  pnds        Particle number density vectors (multiple locations).

  \author Jana Mendrok
  \date   2018-01-16
*/
void opt_prop_ScatSpecBulk(   //Output
    ArrayOfTensor5& ext_mat,  // [nss](nf,nT,ndir,nst,nst)
    ArrayOfTensor4& abs_vec,  // [nss](nf,nT,ndir,nst)
    ArrayOfIndex& ptype,
    //Input
    const ArrayOfArrayOfTensor5& ext_mat_se,  // [nss][nse](nf,nT,ndir,nst,nst)
    const ArrayOfArrayOfTensor4& abs_vec_se,  // [nss][nse](nf,nT,ndir,nst)
    const ArrayOfArrayOfIndex& ptypes_se,
    ConstMatrixView pnds,
    ConstMatrixView t_ok) {
  ARTS_ASSERT(t_ok.nrows() == pnds.nrows());
  ARTS_ASSERT(t_ok.ncols() == pnds.ncols());
  ARTS_ASSERT(TotalNumberOfElements(ext_mat_se) == pnds.nrows());
  ARTS_ASSERT(TotalNumberOfElements(abs_vec_se) == pnds.nrows());
  ARTS_ASSERT(ext_mat_se.nelem() == abs_vec_se.nelem());

  const Index nT = pnds.ncols();
  const Index nf = abs_vec_se[0][0].nbooks();
  const Index nDir = abs_vec_se[0][0].nrows();
  const Index stokes_dim = abs_vec_se[0][0].ncols();

  const Index nss = ext_mat_se.nelem();
  ext_mat.resize(nss);
  abs_vec.resize(nss);
  ptype.resize(nss);
  Tensor4 ext_tmp;
  Tensor3 abs_tmp;

  Index i_se_flat = 0;

  for (Index i_ss = 0; i_ss < nss; i_ss++) {
    ARTS_ASSERT(ext_mat_se[i_ss].nelem() == abs_vec_se[i_ss].nelem());
    ARTS_ASSERT(nT == ext_mat_se[i_ss][0].nbooks());
    ARTS_ASSERT(nT == abs_vec_se[i_ss][0].npages());

    ext_mat[i_ss].resize(nf, nT, nDir, stokes_dim, stokes_dim);
    ext_mat[i_ss] = 0.;
    abs_vec[i_ss].resize(nf, nT, nDir, stokes_dim);
    abs_vec[i_ss] = 0.;

    for (Index i_se = 0; i_se < ext_mat_se[i_ss].nelem(); i_se++) {
      ARTS_ASSERT(nT == ext_mat_se[i_ss][i_se].nbooks());
      ARTS_ASSERT(nT == abs_vec_se[i_ss][i_se].npages());

      for (Index Tind = 0; Tind < nT; Tind++) {
        if (pnds(i_se_flat, Tind) != 0.) {
          if (t_ok(i_se_flat, Tind) > 0.) {
            ext_tmp = ext_mat_se[i_ss][i_se](joker, Tind, joker, joker, joker);
            ext_tmp *= pnds(i_se_flat, Tind);
            ext_mat[i_ss](joker, Tind, joker, joker, joker) += ext_tmp;

            abs_tmp = abs_vec_se[i_ss][i_se](joker, Tind, joker, joker);
            abs_tmp *= pnds(i_se_flat, Tind);
            abs_vec[i_ss](joker, Tind, joker, joker) += abs_tmp;
          } else {
            ARTS_USER_ERROR (
              "Interpolation error for (flat-array) scattering element #",
              i_se_flat, "\n"
              "at location/temperature point #", Tind, "\n")
          }
        }
      }
      i_se_flat++;
    }
    ptype[i_ss] = max(ptypes_se[i_ss]);
  }
}

//! Extinction and absorption from all scattering elements.
/*!
  Derives temperature and direction interpolated extinction matrices and
  absorption vectors for all scattering elements present in scat_data.

  ATTENTION:
  If scat_data has only one freq point, f_index=-1 (i.e. all) extracts only this
  one freq point. To duplicate that as needed if f_grid has more freqs is TASK
  of the CALLING METHOD!

  Loops over opt_prop_1ScatElem and packs its output into all-scat-elements
  containers.

  \param[out] ext_mat    Extinction matrix (over scat elements, freq, temp,
                           propagation direction).
  \param[out] abs_vec    Absorption vector (over scat elements, freq, temp,
                           propagation direction).
  \param[out] ptypes     Scattering element types.
  \param[out] t_ok       Flag whether T-interpol valid (over scat elements, temp).
  \param[in]  scat_data  as the WSV.
  \param[in]  stokes_dim as the WSV.
  \param[in]  T_array    Temperatures to extract ext/abs for.
  \param[in]  dir_array  Propagation directions to extract ext/abs for (as
                           pairs of zenith and azimuth angle per direction).
  \param[in]  f_index    Index of frequency to extract. -1 extracts data for all
                           freqs available in ssd.
  \param[in]  t_interp_order  Temperature interpolation order.

  \author Jana Mendrok
  \date   2018-01-16
*/
void opt_prop_NScatElems(            //Output
    ArrayOfArrayOfTensor5& ext_mat,  // [nss][nse](nf,nT,ndir,nst,nst)
    ArrayOfArrayOfTensor4& abs_vec,  // [nss][nse](nf,nT,ndir,nst)
    ArrayOfArrayOfIndex& ptypes,
    Matrix& t_ok,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& stokes_dim,
    const Vector& T_array,
    const Matrix& dir_array,
    const Index& f_index,
    const Index& t_interp_order) {
  Index f_start, nf;
  if (f_index < 0) {
    nf = scat_data[0][0].ext_mat_data.nshelves();
    f_start = 0;
    //f_end = f_start+nf;
  } else {
    nf = 1;
    if (scat_data[0][0].ext_mat_data.nshelves() == 1)
      f_start = 0;
    else
      f_start = f_index;
    //f_end = f_start+nf;
  }

  const Index nT = T_array.nelem();
  const Index nDir = dir_array.nrows();

  const Index nss = scat_data.nelem();
  ext_mat.resize(nss);
  abs_vec.resize(nss);
  ptypes.resize(nss);

  const Index Nse_all = TotalNumberOfElements(scat_data);
  t_ok.resize(Nse_all, nT);
  Index i_se_flat = 0;

  for (Index i_ss = 0; i_ss < nss; i_ss++) {
    Index nse = scat_data[i_ss].nelem();
    ext_mat[i_ss].resize(nse);
    abs_vec[i_ss].resize(nse);
    ptypes[i_ss].resize(nse);

    for (Index i_se = 0; i_se < nse; i_se++) {
      ext_mat[i_ss][i_se].resize(nf, nT, nDir, stokes_dim, stokes_dim);
      abs_vec[i_ss][i_se].resize(nf, nT, nDir, stokes_dim);

      opt_prop_1ScatElem(ext_mat[i_ss][i_se],
                         abs_vec[i_ss][i_se],
                         ptypes[i_ss][i_se],
                         t_ok(i_se_flat, joker),
                         scat_data[i_ss][i_se],
                         T_array,
                         dir_array,
                         f_start,
                         t_interp_order);
      i_se_flat++;
    }
  }
  ARTS_ASSERT(i_se_flat == Nse_all);
}

//! Determine T-interpol parameters for a specific scattering element.
/*!
  Determine T-interpol order as well as interpol positions and weights (they
  are the same for all directions (and freqs), ie it is sufficient to
  calculate them once).

  \param[out] t_ok       Flag whether T-interpol valid (length of T_array).
  \param[out] this_T_interp_order  Temperature interpolation order adjusted
                           according to data availability.
  \param[out] T_lag      GridPos data for temperature interpolation.
  \param[out] T_itw      Interpolation weights for temperature interpolation.
  \param[in]  T_grid     Single scattering data temperature grid.
  \param[in]  T_array    Temperatures to extract single scattering data for.
  \param[in]  t_interp_order  Requested temperature interpolation order.

  \author Jana Mendrok
  \date   2018-05-01
*/
ArrayOfLagrangeInterpolation ssd_tinterp_parameters(  //Output
    VectorView t_ok,
    Index& this_T_interp_order,
    //Input
    ConstVectorView T_grid,
    const Vector& T_array,
    const Index& t_interp_order) {
  const Index nTin = T_grid.nelem();
  const Index nTout = T_array.nelem();

  this_T_interp_order = -1;


  if (nTin > 1) {
    this_T_interp_order = min(t_interp_order, nTin - 1);

    // we need to check T-grid exceedance. and catch these cases (because T
    // is assumed to correspond to a location and T-exceedance is ok when pnd=0
    // there. however, here we do not know pnd.) and report them to have the
    // calling method dealing with it.

    // we set the extrapolfax explicitly and use it here as well as in
    // gridpos_poly below.
    const Numeric extrapolfac = 0.5;
    const Numeric lowlim = T_grid[0] - extrapolfac * (T_grid[1] - T_grid[0]);
    const Numeric uplim =
        T_grid[nTin - 1] + extrapolfac * (T_grid[nTin - 1] - T_grid[nTin - 2]);

    bool any_T_exceed = false;
    for (Index Tind = 0; Tind < nTout; Tind++) {
      if (T_array[Tind] < lowlim || T_array[Tind] > uplim) {
        t_ok[Tind] = -1.;
        any_T_exceed = true;
      } else
        t_ok[Tind] = 1.;
    }

    if (any_T_exceed) {
      // Reserve output
      ArrayOfLagrangeInterpolation T_lag;
      T_lag.reserve(nTout);

      bool grid_unchecked = true;

      for (Index iT = 0; iT < nTout; iT++) {
        if (t_ok[iT] < 0) {
          T_lag.emplace_back();
        } else {
          if (grid_unchecked) {
            chk_interpolation_grids(
                "Temperature interpolation in pha_mat_1ScatElem",
                T_grid,
                T_array[Range(iT, 1)],
                this_T_interp_order);
            grid_unchecked = false;
          }
          T_lag.emplace_back(0, T_array[iT], T_grid, this_T_interp_order, false);
        }
      }
      return T_lag;
    } else {
      return Interpolation::LagrangeVector(T_array, T_grid, this_T_interp_order, extrapolfac);
    }
  } else {
    t_ok = 1.;
    return {};
  }
}

//! Preparing extinction and absorption from one scattering element.
/*!
  Extracts and prepares extinction matrix and absorption vector data for one
  scattering element for one or all frequencies from the single scattering data.
  Includes interpolation in temperature and to propagation direction. Handles
  multiple output temperatures and propagation directions at a time. Temperature
  interpolation order can be chosen.

  \param[out] ext_mat    1-scattering element extinction matrix (over freq,
                           temp, propagation direction).
  \param[out] abs_vec    1-scattering element absorption vector (over freq,
                           temp, propagation direction).
  \param[out] ptype      Type of scattering element.
  \param[out] t_ok       Flag whether T-interpol valid (length of T_array).
  \param[in]  ssd        Single scattering data of one scattering element.
  \param[in]  T_array    Temperatures to extract ext/abs for.
  \param[in]  dir_array  Propagation directions to extract ext/abs for (as
                           pairs of zenith and azimuth angle per direction).
  \param[in]  f_start    Start index of frequency/ies to extract.
  \param[in]  t_interp_order  Temperature interpolation order.

  \author Jana Mendrok
  \date   2018-01-15
*/
void opt_prop_1ScatElem(  //Output
    Tensor5View ext_mat,  // nf, nT, ndir, nst, nst
    Tensor4View abs_vec,  // nf, nT, ndir, nst
    Index& ptype,
    VectorView t_ok,
    //Input
    const SingleScatteringData& ssd,
    const Vector& T_array,
    const Matrix& dir_array,
    const Index& f_start,
    const Index& t_interp_order) {
  // FIXME: this is prob best done in scat_data_checkedCalc (or
  // cloudbox_checkedCalc) to have it done once and for all. Here at max ARTS_ASSERT.

  // At very first check validity of the scatt elements ptype (so far we only
  // handle PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND).
  /*
  if( ssd.ptype != PTYPE_TOTAL_RND and ssd.ptype != PTYPE_AZIMUTH_RND )
  {
    ostringstream os;
    os << "Only ptypes " << PTYPE_TOTAL_RND << " and " << PTYPE_AZIMUTH_RND
       << " can be handled.\n"
       << "Encountered scattering element with ptype " << ssd.ptype
       << ", though.";
    throw runtime_error( os.str() );
  }
  */

  ARTS_ASSERT(ssd.ptype == PTYPE_TOTAL_RND or ssd.ptype == PTYPE_AZIMUTH_RND);

  const Index nf = ext_mat.nshelves();
  ARTS_ASSERT(abs_vec.nbooks() == nf);
  if (nf > 1) {
    ARTS_ASSERT(nf == ssd.f_grid.nelem());
  }

  const Index nTout = T_array.nelem();
  ARTS_ASSERT(ext_mat.nbooks() == nTout);
  ARTS_ASSERT(abs_vec.npages() == nTout);
  ARTS_ASSERT(t_ok.nelem() == nTout);

  const Index nDir = dir_array.nrows();
  ARTS_ASSERT(ext_mat.npages() == nDir);
  ARTS_ASSERT(abs_vec.nrows() == nDir);

  const Index stokes_dim = abs_vec.ncols();
  ARTS_ASSERT(ext_mat.nrows() == stokes_dim);
  ARTS_ASSERT(ext_mat.ncols() == stokes_dim);

  ptype = ssd.ptype;

  // Determine T-interpol order as well as interpol positions and weights (they
  // are the same for all directions (and freqs), ie it is sufficient to
  // calculate them once).
  const Index nTin = ssd.T_grid.nelem();
  Index this_T_interp_order;
  const auto T_lag =  ssd_tinterp_parameters(t_ok,
                                             this_T_interp_order,
                                             ssd.T_grid,
                                             T_array,
                                             t_interp_order);
  const auto T_itw_lag = interpweights(T_lag);

  // Now loop over requested directions (and apply simultaneously for all freqs):
  // 1) extract/interpolate direction (but not for tot.random)
  // 2) apply T-interpol
  if (ptype == PTYPE_TOTAL_RND)
    if (this_T_interp_order < 0)  // just extract (and unpack) ext and abs data
                                  // for the fs, the Tin, and the dirin and sort
                                  // (copy) into the output fs, Ts, and dirs.
    {
      Tensor3 ext_mat_tmp(nf, stokes_dim, stokes_dim);
      Matrix abs_vec_tmp(nf, stokes_dim);
      for (Index find = 0; find < nf; find++) {
        ext_mat_SSD2Stokes(ext_mat_tmp(find, joker, joker),
                           ssd.ext_mat_data(find + f_start, 0, 0, 0, joker),
                           stokes_dim,
                           ptype);
        abs_vec_SSD2Stokes(abs_vec_tmp(find, joker),
                           ssd.abs_vec_data(find + f_start, 0, 0, 0, joker),
                           stokes_dim,
                           ptype);
      }
      for (Index Tind = 0; Tind < nTout; Tind++)
        for (Index dind = 0; dind < nDir; dind++) {
          ext_mat(joker, Tind, dind, joker, joker) = ext_mat_tmp;
          abs_vec(joker, Tind, dind, joker) = abs_vec_tmp;
        }
    } else  // T-interpolation required (but not dir). To be done on the compact
            // ssd format.
    {
      Tensor4 ext_mat_tmp(nf, nTout, stokes_dim, stokes_dim);
      Tensor3 abs_vec_tmp(nf, nTout, stokes_dim);
      Matrix ext_mat_tmp_ssd(nTout, ssd.ext_mat_data.ncols());
      Matrix abs_vec_tmp_ssd(nTout, ssd.abs_vec_data.ncols());
      for (Index find = 0; find < nf; find++) {
        for (Index nst = 0; nst < ext_mat_tmp_ssd.ncols(); nst++) {
          reinterp(ext_mat_tmp_ssd(joker, nst),
                   ssd.ext_mat_data(find + f_start, joker, 0, 0, nst),
                   T_itw_lag, T_lag);
        }
        for (Index Tind = 0; Tind < nTout; Tind++)
          if (t_ok[Tind] > 0.)
            ext_mat_SSD2Stokes(ext_mat_tmp(find, Tind, joker, joker),
                               ext_mat_tmp_ssd(Tind, joker),
                               stokes_dim,
                               ptype);
          else
            ext_mat_tmp(find, Tind, joker, joker) = 0.;

        for (Index nst = 0; nst < abs_vec_tmp_ssd.ncols(); nst++) {
          reinterp(abs_vec_tmp_ssd(joker, nst),
                   ssd.abs_vec_data(find + f_start, joker, 0, 0, nst),
                   T_itw_lag,
                   T_lag);
        }
        for (Index Tind = 0; Tind < nTout; Tind++)
          if (t_ok[Tind] > 0.)
            abs_vec_SSD2Stokes(abs_vec_tmp(find, Tind, joker),
                               abs_vec_tmp_ssd(Tind, joker),
                               stokes_dim,
                               ptype);
          else
            abs_vec_tmp(find, Tind, joker) = 0.;
      }

      for (Index dind = 0; dind < nDir; dind++) {
        ext_mat(joker, joker, dind, joker, joker) = ext_mat_tmp;
        abs_vec(joker, joker, dind, joker) = abs_vec_tmp;
      }
    }

  else  // dir-interpolation for non-tot-random particles
  {
    // as we don't allow other than az.random here, ext and abs will only vary
    // with za, not with aa. Hence, we could be smart here and calc data only
    // for unique za (and copy the rest). however, this smartness might cost as
    // well. so for now, we leave that and blindly loop over the given direction
    // array, and leave smartness in setting up directional array to the calling
    // methods.

    // derive the direction interpolation weights.
    ArrayOfGridPos dir_gp(nDir);
    gridpos(dir_gp, ssd.za_grid, dir_array(joker, 0));
    Matrix dir_itw(nDir, 2);  // only interpolating in za, ie 1D linear interpol
    interpweights(dir_itw, dir_gp);

    const Index next = ssd.ext_mat_data.ncols();
    const Index nabs = ssd.abs_vec_data.ncols();

    if (this_T_interp_order < 0)  // T only needs to be extracted.
    {
      Matrix ext_mat_tmp_ssd(nDir, next);
      Matrix abs_vec_tmp_ssd(nDir, nabs);
      Tensor4 ext_mat_tmp(nf, nDir, stokes_dim, stokes_dim);
      Tensor3 abs_vec_tmp(nf, nDir, stokes_dim);
      for (Index find = 0; find < nf; find++) {
        for (Index nst = 0; nst < next; nst++)
          interp(ext_mat_tmp_ssd(joker, nst),
                 dir_itw,
                 ssd.ext_mat_data(find + f_start, 0, joker, 0, nst),
                 dir_gp);
        for (Index Dind = 0; Dind < nDir; Dind++)
          ext_mat_SSD2Stokes(ext_mat_tmp(find, Dind, joker, joker),
                             ext_mat_tmp_ssd(Dind, joker),
                             stokes_dim,
                             ptype);

        for (Index nst = 0; nst < nabs; nst++)
          interp(abs_vec_tmp_ssd(joker, nst),
                 dir_itw,
                 ssd.abs_vec_data(find + f_start, 0, joker, 0, nst),
                 dir_gp);
        for (Index Dind = 0; Dind < nDir; Dind++)
          abs_vec_SSD2Stokes(abs_vec_tmp(find, Dind, joker),
                             abs_vec_tmp_ssd(Dind, joker),
                             stokes_dim,
                             ptype);
      }

      for (Index Tind = 0; Tind < nTout; Tind++) {
        ext_mat(joker, Tind, joker, joker, joker) = ext_mat_tmp;
        abs_vec(joker, Tind, joker, joker) = abs_vec_tmp;
      }
    } else  // T- and dir-interpolation required. To be done on the compact ssd
            // format.
    {
      Tensor3 ext_mat_tmp_ssd(nTin, nDir, next);
      Tensor3 abs_vec_tmp_ssd(nTin, nDir, nabs);
      Matrix ext_mat_tmp(nTout, next);
      Matrix abs_vec_tmp(nTout, nabs);

      for (Index find = 0; find < nf; find++) {
        for (Index Tind = 0; Tind < nTin; Tind++) {
          for (Index nst = 0; nst < next; nst++)
            interp(ext_mat_tmp_ssd(Tind, joker, nst),
                   dir_itw,
                   ssd.ext_mat_data(find + f_start, Tind, joker, 0, nst),
                   dir_gp);
          for (Index nst = 0; nst < nabs; nst++)
            interp(abs_vec_tmp_ssd(Tind, joker, nst),
                   dir_itw,
                   ssd.abs_vec_data(find + f_start, Tind, joker, 0, nst),
                   dir_gp);
        }

        for (Index Dind = 0; Dind < nDir; Dind++) {
          for (Index nst = 0; nst < next; nst++) {
            reinterp(ext_mat_tmp(joker, nst),
                     ext_mat_tmp_ssd(joker, Dind, nst),
                     T_itw_lag,
                     T_lag);
          }
          for (Index Tind = 0; Tind < nTout; Tind++)
            ext_mat_SSD2Stokes(ext_mat(find, Tind, Dind, joker, joker),
                               ext_mat_tmp(Tind, joker),
                               stokes_dim,
                               ptype);

          for (Index nst = 0; nst < nabs; nst++) {
            reinterp(abs_vec_tmp(joker, nst),
                     abs_vec_tmp_ssd(joker, Dind, nst),
                     T_itw_lag,
                     T_lag);
          }
          for (Index Tind = 0; Tind < nTout; Tind++)
            abs_vec_SSD2Stokes(abs_vec(find, Tind, Dind, joker),
                               abs_vec_tmp(Tind, joker),
                               stokes_dim,
                               ptype);
        }
      }
    }
  }
}

//! Extinction matrix scat_data to stokes format conversion.
/*!
  Converts extinction matrix from scat_data ptype-dependent compact format to
  Stokes-notation matrix.

  \param[out] ext_mat_stokes  extmat in stokes notation (stokes_dim,stokes_dim).
  \param[in]  ext_mat_ssd     extmat at 1f, 1T, 1dir in scat_data compact
                                (vector) format.
  \param[in]  stokes_dim      as the WSM.
  \param[in]  ptype           type of scattering element.

  \author Jana Mendrok
  \date   2018-01-15
*/
void ext_mat_SSD2Stokes(  //Output
    MatrixView ext_mat_stokes,
    //Input
    ConstVectorView ext_mat_ssd,
    const Index& stokes_dim,
    const Index& ptype) {
  // for now, no handling of PTYPE_GENERAL. should be ensured somewhere in the
  // calling methods, though.
  ARTS_ASSERT(ptype <= PTYPE_AZIMUTH_RND);

  ext_mat_stokes = 0.;

  for (Index ist = 0; ist < stokes_dim; ist++) {
    ext_mat_stokes(ist, ist) = ext_mat_ssd[0];
  }

  if (ptype > PTYPE_TOTAL_RND) {
    switch (stokes_dim) {
      case 4: {
        ext_mat_stokes(2, 3) = ext_mat_ssd[2];
        ext_mat_stokes(3, 2) = -ext_mat_ssd[2];
      } /* FALLTHROUGH */
      case 3: {
        // nothing to be done here. but we need this for executing the below
        // also in case of stokes_dim==3.
      }
      case 2: {
        ext_mat_stokes(0, 1) = ext_mat_ssd[1];
        ext_mat_stokes(1, 0) = ext_mat_ssd[1];
      }
    }
  }
}

//! Absorption vector scat_data to stokes format conversion.
/*!
  Converts absorption vector from scat_data ptype-dependent compact format to
  Stokes-notation matrix.

  \param[out] abs_vec_stokes  absvec in stokes notation (stokes_dim).
  \param[in]  abs_vec_ssd     absvec at 1f, 1T, 1dir in scat_data compact
                                (vector) format.
  \param[in]  stokes_dim      as the WSM.
  \param[in]  ptype           Type of scattering element.

  \author Jana Mendrok
  \date   2018-01-15
*/
void abs_vec_SSD2Stokes(  //Output
    VectorView abs_vec_stokes,
    //Input
    ConstVectorView abs_vec_ssd,
    const Index& stokes_dim,
    const Index& ptype) {
  // for now, no handling of PTYPE_GENERAL. should be ensured somewhere in the
  // calling methods, though.
  ARTS_ASSERT(ptype <= PTYPE_AZIMUTH_RND);

  abs_vec_stokes = 0.;

  abs_vec_stokes[0] = abs_vec_ssd[0];

  if (ptype > PTYPE_TOTAL_RND and stokes_dim > 1) {
    abs_vec_stokes[1] = abs_vec_ssd[1];
  }
}

//! Scattering species bulk phase matrix.
/*!
  Derives bulk properties from per-scat-species bulk properties.

  Ptype is defined by the most complex ptype of the individual scattering
  elements.

  \param[out] pha_mat     Bulk phase matrix (over freq, temp/location,
                            propagation direction, incident direction).
  \param[out] ptype       Bulk ptype.
  \param[in]  pha_mat_ss  Bulk phase matrix per scattering species (over freq,
                            temp/location, propagation direction, incident
                            direction).
  \param[in]  ptypes_ss   Scattering species ptypes.

  \author Jana Mendrok
  \date   2018-03-24
*/
void pha_mat_Bulk(     //Output
    Tensor6& pha_mat,  // (nf,nT,npdir,nidir,nst,nst)
    Index& ptype,
    //Input
    const ArrayOfTensor6& pha_mat_ss,  // [nss](nf,nT,npdir,nidir,nst,nst)
    const ArrayOfIndex& ptypes_ss) {
  pha_mat = pha_mat_ss[0];

  for (Index i_ss = 1; i_ss < pha_mat_ss.nelem(); i_ss++)
    pha_mat += pha_mat_ss[i_ss];

  ptype = max(ptypes_ss);
}

//! Scattering species bulk phase matrices.
/*!
  Derives bulk properties separately per scattering species from (input)
  per-scattering-elements phase matrices given over frequency, temperature and
  propagation direction.

  Temperature dimension in per-scattering-element properties is assumed to
  correspond to column dimension of pnds (leading dimension, i.e. row,
  corresponds to scattering elements).

  Ptype of each scattering species is defined by the most complex ptype of the
  individual scattering elements within this species.

  \param[out] pha_mat     Bulk phase matrix per scattering species (over freq,
                            temp/location, propagation direction, incident
                            direction).
  \param[out] ptype       Scattering species ptypes.
  \param[in]  pha_mat_se  Phase matrix per scattering element.
  \param[in]  ptypes_se   Scattering element types.
  \param[in]  pnds        Particle number density vectors (multiple locations).

  \author Jana Mendrok
  \date   2018-03-24
*/
void pha_mat_ScatSpecBulk(    //Output
    ArrayOfTensor6& pha_mat,  // [nss](nf,nT,npdir,nidir,nst,nst)
    ArrayOfIndex& ptype,
    //Input
    const ArrayOfArrayOfTensor6&
        pha_mat_se,  // [nss][nse](nf,nT,npdir,nidir,nst,nst)
    const ArrayOfArrayOfIndex& ptypes_se,
    ConstMatrixView pnds,
    ConstMatrixView t_ok) {
  ARTS_ASSERT(t_ok.nrows() == pnds.nrows());
  ARTS_ASSERT(t_ok.ncols() == pnds.ncols());
  ARTS_ASSERT(TotalNumberOfElements(pha_mat_se) == pnds.nrows());

  const Index nT = pnds.ncols();
  const Index nf = pha_mat_se[0][0].nvitrines();
  const Index npDir = pha_mat_se[0][0].nbooks();
  const Index niDir = pha_mat_se[0][0].npages();
  const Index stokes_dim = pha_mat_se[0][0].ncols();

  const Index nss = pha_mat_se.nelem();
  pha_mat.resize(nss);
  ptype.resize(nss);
  Tensor5 pha_tmp;

  Index i_se_flat = 0;

  for (Index i_ss = 0; i_ss < nss; i_ss++) {
    ARTS_ASSERT(nT == pha_mat_se[i_ss][0].nshelves());

    pha_mat[i_ss].resize(nf, nT, npDir, niDir, stokes_dim, stokes_dim);
    pha_mat[i_ss] = 0.;

    for (Index i_se = 0; i_se < pha_mat_se[i_ss].nelem(); i_se++) {
      ARTS_ASSERT(nT == pha_mat_se[i_ss][i_se].nshelves());

      for (Index Tind = 0; Tind < nT; Tind++) {
        if (pnds(i_se_flat, Tind) != 0.) {
          if (t_ok(i_se_flat, Tind) > 0.) {
            pha_tmp =
                pha_mat_se[i_ss][i_se](joker, Tind, joker, joker, joker, joker);
            pha_tmp *= pnds(i_se_flat, Tind);
            pha_mat[i_ss](joker, Tind, joker, joker, joker, joker) += pha_tmp;
          } else {
            ARTS_USER_ERROR (
              "Interpolation error for (flat-array) scattering element #",
              i_se_flat, "\n"
              "at location/temperature point #", Tind, "\n")
          }
        }
      }
      i_se_flat++;
    }
    ptype[i_ss] = max(ptypes_se[i_ss]);
  }
}

//! Phase matrices from all scattering elements.
/*!
  Derives temperature and direction interpolated phase matrices for all
  scattering elements present in scat_data.

  ATTENTION:
  If scat_data has only one freq point, f_index=-1 (i.e. all) extracts only this
  one freq point. To duplicate that as needed if f_grid has more freqs is TASK
  of the CALLING METHOD!

  Loops over pha_mat_1ScatElem and packs its output into all-scat-elements
  containers.

  \param[out] pha_mat    Phase matrix (over scat elements, freq, temp,
                           propagation direction, incident direction).
  \param[out] ptypes     Scattering element types.
  \param[out] t_ok       Flag whether T-interpol valid (over scat elements, temp).
  \param[in]  scat_data  as the WSV.
  \param[in]  stokes_dim as the WSV.
  \param[in]  T_array    Temperatures to extract pha for.
  \param[in]  pdir_array Propagation directions to extract pha for (as pairs of
                           zenith and azimuth angle per direction).
  \param[in]  idir_array Incident directions to extract pha for (as pairs of
                           zenith and azimuth angle per direction).
  \param[in]  f_index    Index of frequency to extract. -1 extracts data for all
                           freqs available in ssd.
  \param[in]  t_interp_order  Temperature interpolation order.

  \author Jana Mendrok
  \date   2018-03-24
*/
void pha_mat_NScatElems(             //Output
    ArrayOfArrayOfTensor6& pha_mat,  // [nss][nse](nf,nT,npdir,nidir,nst,nst)
    ArrayOfArrayOfIndex& ptypes,
    Matrix& t_ok,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& stokes_dim,
    const Vector& T_array,
    const Matrix& pdir_array,
    const Matrix& idir_array,
    const Index& f_index,
    const Index& t_interp_order) {
  Index f_start, nf;
  if (f_index < 0) {
    nf = scat_data[0][0].pha_mat_data.nlibraries();
    f_start = 0;
    //f_end = f_start+nf;
  } else {
    nf = 1;
    if (scat_data[0][0].pha_mat_data.nlibraries() == 1)
      f_start = 0;
    else
      f_start = f_index;
    //f_end = f_start+nf;
  }

  const Index nT = T_array.nelem();
  const Index npDir = pdir_array.nrows();
  const Index niDir = idir_array.nrows();

  const Index nss = scat_data.nelem();
  pha_mat.resize(nss);
  ptypes.resize(nss);

  const Index Nse_all = TotalNumberOfElements(scat_data);
  t_ok.resize(Nse_all, nT);
  Index i_se_flat = 0;

  for (Index i_ss = 0; i_ss < nss; i_ss++) {
    Index nse = scat_data[i_ss].nelem();
    pha_mat[i_ss].resize(nse);
    ptypes[i_ss].resize(nse);

    for (Index i_se = 0; i_se < nse; i_se++) {
      pha_mat[i_ss][i_se].resize(nf, nT, npDir, niDir, stokes_dim, stokes_dim);

      pha_mat_1ScatElem(pha_mat[i_ss][i_se],
                        ptypes[i_ss][i_se],
                        t_ok(i_se_flat, joker),
                        scat_data[i_ss][i_se],
                        T_array,
                        pdir_array,
                        idir_array,
                        f_start,
                        t_interp_order);
      i_se_flat++;
    }
  }
  ARTS_ASSERT(i_se_flat == Nse_all);
}

//! Preparing phase matrix from one scattering element.
/*!
  Extracts and prepares phase matrix data for one scattering element for one or
  all frequencies from the single scattering data. Includes interpolation in
  temperature as well as in incident and to propagation direction. Handles
  multiple output temperatures, propagation and incident directions at a time.
  Temperature interpolation order can be chosen.

  \param[out] pha_mat    1-scattering element phase matrix (over freq, temp,
                           propagation dir, incident dir).
  \param[out] ptype      Type of scattering element.
  \param[out] t_ok       Flag whether T-interpol valid (length of T_array).
  \param[in]  ssd        Single scattering data of one scattering element.
  \param[in]  T_array    Temperatures to extract pha for.
  \param[in]  pdir_array Propagation directions to extract pha for (as pairs of
                           zenith and azimuth angle per direction).
  \param[in]  idir_array Inident directions to extract pha for (as pairs of
                           zenith and azimuth angle per direction).
  \param[in]  f_start    Start index of frequency/ies to extract.
  \param[in]  t_interp_order  Temperature interpolation order.

  \author Jana Mendrok
  \date   2018-03-22
*/
void pha_mat_1ScatElem(   //Output
    Tensor6View pha_mat,  // nf, nT, npdir, nidir, nst, nst
    Index& ptype,
    VectorView t_ok,
    //Input
    const SingleScatteringData& ssd,
    const Vector& T_array,
    const Matrix& pdir_array,
    const Matrix& idir_array,
    const Index& f_start,
    const Index& t_interp_order) {
  ARTS_ASSERT(ssd.ptype == PTYPE_TOTAL_RND or ssd.ptype == PTYPE_AZIMUTH_RND);

  const Index nf = pha_mat.nvitrines();
  if (nf > 1) {
    ARTS_ASSERT(nf == ssd.f_grid.nelem());
  }

  const Index nTout = T_array.nelem();
  ARTS_ASSERT(pha_mat.nshelves() == nTout);
  ARTS_ASSERT(t_ok.nelem() == nTout);

  const Index npDir = pdir_array.nrows();
  ARTS_ASSERT(pha_mat.nbooks() == npDir);

  const Index niDir = idir_array.nrows();
  ARTS_ASSERT(pha_mat.npages() == niDir);

  const Index stokes_dim = pha_mat.ncols();
  ARTS_ASSERT(pha_mat.nrows() == stokes_dim);

  ptype = ssd.ptype;

  // Determine T-interpol order as well as interpol positions and weights (they
  // are the same for all directions (and freqs), ie it is sufficient to
  // calculate them once).
  const Index nTin = ssd.T_grid.nelem();
  Index this_T_interp_order;
  const auto T_lag = ssd_tinterp_parameters(t_ok,
                                            this_T_interp_order,
                                            ssd.T_grid,
                                            T_array,
                                            t_interp_order);
  const auto T_itw_lag = interpweights(T_lag);

  // Now loop over requested directions (and apply simultaneously for all freqs):
  // 1) interpolate direction
  // 2) apply T-interpol
  if (ptype == PTYPE_TOTAL_RND) {
    // determine how many of the compact stokes elements we will need.
    // restrict interpolations to those.
    Index npha;
    if (stokes_dim == 1)
      npha = 1;
    else if (stokes_dim < 4)  // stokes_dim==2 || stokes_dim==3
      npha = 4;
    else
      npha = 6;
    if (this_T_interp_order < 0)  // just extract (and unpack) pha data for the
    // fs and Tin, and interpolate in sca-angs, and
    // sort (copy) into the output fs, Ts, and dirs.
    {
      for (Index pdir = 0; pdir < npDir; pdir++)
        for (Index idir = 0; idir < niDir; idir++) {
          // calc scat ang theta from incident and prop dirs
          Numeric theta = scat_angle(pdir_array(pdir, 0),
                                     pdir_array(pdir, 1),
                                     idir_array(idir, 0),
                                     idir_array(idir, 1));

          // get scat angle interpolation weights
          GridPos dir_gp;
          gridpos(dir_gp, ssd.za_grid, theta * RAD2DEG);
          Vector dir_itw(2);
          interpweights(dir_itw, dir_gp);

          Vector pha_mat_int(npha, 0.);
          Matrix pha_mat_tmp(stokes_dim, stokes_dim);
          for (Index find = 0; find < nf; find++) {
            // perform the scat angle interpolation
            for (Index nst = 0; nst < npha; nst++)
              pha_mat_int[nst] = interp(
                  dir_itw,
                  ssd.pha_mat_data(find + f_start, 0, joker, 0, 0, 0, nst),
                  dir_gp);

            // convert from scat to lab frame
            pha_mat_labCalc(pha_mat_tmp(joker, joker),
                            pha_mat_int,
                            pdir_array(pdir, 0),
                            pdir_array(pdir, 1),
                            idir_array(idir, 0),
                            idir_array(idir, 1),
                            theta);

            for (Index Tind = 0; Tind < nTout; Tind++)
              pha_mat(find, Tind, pdir, idir, joker, joker) = pha_mat_tmp;
          }
        }
    } else  // T-interpolation required. To be done on the compact ssd format.
    {
      for (Index pdir = 0; pdir < npDir; pdir++)
        for (Index idir = 0; idir < niDir; idir++) {
          // calc scat ang theta from incident and prop dirs
          Numeric theta = scat_angle(pdir_array(pdir, 0),
                                     pdir_array(pdir, 1),
                                     idir_array(idir, 0),
                                     idir_array(idir, 1));

          // get scat angle interpolation weights
          GridPos dir_gp;
          gridpos(dir_gp, ssd.za_grid, theta * RAD2DEG);
          Vector dir_itw(2);
          interpweights(dir_itw, dir_gp);

          Matrix pha_mat_int(nTin, npha, 0.);
          Matrix pha_mat_tmp(nTout, npha, 0.);
          for (Index find = 0; find < nf; find++) {
            for (Index Tind = 0; Tind < nTin; Tind++)
              // perform the scat angle interpolation
              for (Index nst = 0; nst < npha; nst++) {
                pha_mat_int(Tind, nst) = interp(
                    dir_itw,
                    ssd.pha_mat_data(find + f_start, Tind, joker, 0, 0, 0, nst),
                    dir_gp);
              }
            // perform the T-interpolation
            for (Index nst = 0; nst < npha; nst++) {
              reinterp(pha_mat_tmp(joker, nst),
                       pha_mat_int(joker, nst),
                       T_itw_lag,
                       T_lag);
            }
            // FIXME: it's probably better to do the frame conversion first,
            // then the T-interpolation (which is better depends on how many T-
            // aka vertical grid points we have. for single point, T-interpol
            // first is better, for a full vertical grid, frame conversion first
            // should be faster.)
            for (Index Tind = 0; Tind < nTout; Tind++) {
              // convert from scat to lab frame
              pha_mat_labCalc(pha_mat(find, Tind, pdir, idir, joker, joker),
                              pha_mat_tmp(Tind, joker),
                              pdir_array(pdir, 0),
                              pdir_array(pdir, 1),
                              idir_array(idir, 0),
                              idir_array(idir, 1),
                              theta);
            }
          }
        }
    }
  } else  // dir-interpolation for non-tot-random particles
          // Data is already stored in the laboratory frame,
          // but it is compressed a little. Details elsewhere.
  {
    Index nDir = npDir * niDir;
    Vector adelta_aa(nDir);
    Matrix delta_aa(npDir, niDir);
    ArrayOfGridPos daa_gp(nDir), pza_gp(nDir), iza_gp(nDir);
    ArrayOfGridPos pza_gp_tmp(npDir), iza_gp_tmp(niDir);

    gridpos(pza_gp_tmp, ssd.za_grid, pdir_array(joker, 0));
    gridpos(iza_gp_tmp, ssd.za_grid, idir_array(joker, 0));

    Index j = 0;
    for (Index pdir = 0; pdir < npDir; pdir++) {
      for (Index idir = 0; idir < niDir; idir++) {
        delta_aa(pdir, idir) =
            pdir_array(pdir, 1) - idir_array(idir, 1) +
            (pdir_array(pdir, 1) - idir_array(idir, 1) < -180) * 360 -
            (pdir_array(pdir, 1) - idir_array(idir, 1) > 180) * 360;
        adelta_aa[j] = abs(delta_aa(pdir, idir));
        pza_gp[j] = pza_gp_tmp[pdir];
        iza_gp[j] = iza_gp_tmp[idir];
        j++;
      }
    }

    gridpos(daa_gp, ssd.aa_grid, adelta_aa);

    Matrix dir_itw(nDir, 8);
    interpweights(dir_itw, pza_gp, daa_gp, iza_gp);

    if (this_T_interp_order < 0)  // T only needs to be extracted.
    {
      Tensor3 pha_mat_int(nDir, stokes_dim, stokes_dim, 0.);
      Tensor4 pha_mat_tmp(npDir, niDir, stokes_dim, stokes_dim);

      for (Index find = 0; find < nf; find++) {
        // perform the (tri-linear) angle interpolation. but only for the
        // pha_mat elements that we actually need.

        for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
          for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
            interp(
                pha_mat_int(joker, ist1, ist2),
                dir_itw,
                ssd.pha_mat_data(
                    find + f_start, 0, joker, joker, joker, 0, ist1 * 4 + ist2),
                pza_gp,
                daa_gp,
                iza_gp);

        // sort direction-combined 1D-array back into prop and incident
        // direction matrix
        Index i = 0;
        for (Index pdir = 0; pdir < npDir; pdir++)
          for (Index idir = 0; idir < niDir; idir++) {
            pha_mat_tmp(pdir, idir, joker, joker) =
                pha_mat_int(i, joker, joker);
            i++;
          }

        if (stokes_dim > 2) {
          for (Index pdir = 0; pdir < npDir; pdir++)
            for (Index idir = 0; idir < niDir; idir++)
              if (delta_aa(pdir, idir) < 0.) {
                pha_mat_tmp(pdir, idir, 0, 2) *= -1;
                pha_mat_tmp(pdir, idir, 1, 2) *= -1;
                pha_mat_tmp(pdir, idir, 2, 0) *= -1;
                pha_mat_tmp(pdir, idir, 2, 1) *= -1;
              }
        }

        if (stokes_dim > 3) {
          for (Index pdir = 0; pdir < npDir; pdir++)
            for (Index idir = 0; idir < niDir; idir++)
              if (delta_aa(pdir, idir) < 0.) {
                pha_mat_tmp(pdir, idir, 0, 3) *= -1;
                pha_mat_tmp(pdir, idir, 1, 3) *= -1;
                pha_mat_tmp(pdir, idir, 3, 0) *= -1;
                pha_mat_tmp(pdir, idir, 3, 1) *= -1;
              }
        }

        for (Index Tind = 0; Tind < nTout; Tind++)
          pha_mat(find, Tind, joker, joker, joker, joker) = pha_mat_tmp;
      }
    }

    else  // T- and dir-interpolation required. To be done on the compact ssd
          // format.
    {
      Tensor4 pha_mat_int(nTin, nDir, stokes_dim, stokes_dim, 0.);

      for (Index find = 0; find < nf; find++) {
        // perform the (tri-linear) angle interpolation. but only for the
        // pha_mat elements that we actually need.
        for (Index Tind = 0; Tind < nTin; Tind++) {
          for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
            for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
              interp(pha_mat_int(Tind, joker, ist1, ist2),
                     dir_itw,
                     ssd.pha_mat_data(find + f_start,
                                      Tind,
                                      joker,
                                      joker,
                                      joker,
                                      0,
                                      ist1 * 4 + ist2),
                     pza_gp,
                     daa_gp,
                     iza_gp);
        }

        // perform the T-interpolation and simultaneously sort
        // direction-combined 1D-array back into prop and incident direction
        // matrix.
        Index i = 0;
        for (Index pdir = 0; pdir < npDir; pdir++)
          for (Index idir = 0; idir < niDir; idir++) {
            for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
              for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
                reinterp(pha_mat(find, joker, pdir, idir, ist1, ist2),
                         pha_mat_int(joker, i, ist1, ist2),
                         T_itw_lag,
                         T_lag);
            i++;
          }

        if (stokes_dim > 2) {
          for (Index pdir = 0; pdir < npDir; pdir++)
            for (Index idir = 0; idir < niDir; idir++)
              if (delta_aa(pdir, idir) < 0.) {
                pha_mat(find, joker, pdir, idir, 0, 2) *= -1;
                pha_mat(find, joker, pdir, idir, 1, 2) *= -1;
                pha_mat(find, joker, pdir, idir, 2, 0) *= -1;
                pha_mat(find, joker, pdir, idir, 2, 1) *= -1;
              }
        }

        if (stokes_dim > 3) {
          for (Index pdir = 0; pdir < npDir; pdir++)
            for (Index idir = 0; idir < niDir; idir++)
              if (delta_aa(pdir, idir) < 0.) {
                pha_mat(find, joker, pdir, idir, 0, 3) *= -1;
                pha_mat(find, joker, pdir, idir, 1, 3) *= -1;
                pha_mat(find, joker, pdir, idir, 3, 0) *= -1;
                pha_mat(find, joker, pdir, idir, 3, 1) *= -1;
              }
        }
      }
    }
  }
}

//! Preparing phase matrix fourier series components for one scattering element.
/*!
  Calculates phase matrix fourier series components (currently only mode 0) for
  one scattering element for one or all frequencies from the single scattering
  data. Includes interpolation in temperature as well as in incident and
  scattered stream directions. Handles multiple output temperatures, propagation
  and incident directions at a time.
  Temperature interpolation order can be chosen.

  \param[out] pha_mat_fou  1-scattering element phase matrix Fourier components
                           (over freq, temp, propagation dir, incident dir).
  \param[out] ptype      Type of scattering element.
  \param[out] t_ok       Flag whether T-interpol valid (length of T_array).
  \param[in]  ssd        Single scattering data of one scattering element.
  \param[in]  T_array    Temperatures to extract pha for.
  \param[in]  pdir_array  Propagation directions (polar only) to extract Fourier
                            components for.
  \param[in]  idir_array  Incident directions (polar only) to extract Fourier
                            components for.
  \param[in]  f_start    Start index of frequency/ies to extract.
  \param[in]  t_interp_order  Temperature interpolation order.
  \param[in]  naa_totran Number of (equidistant) azimuth directions to consider
                           in calculation of Fourier series for totally random
                           orientation elements.

  \author Jana Mendrok
  \date   2018-05-01
*/
void FouComp_1ScatElem(       //Output
    Tensor7View pha_mat_fou,  // nf, nT, npdir, nidir, nst, nst, m
    Index& ptype,
    VectorView t_ok,
    //Input
    const SingleScatteringData& ssd,
    const Vector& T_array,
    const Vector& pdir_array,
    const Vector& idir_array,
    const Index& f_start,
    const Index& t_interp_order,
    const Index& naa_totran) {
  ARTS_ASSERT(ssd.ptype == PTYPE_TOTAL_RND or ssd.ptype == PTYPE_AZIMUTH_RND);

  const Index nf = pha_mat_fou.nlibraries();
  if (nf > 1) {
    ARTS_ASSERT(nf == ssd.f_grid.nelem());
  }

  const Index nTout = T_array.nelem();
  ARTS_ASSERT(pha_mat_fou.nvitrines() == nTout);
  ARTS_ASSERT(t_ok.nelem() == nTout);

  const Index npDir = pdir_array.nelem();
  ARTS_ASSERT(pha_mat_fou.nshelves() == npDir);
  const Index niDir = idir_array.nelem();
  ARTS_ASSERT(pha_mat_fou.nbooks() == niDir);

  const Index stokes_dim = pha_mat_fou.nrows();
  ARTS_ASSERT(pha_mat_fou.npages() == stokes_dim);
  // currently code is only prepared for stokes_dim up to 2 (nothing else needed in
  // RT4 and generally in azimuth-symmetrical system)
  ARTS_ASSERT(stokes_dim < 3);

  const Index nmodes = pha_mat_fou.ncols();
  // currently code is only prepared for fourier mode 0 (nothing else needed in
  // RT4 and generally in azimuth-symmetrical system)
  ARTS_ASSERT(nmodes == 0);

  ptype = ssd.ptype;

  // Determine T-interpol order as well as interpol positions and weights (they
  // are the same for all directions (and freqs), ie it is sufficient to
  // calculate them once).
  const Index nTin = ssd.T_grid.nelem();
  Index this_T_interp_order;
  const auto T_lag = ssd_tinterp_parameters(t_ok,
                                            this_T_interp_order,
                                            ssd.T_grid,
                                            T_array,
                                            t_interp_order);
  const auto T_itw_lag = interpweights(T_lag);

  // 1) derive Fourier component(s)
  // 2) apply T-interpol
  if (ptype == PTYPE_TOTAL_RND) {
    // DCalculate azimuth angles and their integration weights for Fourier
    // component derivation (they are only determined by naa_totran).
    ARTS_USER_ERROR_IF (naa_totran < 3,
        "Azimuth grid size for scatt matrix extraction"
        " (*naa_totran*) must be >3.\n"
        "Yours is ", naa_totran, ".\n")
    Vector aa_grid;
    nlinspace(aa_grid, 0, 180, naa_totran);
    Numeric daa_totran =
        1. / float(naa_totran - 1);  // 2*180./360./(naa_totran-1)
    Vector theta(naa_totran);
    ArrayOfGridPos theta_gp;
    Matrix theta_itw(naa_totran, 2);

    Index npha;
    if (stokes_dim == 1)
      npha = 1;
    else if (stokes_dim < 4)  // stokes_dim==2 || stokes_dim==3
      npha = 4;
    else
      npha = 6;

    Matrix pha_mat(stokes_dim, stokes_dim);
    Matrix pha_mat_angint(naa_totran, npha);
    Tensor3 Fou_int(nTin, stokes_dim, stokes_dim);

    for (Index idir = 0; idir < niDir; idir++)
      for (Index pdir = 0; pdir < npDir; pdir++) {
        for (Index iaa = 0; iaa < naa_totran; iaa++)
          // calc scat ang theta from incident and prop dirs
          theta[iaa] =
              scat_angle(pdir_array[pdir], aa_grid[iaa], idir_array[idir], 0.);

        // get scat angle interpolation weights
        gridpos(theta_gp, ssd.za_grid, theta * RAD2DEG);
        interpweights(theta_itw, theta_gp);

        for (Index find = 0; find < nf; find++) {
          Fou_int = 0.;
          for (Index Tind = 0; Tind < nTin; Tind++) {
            // perform the scat angle interpolation
            for (Index nst = 0; nst < npha; nst++)
              interp(
                  pha_mat_angint(joker, nst),
                  theta_itw,
                  ssd.pha_mat_data(find + f_start, Tind, joker, 0, 0, 0, nst),
                  theta_gp);
            for (Index iaa = 0; iaa < naa_totran; iaa++) {
              // convert from scat to lab frame
              pha_mat_labCalc(pha_mat,
                              pha_mat_angint(iaa, joker),
                              pdir_array[pdir],
                              aa_grid[iaa],
                              idir_array[idir],
                              0.,
                              theta[iaa]);
              // and sum up/integrate
              // FIXME: can the integration probably be done in lab frame?
              // test!
              if (iaa == 0 || iaa == naa_totran - 1)
                pha_mat *= (daa_totran / 2.);
              else
                pha_mat *= daa_totran;
              Fou_int(Tind, joker, joker) += pha_mat;
            }
          }

          if (this_T_interp_order <
              0)  // T only needs to be sorted into pha_mat_fou.
          {
            for (Index Tind = 0; Tind < nTout; Tind++)
              pha_mat_fou(find, Tind, pdir, idir, joker, joker, 0) =
                  Fou_int(0, joker, joker);
          } else {
            for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
              for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
                for (Index im = 0; im < nmodes; im++)
                  reinterp(pha_mat_fou(find, joker, pdir, idir, ist1, ist2, im),
                           Fou_int(joker, ist1, ist2),
                           T_itw_lag,
                           T_lag);
          }
        }
      }
  } else  // derive Fourier component(s) on scattering elements own za_grid (using
          // given aa_grid), then 2D-interpolate inc and sca polar directions.
  // Do calcs only for required pha_mat components depending on stokes_dim.
  {
    Index nza = ssd.za_grid.nelem();
    Index naa = ssd.aa_grid.nelem();
    ConstVectorView za_datagrid = ssd.za_grid;
    ConstVectorView aa_datagrid = ssd.aa_grid;
    ARTS_ASSERT(aa_datagrid[0] == 0.);
    ARTS_ASSERT(aa_datagrid[naa - 1] == 180.);
    Vector daa(naa);

    // Precalculate azimuth integration weights for this azimuthally randomly
    // oriented scat element (need to do this per scat element as ssd.aa_grid is
    // scat element specific (might change between elements) and need to do this
    // on actual grid instead of grid number since the grid can, at least
    // theoretically be non-equidistant).
    daa[0] = (aa_datagrid[1] - aa_datagrid[0]) / 360.;
    for (Index iaa = 1; iaa < naa - 1; iaa++)
      daa[iaa] = (aa_datagrid[iaa + 1] - aa_datagrid[iaa - 1]) / 360.;
    daa[naa - 1] = (aa_datagrid[naa - 1] - aa_datagrid[naa - 2]) / 360.;

    // Precalculate polar angle interpolation grid positions and weights
    ArrayOfGridPos pdir_za_gp, idir_za_gp;
    gridpos(pdir_za_gp, za_datagrid, pdir_array);
    gridpos(idir_za_gp, za_datagrid, idir_array);
    Tensor3 dir_itw(npDir, niDir, 4);
    interpweights(dir_itw, pdir_za_gp, idir_za_gp);

    Tensor4 Fou_ssd(nza, nza, stokes_dim, stokes_dim);
    Tensor5 Fou_angint(nTin, npDir, niDir, stokes_dim, stokes_dim);

    for (Index find = 0; find < nf; find++) {
      for (Index Tind = 0; Tind < nTin; Tind++) {
        // first, extract the phase matrix at the scatt elements own polar angle
        // grid and integrate over azimuth deriving their respective azimuthal
        // (Fourier series) 0-mode
        Fou_ssd = 0.;
        for (Index iza = 0; iza < nza; iza++)
          for (Index sza = 0; sza < nza; sza++) {
            for (Index iaa = 0; iaa < naa; iaa++) {
              // FIXME: are there any possible shortcuts for specific stokes
              // components (specifically, some where F0(ist1,ist2)==0?)
              for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
                for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
                  Fou_ssd(sza, iza, ist1, ist2) +=
                      daa[iaa] * ssd.pha_mat_data(find + f_start,
                                                  Tind,
                                                  sza,
                                                  iaa,
                                                  iza,
                                                  0,
                                                  ist1 * 4 + ist2);
            }
          }

        // second, interpolate the extracted azimuthal mode to the stream directions
        for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
          for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
            // FIXME: do we need to apply any sign changes?
            interp(Fou_angint(Tind, joker, joker, ist1, ist2),
                   dir_itw,
                   Fou_ssd(joker, joker, ist1, ist2),
                   pdir_za_gp,
                   idir_za_gp);
      }

      if (this_T_interp_order <
          0)  // T only needs to be sorted into pha_mat_fou.
      {
        for (Index Tind = 0; Tind < nTout; Tind++)
          pha_mat_fou(find, Tind, joker, joker, joker, joker, 0) =
              Fou_angint(0, joker, joker, joker, joker);
      } else {
        for (Index pdir = 0; pdir < npDir; pdir++)
          for (Index idir = 0; idir < niDir; idir++)
            for (Index ist1 = 0; ist1 < stokes_dim; ist1++)
              for (Index ist2 = 0; ist2 < stokes_dim; ist2++)
                for (Index im = 0; im < nmodes; im++)
                  reinterp(pha_mat_fou(find, joker, pdir, idir, ist1, ist2, im),
                           Fou_angint(joker, pdir, idir, ist1, ist2),
                           T_itw_lag,
                           T_lag);
      }
    }
  }
}

//! Transformation of absorption vector.
/*!
  In the single scattering database the data of the absorption vector is
  stored in different coordinate systems, depending on the type (ptype) of
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  Output and Input:
  \param abs_vec_lab Absorption vector in Laboratory frame.
  Input:
  \param abs_vec_data Absorption vector in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Type of scattering element.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.

  \author Claudia Emde
  \date   2003-05-24
*/
void abs_vecTransform(  //Output and Input
    StokesVector& abs_vec_lab,
    //Input
    ConstTensor3View abs_vec_data,
    ConstVectorView za_datagrid,
    ConstVectorView aa_datagrid _U_,
    const PType& ptype,
    const Numeric& za_sca _U_,
    const Numeric& aa_sca _U_,
    const Verbosity& verbosity) {
  const Index stokes_dim = abs_vec_lab.StokesDimensions();
  ARTS_ASSERT(abs_vec_lab.NumberOfFrequencies() == 1);

  ARTS_USER_ERROR_IF (stokes_dim > 4 || stokes_dim < 1,
        "The dimension of the stokes vector \n"
        "must be 1,2,3 or 4");

  switch (ptype) {
    case PTYPE_GENERAL: {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n";
      break;
    }
    case PTYPE_TOTAL_RND: {
      // The first element of the vector corresponds to the absorption
      // coefficient which is stored in the database, the others are 0.

      abs_vec_lab.SetZero();

      abs_vec_lab.Kjj()[0] = abs_vec_data(0, 0, 0);
      break;
    }

    case PTYPE_AZIMUTH_RND:  //Added by Cory Davis 9/12/03
    {
      ARTS_ASSERT(abs_vec_data.ncols() == 2);

      // In the case of azimuthally randomly oriented particles, only the first
      // two elements of the absorption coefficient vector are non-zero.
      // These values are dependent on the zenith angle of propagation.

      // 1st interpolate data by za_sca
      GridPos gp;
      Vector itw(2);

      gridpos(gp, za_datagrid, za_sca);
      interpweights(itw, gp);

      abs_vec_lab.SetZero();

      abs_vec_lab.Kjj()[0] = interp(itw, abs_vec_data(Range(joker), 0, 0), gp);

      if (stokes_dim == 1) {
        break;
      }
      abs_vec_lab.K12()[0] = interp(itw, abs_vec_data(Range(joker), 0, 1), gp);
      break;
    }
    default: {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }
}

//! Transformation of extinction matrix.
/*!
  In the single scattering database the data of the extinction matrix is
  stored in different coordinate systems, depending on the type (ptype) of
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  Output and Input:
  \param ext_mat_lab Extinction matrix in Laboratory frame.
  Input:
  \param ext_mat_data Extinction matrix in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Type of scattering element.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.

  \author Claudia Emde
  \date   2003-05-24
*/
void ext_matTransform(  //Output and Input
    PropagationMatrix& ext_mat_lab,
    //Input
    ConstTensor3View ext_mat_data,
    ConstVectorView za_datagrid,
    ConstVectorView aa_datagrid _U_,
    const PType& ptype,
    const Numeric& za_sca,
    const Numeric& aa_sca _U_,
    const Verbosity& verbosity) {
  const Index stokes_dim = ext_mat_lab.StokesDimensions();
  ARTS_ASSERT(ext_mat_lab.NumberOfFrequencies() == 1);

  ARTS_USER_ERROR_IF (stokes_dim > 4 || stokes_dim < 1,
        "The dimension of the stokes vector \n"
        "must be 1,2,3 or 4");

  switch (ptype) {
    case PTYPE_GENERAL: {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n";
      break;
    }
    case PTYPE_TOTAL_RND: {
      ARTS_ASSERT(ext_mat_data.ncols() == 1);

      // In the case of randomly oriented particles the extinction matrix is
      // diagonal. The value of each element of the diagonal is the
      // extinction cross section, which is stored in the database.

      ext_mat_lab.SetZero();

      ext_mat_lab.Kjj() = ext_mat_data(0, 0, 0);
      break;
    }

    case PTYPE_AZIMUTH_RND:  //Added by Cory Davis 9/12/03
    {
      ARTS_ASSERT(ext_mat_data.ncols() == 3);

      // In the case of azimuthally randomly oriented particles, the extinction
      // matrix has only 3 independent non-zero elements Kjj, K12=K21, and K34=-K43.
      // These values are dependent on the zenith angle of propagation.

      // 1st interpolate data by za_sca
      GridPos gp;
      Vector itw(2);
      Numeric Kjj;
      Numeric K12;
      Numeric K34;

      gridpos(gp, za_datagrid, za_sca);
      interpweights(itw, gp);

      ext_mat_lab.SetZero();

      Kjj = interp(itw, ext_mat_data(Range(joker), 0, 0), gp);
      ext_mat_lab.Kjj() = Kjj;

      if (stokes_dim < 2) {
        break;
      }

      K12 = interp(itw, ext_mat_data(Range(joker), 0, 1), gp);
      ext_mat_lab.K12()[0] = K12;

      if (stokes_dim < 4) {
        break;
      }

      K34 = interp(itw, ext_mat_data(Range(joker), 0, 2), gp);
      ext_mat_lab.K34()[0] = K34;
      break;
    }
    default: {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }
}

//! Transformation of phase matrix.
/*!
  In the single scattering database the data of the phase matrix is
  stored in different coordinate systems, depending on the type (ptype) of
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  \param[out,in] pha_mat_lab   Phase matrix in Laboratory frame.
  \param[in]     pha_mat_data  Phase matrix in database.
  \param[in]     za_datagrid   Zenith angle grid in the database.
  \param[in]     aa_datagrid   Zenith angle grid in the database.
  \param[in]     ptype Type of scattering element.
  \param[in]     za_sca_idx    Index of zenith angle of scattered direction
                                 within za_grid.
  \param[in]     aa_sca_idx    Index of azimuth angle of scattered direction
                                 within aa_grid.
  \param[in]     za_inc_idx    Index of zenith angle of incoming direction
                                 within za_grid.
  \param[in]     aa_inc_idx    Index of azimuth angle of incoming direction
                                 within aa_grid.
  \param[in]     za_grid  Grid of zenith angles to extract pha_mat for.
  \param[in]     aa_grid  Grid of azimuth angles to extract pha_mat for.

  \author Claudia Emde
  \date   2003-08-19
*/
void pha_matTransform(  //Output
    MatrixView pha_mat_lab,
    //Input
    ConstTensor5View pha_mat_data,
    ConstVectorView za_datagrid,
    ConstVectorView aa_datagrid,
    const PType& ptype,
    const Index& za_sca_idx,
    const Index& aa_sca_idx,
    const Index& za_inc_idx,
    const Index& aa_inc_idx,
    ConstVectorView za_grid,
    ConstVectorView aa_grid,
    const Verbosity& verbosity) {
  const Index stokes_dim = pha_mat_lab.ncols();

  Numeric za_sca = za_grid[za_sca_idx];
  Numeric aa_sca = aa_grid[aa_sca_idx];
  Numeric za_inc = za_grid[za_inc_idx];
  Numeric aa_inc = aa_grid[aa_inc_idx];

  ARTS_USER_ERROR_IF (stokes_dim > 4 || stokes_dim < 1,
        "The dimension of the stokes vector \n"
        "must be 1,2,3 or 4");

  switch (ptype) {
    case PTYPE_GENERAL: {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n";
      break;
    }
    case PTYPE_TOTAL_RND: {
      // Calculate the scattering angle and interpolate the data in it:
      Numeric theta_rad = scat_angle(za_sca, aa_sca, za_inc, aa_inc);
      const Numeric theta = RAD2DEG * theta_rad;

      // Interpolation of the data on the scattering angle:
      Vector pha_mat_int(6);
      interpolate_scat_angle(pha_mat_int, pha_mat_data, za_datagrid, theta);

      // Calculate the phase matrix in the laboratory frame:
      pha_mat_labCalc(
          pha_mat_lab, pha_mat_int, za_sca, aa_sca, za_inc, aa_inc, theta_rad);

      break;
    }

    case PTYPE_AZIMUTH_RND:  //Added by Cory Davis
                             //Data is already stored in the laboratory frame,
      //but it is compressed a little.  Details elsewhere.
      {
        ARTS_ASSERT(pha_mat_data.ncols() == 16);
        ARTS_ASSERT(pha_mat_data.npages() == za_datagrid.nelem());
        Numeric delta_aa = aa_sca - aa_inc + (aa_sca - aa_inc < -180) * 360 -
                           (aa_sca - aa_inc > 180) *
                               360;  //delta_aa corresponds to the "books"
                                     //dimension of pha_mat_data
        GridPos za_sca_gp;
        GridPos delta_aa_gp;
        GridPos za_inc_gp;
        Vector itw(8);

        gridpos(delta_aa_gp, aa_datagrid, abs(delta_aa));
        gridpos(za_inc_gp, za_datagrid, za_inc);
        gridpos(za_sca_gp, za_datagrid, za_sca);

        interpweights(itw, za_sca_gp, delta_aa_gp, za_inc_gp);

        pha_mat_lab(0, 0) =
            interp(itw,
                   pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 0),
                   za_sca_gp,
                   delta_aa_gp,
                   za_inc_gp);
        if (stokes_dim == 1) {
          break;
        }
        pha_mat_lab(0, 1) =
            interp(itw,
                   pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 1),
                   za_sca_gp,
                   delta_aa_gp,
                   za_inc_gp);
        pha_mat_lab(1, 0) =
            interp(itw,
                   pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 4),
                   za_sca_gp,
                   delta_aa_gp,
                   za_inc_gp);
        pha_mat_lab(1, 1) =
            interp(itw,
                   pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 5),
                   za_sca_gp,
                   delta_aa_gp,
                   za_inc_gp);
        if (stokes_dim == 2) {
          break;
        }
        if (delta_aa >= 0) {
          pha_mat_lab(0, 2) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 2),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(1, 2) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 6),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(2, 0) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 8),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(2, 1) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 9),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
        } else {
          pha_mat_lab(0, 2) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 2),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(1, 2) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 6),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(2, 0) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 8),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(2, 1) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 9),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
        }
        pha_mat_lab(2, 2) = interp(
            itw,
            pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 10),
            za_sca_gp,
            delta_aa_gp,
            za_inc_gp);
        if (stokes_dim == 3) {
          break;
        }
        if (delta_aa >= 0) {
          pha_mat_lab(0, 3) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 3),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(1, 3) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 7),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(3, 0) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 12),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(3, 1) = interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 13),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
        } else {
          pha_mat_lab(0, 3) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 3),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(1, 3) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 7),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(3, 0) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 12),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
          pha_mat_lab(3, 1) = -interp(
              itw,
              pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 13),
              za_sca_gp,
              delta_aa_gp,
              za_inc_gp);
        }
        pha_mat_lab(2, 3) = interp(
            itw,
            pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 11),
            za_sca_gp,
            delta_aa_gp,
            za_inc_gp);
        pha_mat_lab(3, 2) = interp(
            itw,
            pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 14),
            za_sca_gp,
            delta_aa_gp,
            za_inc_gp);
        pha_mat_lab(3, 3) = interp(
            itw,
            pha_mat_data(Range(joker), Range(joker), Range(joker), 0, 15),
            za_sca_gp,
            delta_aa_gp,
            za_inc_gp);
        break;
      }

    default: {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }
}

//! Derive extinction matrix from absorption vector.
/*!
  In case, when only absorption of a scattering element shall be considered, and
  the scattering is neglected, the extinction matrix is set from the absorption
  vector only.

  Extinction matrix is set the following way:

  K11 = K22 = K33 = K44 = a1
  K12 = K21 = a2
  K13 = K31 = a3
  K14 = K41 = a4

  Other elements remain 0.
  However, note that the other elements might be supposed to contain non-zero
  values as well. We couldn't find an appropriate solution yet, though, and
  these elements are expected to be rather small. Also, using the absorption
  part only is anyway only a rough approximation. Hence, we deem the assumption
  of setting the remaining elements to zero reasonable.

  Output and Input:
  \param ext_mat     Extinction matrix.
  Input:
  \param abs_vec     Absorption vector.
  \param stokes_dim  as the WSV.

  \author Jana Mendrok
  \date   2013-04-30
*/
void ext_matFromabs_vec(  //Output
    MatrixView ext_mat,
    //Input
    ConstVectorView abs_vec,
    const Index& stokes_dim) {
  ARTS_ASSERT(stokes_dim >= 1 && stokes_dim <= 4);
  ARTS_ASSERT(ext_mat.nrows() == stokes_dim);
  ARTS_ASSERT(ext_mat.ncols() == stokes_dim);
  ARTS_ASSERT(abs_vec.nelem() == stokes_dim);

  // first: diagonal elements
  for (Index is = 0; is < stokes_dim; is++) {
    ext_mat(is, is) += abs_vec[0];
  }
  // second: off-diagonal elements, namely first row and column
  for (Index is = 1; is < stokes_dim; is++) {
    ext_mat(0, is) += abs_vec[is];
    ext_mat(is, 0) += abs_vec[is];
  }
}

//! Calculates the scattering angle.
/*!
  The scattering angle is calculated from the angles defining
  the directions of the incoming and scattered radiation.

  \param[in]  za_sca       Zenith angle of scattered direction [deg].
  \param[in]  aa_sca       Azimuth angle of scattered direction [deg].
  \param[in]  za_inc       Zenith angle of incoming direction [deg].
  \param[in]  aa_inc       Azimuth angle of incoming direction [deg].
  \return Scattering angle [rad].

  \author Jana Mendrok (moved out from interpolate_scat_angle by C.Emde)
  \date   2018-03-23
*/
Numeric scat_angle(const Numeric& za_sca,
                   const Numeric& aa_sca,
                   const Numeric& za_inc,
                   const Numeric& aa_inc) {
  Numeric theta_rad;
  Numeric ANG_TOL = 1e-7;

  // CPD 5/10/03.
  // Two special cases are implemented here to avoid NaNs that can sometimes
  // occur in in the acos() formula in forward and backscattering cases.
  //
  // GH 2011-05-31
  // Consider not only aa_sca-aa_inc ~= 0, but also aa_sca-aa_inc ~= 360.

  if ((abs(aa_sca - aa_inc) < ANG_TOL) ||
      (abs(abs(aa_sca - aa_inc) - 360) < ANG_TOL)) {
    theta_rad = Conversion::deg2rad(abs(za_sca - za_inc));
  } else if (abs(abs(aa_sca - aa_inc) - 180) < ANG_TOL) {
    theta_rad = Conversion::deg2rad(za_sca + za_inc);
    if (theta_rad > PI) {
      theta_rad = 2 * PI - theta_rad;
    }
  } else {theta_rad =
      acos(Conversion::cosd(za_sca) * Conversion::cosd(za_inc) +
           Conversion::sind(za_sca) * Conversion::sind(za_inc) *
           Conversion::cosd(aa_sca - aa_inc));
  }
  return theta_rad;
}

//! Interpolate data on the scattering angle.
/*!
  This function is used for the transformation of the phase matrix
  from scattering frame to the laboratory frame for randomly oriented
  scattering media (case PTYPE_TOTAL_RND).

  The scattering angle is calculated from the angles defining
  the directions of the incoming and scattered radiation.
  After that the data (which is stored in the data files as a function
  of the scattering angle) is interpolated on the calculated
  scattering angle.

  \param[out] pha_mat_int  Interpolated phase matrix.
  \param[out] theta_rad    Scattering angle [rad].
  \param[in]  pha_mat_data Phase matrix in database.
  \param[in]  za_datagrid  Zenith angle grid in the database.
  \param[in]  za_sca       Zenith angle of scattered direction [rad].
  \param[in]  aa_sca       Azimuth angle of scattered direction [rad].
  \param[in]  za_inc       Zenith angle of incoming direction [rad].
  \param[in]  aa_inc       Azimuth angle of incoming direction [rad].

  \author Claudia Emde
  \date   2003-08-19
*/
void interpolate_scat_angle(  //Output:
    VectorView pha_mat_int,
    //Input:
    ConstTensor5View pha_mat_data,
    ConstVectorView za_datagrid,
    const Numeric theta) {
  GridPos thet_gp;
  gridpos(thet_gp, za_datagrid, theta);

  Vector itw(2);
  interpweights(itw, thet_gp);

  ARTS_ASSERT(pha_mat_data.ncols() == 6);
  for (Index i = 0; i < 6; i++) {
    pha_mat_int[i] = interp(itw, pha_mat_data(joker, 0, 0, 0, i), thet_gp);
  }
}

//! Calculate phase matrix in laboratory coordinate system.
/*!
  Transformation function for the phase matrix for the case of
  randomly oriented particles (case PTYPE_TOTAL_RND).

  Some of the formulas can be found in

  Mishchenkho: "Scattering, Absorption and Emission of Light
  by Small Particles", Cambridge University Press, 2002
  Capter 4

  The full set of formulas will be documented in AUG.

  Output and Input:
  \param pha_mat_lab Phase matrix in laboratory frame.
  Input:
  \param pha_mat_int Interpolated phase matrix.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
  \param za_inc Zenith angle of incoming direction.
  \param aa_inc Azimuth angle of incoming direction.
  \param theta_rad Scattering angle [rad].

  \author Claudia Emde
  \date   2003-05-13
*/
void pha_mat_labCalc(  //Output:
    MatrixView pha_mat_lab,
    //Input:
    ConstVectorView pha_mat_int,
    const Numeric& za_sca,
    const Numeric& aa_sca,
    const Numeric& za_inc,
    const Numeric& aa_inc,
    const Numeric& theta_rad) {
  const Index stokes_dim = pha_mat_lab.ncols();

  ARTS_USER_ERROR_IF (std::isnan(F11),
        "NaN value(s) detected in *pha_mat_labCalc* (0,0). Could the "
        "input data contain NaNs? Please check with *scat_dataCheck*. If "
        "input data are OK and you critically need the ongoing calculations, "
        "try to change the observation LOS slightly. If you can reproduce "
        "this error, please contact Patrick in order to help tracking down "
        "the reason to this problem. If you see this message occasionally "
        "when doing MC calculations, it should not be critical. This path "
        "sampling will be rejected and replaced with a new one.");

  // For stokes_dim = 1, we only need Z11=F11:
  pha_mat_lab(0, 0) = F11;

  if (stokes_dim > 1) {
    Numeric za_sca_rad = za_sca * DEG2RAD;
    Numeric za_inc_rad = za_inc * DEG2RAD;
    Numeric aa_sca_rad = aa_sca * DEG2RAD;
    Numeric aa_inc_rad = aa_inc * DEG2RAD;

    const Numeric ANGTOL_RAD = 1e-6;  //CPD: this constant is used to adjust
        //zenith angles close to 0 and PI.  This is
        //also used to avoid float == float statements.

    //
    // Several cases have to be considered:
    //

    if ((abs(theta_rad) < ANGTOL_RAD)          // forward scattering
        || (abs(theta_rad - PI) < ANGTOL_RAD)  // backward scattering
        ||
        (abs(aa_inc_rad - aa_sca_rad) < ANGTOL_RAD)  // inc and sca on meridian
        || (abs(abs(aa_inc_rad - aa_sca_rad) - 360.) < ANGTOL_RAD)  //   "
        || (abs(abs(aa_inc_rad - aa_sca_rad) - 180.) < ANGTOL_RAD)  //   "
    ) {
      pha_mat_lab(0, 1) = F12;
      pha_mat_lab(1, 0) = F12;
      pha_mat_lab(1, 1) = F22;

      if (stokes_dim > 2) {
        pha_mat_lab(0, 2) = 0;
        pha_mat_lab(1, 2) = 0;
        pha_mat_lab(2, 0) = 0;
        pha_mat_lab(2, 1) = 0;
        pha_mat_lab(2, 2) = F33;

        if (stokes_dim > 3) {
          pha_mat_lab(0, 3) = 0;
          pha_mat_lab(1, 3) = 0;
          pha_mat_lab(2, 3) = F34;
          pha_mat_lab(3, 0) = 0;
          pha_mat_lab(3, 1) = 0;
          pha_mat_lab(3, 2) = -F34;
          pha_mat_lab(3, 3) = F44;
        }
      }
    }

    else {
      Numeric sigma1;
      Numeric sigma2;

      Numeric s1, s2;

      // In these cases we have to take limiting values.

      if (za_inc_rad < ANGTOL_RAD) {
        sigma1 = PI + aa_sca_rad - aa_inc_rad;
        sigma2 = 0;
      } else if (za_inc_rad > PI - ANGTOL_RAD) {
        sigma1 = aa_sca_rad - aa_inc_rad;
        sigma2 = PI;
      } else if (za_sca_rad < ANGTOL_RAD) {
        sigma1 = 0;
        sigma2 = PI + aa_sca_rad - aa_inc_rad;
      } else if (za_sca_rad > PI - ANGTOL_RAD) {
        sigma1 = PI;
        sigma2 = aa_sca_rad - aa_inc_rad;
      } else {
        s1 = (cos(za_sca_rad) - cos(za_inc_rad) * cos(theta_rad)) /
             (sin(za_inc_rad) * sin(theta_rad));
        s2 = (cos(za_inc_rad) - cos(za_sca_rad) * cos(theta_rad)) /
             (sin(za_sca_rad) * sin(theta_rad));

        sigma1 = acos(s1);
        sigma2 = acos(s2);

        // Arccos is only defined in the range from -1 ... 1
        // Numerical problems can appear for values close to 1 or -1
        // this (also) catches the case when inc and sca are on one meridian
        if (std::isnan(sigma1) || std::isnan(sigma2)) {
          if (abs(s1 - 1) < ANGTOL_RAD) sigma1 = 0;
          if (abs(s1 + 1) < ANGTOL_RAD) sigma1 = PI;
          if (abs(s2 - 1) < ANGTOL_RAD) sigma2 = 0;
          if (abs(s2 + 1) < ANGTOL_RAD) sigma2 = PI;
        }
      }

      const Numeric C1 = cos(2 * sigma1);
      const Numeric C2 = cos(2 * sigma2);

      const Numeric S1 = sin(2 * sigma1);
      const Numeric S2 = sin(2 * sigma2);

      pha_mat_lab(0, 1) = C1 * F12;
      pha_mat_lab(1, 0) = C2 * F12;
      pha_mat_lab(1, 1) = C1 * C2 * F22 - S1 * S2 * F33;

      //ARTS_ASSERT(!std::isnan(pha_mat_lab(0,1)));
      //ARTS_ASSERT(!std::isnan(pha_mat_lab(1,0)));
      //ARTS_ASSERT(!std::isnan(pha_mat_lab(1,1)));
      ARTS_USER_ERROR_IF (std::isnan(pha_mat_lab(0, 1)) || std::isnan(pha_mat_lab(1, 0)) ||
          std::isnan(pha_mat_lab(1, 1)),
            "NaN value(s) detected in *pha_mat_labCalc* (0/1,1). Could the "
            "input data contain NaNs? Please check with *scat_dataCheck*. If "
            "input data are OK  and you critically need the ongoing calculations, "
            "try to change the observation LOS slightly. If you can reproduce "
            "this error, please contact Patrick in order to help tracking down "
            "the reason to this problem. If you see this message occasionally "
            "when doing MC calculations, it should not be critical. This path "
            "sampling will be rejected and replaced with a new one.");

      if (stokes_dim > 2) {
        /*CPD: For skokes_dim > 2 some of the transformation formula
            for each element have a different sign depending on whether or
            not 0<aa_scat-aa_inc<180.  For details see pages 94 and 95 of
            Mishchenkos chapter in :
            Mishchenko, M. I., and L. D. Travis, 2003: Electromagnetic
            scattering by nonspherical particles. In Exploring the Atmosphere
            by Remote Sensing Techniques (R. Guzzi, Ed.), Springer-Verlag,
            Berlin, pp. 77-127.
            This is available at http://www.giss.nasa.gov/~crmim/publications/ */
        Numeric delta_aa = aa_sca - aa_inc + (aa_sca - aa_inc < -180) * 360 -
                           (aa_sca - aa_inc > 180) * 360;
        if (delta_aa >= 0) {
          pha_mat_lab(0, 2) = S1 * F12;
          pha_mat_lab(1, 2) = S1 * C2 * F22 + C1 * S2 * F33;
          pha_mat_lab(2, 0) = -S2 * F12;
          pha_mat_lab(2, 1) = -C1 * S2 * F22 - S1 * C2 * F33;
        } else {
          pha_mat_lab(0, 2) = -S1 * F12;
          pha_mat_lab(1, 2) = -S1 * C2 * F22 - C1 * S2 * F33;
          pha_mat_lab(2, 0) = S2 * F12;
          pha_mat_lab(2, 1) = C1 * S2 * F22 + S1 * C2 * F33;
        }
        pha_mat_lab(2, 2) = -S1 * S2 * F22 + C1 * C2 * F33;

        if (stokes_dim > 3) {
          if (delta_aa >= 0) {
            pha_mat_lab(1, 3) = S2 * F34;
            pha_mat_lab(3, 1) = S1 * F34;
          } else {
            pha_mat_lab(1, 3) = -S2 * F34;
            pha_mat_lab(3, 1) = -S1 * F34;
          }
          pha_mat_lab(0, 3) = 0;
          pha_mat_lab(2, 3) = C2 * F34;
          pha_mat_lab(3, 0) = 0;
          pha_mat_lab(3, 2) = -C1 * F34;
          pha_mat_lab(3, 3) = F44;
        }
      }
    }
  }
}

ostream& operator<<(ostream& os, const SingleScatteringData& /*ssd*/) {
  os << "SingleScatteringData: Output operator not implemented";
  return os;
}

ostream& operator<<(ostream& os, const ArrayOfSingleScatteringData& /*assd*/) {
  os << "ArrayOfSingleScatteringData: Output operator not implemented";
  return os;
}

ostream& operator<<(ostream& os, const ScatteringMetaData& /*ssd*/) {
  os << "ScatteringMetaData: Output operator not implemented";
  return os;
}

ostream& operator<<(ostream& os, const ArrayOfScatteringMetaData& /*assd*/) {
  os << "ArrayOfScatteringMetaData: Output operator not implemented";
  return os;
}

//! Get optical properties from propmat_clearsky
/*!
  This turns propmat_clearsky into the extinction matrix
  and absorption vector for use when these are important.

  Internal function to replace the old opt_prop_gas_agenda.

  Output and Input:
  \param ext_mat Extinction matrix.
  \param abs_vec Absorption vector.
  Input:
  \param propmat_clearsky as the WSV.

  \author Richard Larsson
  \date   2012-07-24
*/
void opt_prop_sum_propmat_clearsky(  //Output:
    PropagationMatrix& ext_mat,
    StokesVector& abs_vec,
    //Input:
    const PropagationMatrix& propmat_clearsky) {
   const Index stokes_dim = propmat_clearsky.StokesDimensions();
   const Index freq_dim = propmat_clearsky.NumberOfFrequencies();

  // old abs_vecInit
  abs_vec = StokesVector(freq_dim, stokes_dim);
  abs_vec.SetZero();

  ext_mat = PropagationMatrix(freq_dim, stokes_dim);
  ext_mat.SetZero();

  // old ext_matAddGas and abs_vecAddGas for 0 vector and matrix
  abs_vec += propmat_clearsky;
  ext_mat += propmat_clearsky;
}

//! Convert ptype name to enum value
/*!
 Returns the PType enum value for the given String.

 This is the conversion for SingleScatteringData version 2.

 \param[in]  ptype_string  Particle type name
 \return     PType enum value

 \author Oliver Lemke
 */
PType PTypeFromString(const String& ptype_string) {
  PType ptype;
  if (ptype_string == "general")
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "totally_random")
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "azimuthally_random")
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR (
      "Unknown ptype: ", ptype_string, "\n"
      "Valid types are: general, totally_random and "
      "azimuthally_random.")
  }

  return ptype;
}

//! Convert ptype name to enum value
/*!
 Returns the PType enum value for the given String.

 This is the conversion for SingleScatteringData version 2.

 \param[in]  ptype_string  Particle type name
 \return     PType enum value

 \author Oliver Lemke
 */
PType PType2FromString(const String& ptype_string) {
  PType ptype;
  if (ptype_string == "general")
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "macroscopically_isotropic")
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "horizontally_aligned")
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR (
      "Unknown ptype: ", ptype_string, "\n"
       "Valid types are: general, macroscopically_isotropic and "
       "horizontally_aligned.")
  }

  return ptype;
}

//! Convert particle type enum value to String.
/*!
 Returns the PType enum value for the given String.

 \param[in]  ptype  Particle type
 \return     String representation of PType

 \author Oliver Lemke
 */
String PTypeToString(const PType& ptype) {
  String ptype_string;

  switch (ptype) {
    case PTYPE_GENERAL:
      ptype_string = "general";
      break;
    case PTYPE_TOTAL_RND:
      ptype_string = "totally_random";
      break;
    case PTYPE_AZIMUTH_RND:
      ptype_string = "azimuthally_random";
      break;
    default:
      ARTS_USER_ERROR (
                          "Internal error: Cannot map PType enum value ",
                          ptype, " to String.")
      break;
  }

  return ptype_string;
}

//! Convert azimuthally-random oriented SingleScatteringData to latest version.
/*!
 Converts SingleScatteringData to version 3.

 \param[in,out]  ssd  SingleScatteringData

 \author Jana Mendrok
*/
void ConvertAzimuthallyRandomSingleScatteringData(SingleScatteringData& ssd) {
  // First check that input fulfills requirements on older data formats:
  // 1) Is za_grid symmetric and includes 90deg?
  Index nza = ssd.za_grid.nelem();
  for (Index i = 0; i < nza / 2; i++) {
    ARTS_USER_ERROR_IF (!is_same_within_epsilon(
            180. - ssd.za_grid[nza - 1 - i], ssd.za_grid[i], 2 * DBL_EPSILON),
        "Zenith grid of azimuthally_random single scattering data\n"
        "is not symmetric with respect to 90degree.")
  }
  ARTS_USER_ERROR_IF (!is_same_within_epsilon(ssd.za_grid[nza / 2], 90., 2 * DBL_EPSILON),
      "Zenith grid of azimuthally_random single scattering data\n"
      "does not contain 90 degree grid point.")

  // 2) Are data sizes correct?
  ostringstream os_pha_mat;
  os_pha_mat << "pha_mat ";
  ostringstream os_ext_mat;
  os_ext_mat << "ext_mat ";
  ostringstream os_abs_vec;
  os_abs_vec << "abs_vec ";
  chk_size(os_pha_mat.str(),
           ssd.pha_mat_data,
           ssd.f_grid.nelem(),
           ssd.T_grid.nelem(),
           ssd.za_grid.nelem(),
           ssd.aa_grid.nelem(),
           ssd.za_grid.nelem() / 2 + 1,
           1,
           16);

  chk_size(os_ext_mat.str(),
           ssd.ext_mat_data,
           ssd.f_grid.nelem(),
           ssd.T_grid.nelem(),
           ssd.za_grid.nelem() / 2 + 1,
           1,
           3);

  chk_size(os_abs_vec.str(),
           ssd.abs_vec_data,
           ssd.f_grid.nelem(),
           ssd.T_grid.nelem(),
           ssd.za_grid.nelem() / 2 + 1,
           1,
           2);

  // Now that we are sure that za_grid is properly symmetric, we just need to
  // copy over the data (ie no interpolation).
  Tensor5 tmpT5 = ssd.abs_vec_data;
  ssd.abs_vec_data.resize(tmpT5.nshelves(),
                          tmpT5.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT5.nrows(),
                          tmpT5.ncols());
  ssd.abs_vec_data(joker, joker, Range(0, nza / 2 + 1), joker, joker) = tmpT5;
  for (Index i = 0; i < nza / 2; i++) {
    ssd.abs_vec_data(joker, joker, nza - 1 - i, joker, joker) =
        tmpT5(joker, joker, i, joker, joker);
  }

  tmpT5 = ssd.ext_mat_data;
  ssd.ext_mat_data.resize(tmpT5.nshelves(),
                          tmpT5.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT5.nrows(),
                          tmpT5.ncols());
  ssd.ext_mat_data(joker, joker, Range(0, nza / 2 + 1), joker, joker) = tmpT5;
  for (Index i = 0; i < nza / 2; i++) {
    ssd.ext_mat_data(joker, joker, nza - 1 - i, joker, joker) =
        tmpT5(joker, joker, i, joker, joker);
  }

  Tensor7 tmpT7 = ssd.pha_mat_data;
  ssd.pha_mat_data.resize(tmpT7.nlibraries(),
                          tmpT7.nvitrines(),
                          tmpT7.nshelves(),
                          tmpT7.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT7.nrows(),
                          tmpT7.ncols());
  ssd.pha_mat_data(
      joker, joker, joker, joker, Range(0, nza / 2 + 1), joker, joker) = tmpT7;

  // scatt. matrix elements 13,23,31,32 and 14,24,41,42 (=elements 2,6,8,9 and
  // 3,7,12,13 in ARTS' flattened format, respectively) change sign.
  tmpT7(joker, joker, joker, joker, joker, joker, Range(2, 2)) *= -1.;
  tmpT7(joker, joker, joker, joker, joker, joker, Range(6, 4)) *= -1.;
  tmpT7(joker, joker, joker, joker, joker, joker, Range(12, 2)) *= -1.;

  // For second half of incident polar angles (>90deg), we need to mirror the
  // original data in both incident and scattered polar angle around 90deg "planes".
  for (Index i = 0; i < nza / 2; i++)
    for (Index j = 0; j < nza; j++)
      ssd.pha_mat_data(
          joker, joker, nza - 1 - j, joker, nza - 1 - i, joker, joker) =
          tmpT7(joker, joker, j, joker, i, joker, joker);
}

//! Convert particle ssd method name to enum value
/*!
 Returns the ParticleSSDMethod enum value for the given String.

 \param[in]  particle_ssdmethod_string  Particle SSD method name
 \return     ParticleSSDMethod enum value

 \author Oliver Lemke
 */
ParticleSSDMethod ParticleSSDMethodFromString(
    const String& particle_ssdmethod_string) {
  ParticleSSDMethod particle_ssdmethod;
  if (particle_ssdmethod_string == "tmatrix")
    particle_ssdmethod = PARTICLE_SSDMETHOD_TMATRIX;
  else {
    ARTS_USER_ERROR (
                        "Unknown particle SSD method: ",
                        particle_ssdmethod_string, "\n"
                        "Valid methods: tmatrix")
  }

  return particle_ssdmethod;
}

//! Convert particle type enum value to String.
/*!
 Returns the PType enum value for the given String.

 \param[in]  ptype  Particle type
 \return     String representation of ParticleSSDMethod

 \author Oliver Lemke
 */
String PTypeToString(const ParticleSSDMethod& particle_ssdmethod) {
  String particle_ssdmethod_string;

  switch (particle_ssdmethod) {
    case PARTICLE_SSDMETHOD_TMATRIX:
      particle_ssdmethod_string = "tmatrix";
      break;
    default:
      ARTS_USER_ERROR (
        "Internal error: Cannot map ParticleSSDMethod enum value ",
        particle_ssdmethod, " to String.")
      break;
  }

  return particle_ssdmethod_string;
}


//! Extinction, absorption and phase function for one scattering species, 1D and TRO
/*!
 This function interpolates scat_data and sums up values to obtain the
 extinction, absorption and phase function for one scattering species, on the
 condition that scat_data only contains data of TRO-type and a 1D atmosphere is
 given.

 For flexibility, the quantities are added to ext, mat and pfun. This means
 that you need to set these variables to zero before calling this function for
 the first or each scattering species.

 The function gives error if not all data are TRO.

 ext, abs, pfun and T_grid shall match the complete atmosphere, while pnd shall
 match the cloudbox.

 \param[in,out]  ext    Extinction [frequency, temperature]
 \param[in,out]  abs    Absorption [frequency, temperature]
 \param[in,out]  pfun   Phase function [frequency, temperature, scattering angle]
 \param[in]  scat_data  Scattering data for one scattering species
 \param[in]  pnd      Particle number density [scattering element, cloudbox level]
 \param[in]  T_grid   Temperatures
 \param[in]  sa_grid  Scattering angles
 \param[in]  cloudbox_limits  Limits in pressure dimension


 \author Patrick Eriksspn
*/
void ext_abs_pfun_from_tro(MatrixView ext,
                           MatrixView abs,
                           Tensor3View pfun,
                           const ArrayOfSingleScatteringData& scat_data,
                           ConstMatrixView pnd_data,
                           ArrayOfIndex& cloudbox_limits,
                           ConstVectorView T_grid,
                           ConstVectorView sa_grid)
{
  // Sizes
  const Index nse = scat_data.nelem();
  const Index nf = scat_data[0].f_grid.nelem();
  const Index nt = T_grid.nelem();
  const Index nsa = sa_grid.nelem();
  const Index ncl = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  

  ARTS_ASSERT(ext.nrows() == nf);
  ARTS_ASSERT(ext.ncols() == nt);
  ARTS_ASSERT(abs.nrows() == nf);
  ARTS_ASSERT(abs.ncols() == nt);
  ARTS_ASSERT(pfun.npages() == nf);
  ARTS_ASSERT(pfun.nrows() == nt);
  ARTS_ASSERT(pfun.ncols() == nsa);
  ARTS_ASSERT(cloudbox_limits.nelem() == 2);
  ARTS_ASSERT(pnd_data.nrows() == nse);
  ARTS_ASSERT(pnd_data.ncols() == ncl);

  // Check that all data is TRO
  {
    bool all_totrand = true;
    for (Index ie = 0; ie < nse; ie++) {
      if (scat_data[ie].ptype != PTYPE_TOTAL_RND)
        all_totrand = false;
    }
    ARTS_USER_ERROR_IF (!all_totrand,
                        "This method demands that all scat_data are TRO");
  }

  // Range for layers inside the cloudbox
  const Index i0 = cloudbox_limits[0];
  Range clevels = Range(i0, ncl);
  
  // Loop scattering elements
  for (Index ie = 0; ie < nse; ie++) {
    // Add any check that pnd_data(ie,:) > 0?

    // Temperature-only interpolation weights
    ArrayOfGridPos gp_t(ncl);
    gridpos(gp_t, scat_data[ie].T_grid, T_grid[clevels]);
    Matrix itw1(ncl, 2);
    interpweights(itw1, gp_t);

    // Temperature + scattering angle interpolation weights
    ArrayOfGridPos gp_sa(nsa);
    gridpos(gp_sa, scat_data[ie].za_grid, sa_grid);
    Tensor3 itw2(ncl, nsa, 4);
    interpweights(itw2, gp_t, gp_sa);

    // Loop frequencies
    for (Index iv = 0; iv < nf; iv++) {
      // Interpolate
      Vector ext1(ncl), abs1(ncl);
      Matrix pfu1(ncl,nsa);
      interp(ext1, itw1, scat_data[ie].ext_mat_data(iv,joker,0,0,0), gp_t);
      interp(abs1, itw1, scat_data[ie].abs_vec_data(iv,joker,0,0,0), gp_t);
      interp(pfu1, itw2, scat_data[ie].pha_mat_data(iv,joker,joker,0,0,0,0),
             gp_t, gp_sa);

      // Add to container variables
      for (Index il = 0; il < ncl; il++) {
        ext(iv,i0+il) += pnd_data(ie,il) * ext1[il];
        abs(iv,i0+il) += pnd_data(ie,il) * abs1[il];
        for (Index ia = 0; ia < nsa; ia++) {
          pfun(iv,i0+il,ia) += pnd_data(ie,il) * pfu1(il,ia);
        }
      }
    }
  }
}
