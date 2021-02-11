/* Copyright (C) 2006-2012 Claudia Emde <claudia.emde@dlr.de>
                      
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
   USA. 
*/

/**
 * @file   disort.cc
 * @author Claudia Emde <claudia.emde@dlr.de>
 * @date   Tue Feb  7 10:08:28 2006
 * 
 * @brief  This file contains functions related to the DISORT interface.
 */

#include "disort.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"

extern "C" {
#include "cdisort.h"
}

#include "disort.h"
#include "interpolation.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "rte.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric COSMIC_BG_TEMP;

void check_disort_input(  // Input
    const Index& cloudbox_on,
    const Index& atmfields_checked,
    const Index& atmgeom_checked,
    const Index& cloudbox_checked,
    const Index& scat_data_checked,
    const Index& atmosphere_dim,
    const Index& stokes_dim,
    const ArrayOfIndex& cloudbox_limits,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    ConstVectorView za_grid,
    const Index& nstreams) {
  if (!cloudbox_on) {
    throw runtime_error(
        "Cloudbox is off, no scattering calculations to be"
        "performed.");
  }

  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  if (atmgeom_checked != 1)
    throw runtime_error(
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  if (cloudbox_checked != 1)
    throw runtime_error(
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  if (scat_data_checked != 1)
    throw runtime_error(
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  if (atmosphere_dim != 1)
    throw runtime_error(
        "For running DISORT, atmospheric dimensionality "
        "must be 1.\n");

  if (stokes_dim < 0 || stokes_dim > 1)
    throw runtime_error(
        "For running DISORT, the dimension of stokes vector "
        "must be 1.\n");

  if (cloudbox_limits.nelem() != 2 * atmosphere_dim)
    throw runtime_error(
        "*cloudbox_limits* is a vector which contains the"
        "upper and lower limit of the cloud for all "
        "atmospheric dimensions. So its dimension must"
        "be 2 x *atmosphere_dim*");

  if (cloudbox_limits[0] != 0) {
    ostringstream os;
    os << "DISORT calculations currently only possible with "
       << "lower cloudbox limit\n"
       << "at 0th atmospheric level "
       << "(assumes surface there, ignoring z_surface).\n";
    throw runtime_error(os.str());
  }

  if (scat_data.empty())
    throw runtime_error(
        "No single scattering data present.\n"
        "See documentation of WSV *scat_data* for options to "
        "make single scattering data available.\n");

  // DISORT requires even number of streams:
  // nstreams is total number of directions, up- and downwelling, and the up-
  // and downwelling directions need to be symmetrically distributed, i.e. same
  // number of directions in both hemispheres is required. horizontal direction
  // (90deg) can not be covered in a plane-parallel atmosphere.
  if (nstreams / 2 * 2 != nstreams) {
    ostringstream os;
    os << "DISORT requires an even number of streams, but yours is " << nstreams
       << ".\n";
    throw runtime_error(os.str());
  }

  // Zenith angle grid.
  Index nza = za_grid.nelem();

  // za_grid here is only relevant to provide an i_field from which the
  // sensor los angles can be interpolated by yCalc; it does not the determine
  // the accuracy of the DISORT output itself at these angles. So we can only
  // apply a very rough test here, whether the grid is appropriate. However, we
  // set the threshold fairly high since calculation costs for a higher number
  // of angles are negligible.
  if (nza < 20) {
    ostringstream os;
    os << "We require size of za_grid to be >= 20, to ensure a\n"
       << "reasonable interpolation of the calculated cloudbox field.\n"
       << "Note that for DISORT additional computation costs for\n"
       << "larger numbers of angles are negligible.";
    throw runtime_error(os.str());
  }

  if (za_grid[0] != 0. || za_grid[nza - 1] != 180.)
    throw runtime_error("The range of *za_grid* must [0 180].");

  if (!is_increasing(za_grid))
    throw runtime_error("*za_grid* must be increasing.");

  Index i = 1;
  while (za_grid[i] <= 90) {
    if (za_grid[i] == 90)
      throw runtime_error("*za_grid* is not allowed to contain the value 90");
    i++;
  }

  // DISORT can only handle randomly oriented particles.
  bool all_totrand = true;
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++)
      if (scat_data[i_ss][i_se].ptype != PTYPE_TOTAL_RND) all_totrand = false;
  if (!all_totrand) {
    ostringstream os;
    os << "DISORT can only handle scattering elements of type "
       << PTYPE_TOTAL_RND << " (" << PTypeToString(PTYPE_TOTAL_RND) << "),\n"
       << "but at least one element of other type (" << PTYPE_AZIMUTH_RND << "="
       << PTypeToString(PTYPE_AZIMUTH_RND) << " or " << PTYPE_GENERAL << "="
       << PTypeToString(PTYPE_GENERAL) << ") is present.\n";
    throw runtime_error(os.str());
  }
}

void init_ifield(  // Output
    Tensor7& cloudbox_field,
    // Input
    const Vector& f_grid,
    const ArrayOfIndex& cloudbox_limits,
    const Index& nang,
    const Index& stokes_dim) {
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  //const Index Nza = za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  cloudbox_field.resize(Nf, Np_cloud, 1, 1, nang, 1, stokes_dim);
  cloudbox_field = NAN;
}

void get_disortsurf_props(  // Output
    Vector& albedo,
    Numeric& btemp,
    // Input
    ConstVectorView f_grid,
    const Numeric& surface_skin_t,
    ConstVectorView surface_scalar_reflectivity) {
  // temperature of surface
  if (surface_skin_t < 0. || surface_skin_t > 1000.) {
    ostringstream os;
    os << "Surface temperature has been set or derived as " << btemp << " K,\n"
       << "which is not considered a meaningful value.\n"
       << "For surface method 'L', *surface_skin_t* needs to\n"
       << "be set and passed explicitly. Maybe you didn't do this?";
    throw runtime_error(os.str());
  }
  btemp = surface_skin_t;

  // surface albedo
  if (surface_scalar_reflectivity.nelem() != f_grid.nelem() &&
      surface_scalar_reflectivity.nelem() != 1) {
    ostringstream os;
    os << "The number of elements in *surface_scalar_reflectivity*\n"
       << "should match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << f_grid.nelem()
       << "\n length of *surface_scalar_reflectivity* : "
       << surface_scalar_reflectivity.nelem() << "\n";
    throw runtime_error(os.str());
  }

  if (min(surface_scalar_reflectivity) < 0 ||
      max(surface_scalar_reflectivity) > 1) {
    throw runtime_error(
        "All values in *surface_scalar_reflectivity*"
        " must be inside [0,1].");
  }

  if (surface_scalar_reflectivity.nelem() > 1)
    for (Index f_index = 0; f_index < f_grid.nelem(); f_index++)
      albedo[f_index] = surface_scalar_reflectivity[f_index];
  else
    for (Index f_index = 0; f_index < f_grid.nelem(); f_index++)
      albedo[f_index] = surface_scalar_reflectivity[0];
}

void get_gasoptprop(Workspace& ws,
                    MatrixView ext_bulk_gas,
                    const Agenda& propmat_clearsky_agenda,
                    ConstVectorView t_profile,
                    ConstMatrixView vmr_profiles,
                    ConstVectorView p_grid,
                    ConstVectorView f_grid) {
  const Index Np = p_grid.nelem();

  ARTS_ASSERT(ext_bulk_gas.nrows() == f_grid.nelem());
  ARTS_ASSERT(ext_bulk_gas.ncols() == Np);

  // Initialization
  ext_bulk_gas = 0.;

  // making gas property output containers and input dummies
  const EnergyLevelMap rtp_nlte_dummy;
  const Vector rtp_mag_dummy(3, 0);
  const Vector ppath_los_dummy;
  StokesVector nlte_dummy;
  ArrayOfStokesVector partial_nlte_dummy;
  ArrayOfPropagationMatrix partial_dummy;

  PropagationMatrix propmat_clearsky_local;
  for (Index ip = 0; ip < Np; ip++) {
    propmat_clearsky_agendaExecute(ws,
                                   propmat_clearsky_local,
                                   nlte_dummy,
                                   partial_dummy,
                                   partial_nlte_dummy,
                                   ArrayOfRetrievalQuantity(0),
                                   f_grid,
                                   rtp_mag_dummy,
                                   ppath_los_dummy,
                                   p_grid[ip],
                                   t_profile[ip],
                                   rtp_nlte_dummy,
                                   vmr_profiles(joker, ip),
                                   propmat_clearsky_agenda);
    ext_bulk_gas(joker, ip) += propmat_clearsky_local.Kjj();
  }
}

void get_paroptprop(MatrixView ext_bulk_par,
                    MatrixView abs_bulk_par,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    ConstMatrixView pnd_profiles,
                    ConstVectorView t_profile,
                    ConstVectorView DEBUG_ONLY(p_grid),
                    const ArrayOfIndex& cloudbox_limits,
                    ConstVectorView f_grid) {
  const Index Np_cloud = pnd_profiles.ncols();
  DEBUG_ONLY(const Index Np = p_grid.nelem());
  const Index nf = f_grid.nelem();

  ARTS_ASSERT(ext_bulk_par.nrows() == nf);
  ARTS_ASSERT(abs_bulk_par.nrows() == nf);
  ARTS_ASSERT(ext_bulk_par.ncols() == Np);
  ARTS_ASSERT(abs_bulk_par.ncols() == Np);

  // Initialization
  ext_bulk_par = 0.;
  abs_bulk_par = 0.;

  // preparing input data
  Vector T_array = t_profile[Range(cloudbox_limits[0], Np_cloud)];
  Matrix dir_array(1, 2, 0.);  // just a dummy. only tot_random allowed, ie.
  // optprop are independent of direction.

  // making particle property output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk;  //nf,nT,ndir,nst,nst
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  // calculate particle optical properties
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      1,
                      T_array,
                      dir_array,
                      -1);
  opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                        abs_vec_ssbulk,
                        ptype_ssbulk,
                        ext_mat_Nse,
                        abs_vec_Nse,
                        ptypes_Nse,
                        pnd_profiles,
                        t_ok);
  opt_prop_Bulk(ext_mat_bulk,
                abs_vec_bulk,
                ptype_bulk,
                ext_mat_ssbulk,
                abs_vec_ssbulk,
                ptype_ssbulk);

  Index f_this = 0;
  bool pf = (abs_vec_bulk.nbooks() != 1);
  for (Index ip = 0; ip < Np_cloud; ip++)
    for (Index f_index = 0; f_index < nf; f_index++) {
      if (pf) f_this = f_index;
      ext_bulk_par(f_index, ip + cloudbox_limits[0]) =
          ext_mat_bulk(f_this, ip, 0, 0, 0);
      abs_bulk_par(f_index, ip + cloudbox_limits[0]) =
          abs_vec_bulk(f_this, ip, 0, 0);
    }
}

void get_dtauc_ssalb(MatrixView dtauc,
                     MatrixView ssalb,
                     ConstMatrixView ext_bulk_gas,
                     ConstMatrixView ext_bulk_par,
                     ConstMatrixView abs_bulk_par,
                     ConstVectorView z_profile) {
  const Index nf = ext_bulk_gas.nrows();
  const Index Np = ext_bulk_gas.ncols();

  ARTS_ASSERT(dtauc.nrows() == nf);
  ARTS_ASSERT(ssalb.nrows() == nf);
  ARTS_ASSERT(dtauc.ncols() == Np - 1);
  ARTS_ASSERT(ssalb.ncols() == Np - 1);

  // Initialization
  dtauc = 0.;
  ssalb = 0.;

  for (Index ip = 0; ip < Np - 1; ip++)
    // Do layer averaging and derive single scattering albedo & optical depth
    for (Index f_index = 0; f_index < nf; f_index++) {
      Numeric ext =
          0.5 * (ext_bulk_gas(f_index, ip) + ext_bulk_par(f_index, ip) +
                 ext_bulk_gas(f_index, ip + 1) + ext_bulk_par(f_index, ip + 1));
      if (ext != 0) {
        Numeric abs =
            0.5 *
            (ext_bulk_gas(f_index, ip) + abs_bulk_par(f_index, ip) +
             ext_bulk_gas(f_index, ip + 1) + abs_bulk_par(f_index, ip + 1));
        ssalb(f_index, Np - 2 - ip) = (ext - abs) / ext;
      }

      dtauc(f_index, Np - 2 - ip) = ext * (z_profile[ip + 1] - z_profile[ip]);
    }
}

void get_angs(Vector& pfct_angs,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Index& Npfct) {
  const Index min_nang = 3;
  Index nang = Npfct;

  if (Npfct < 0) {
    Index this_ss = 0, this_se = 0;
    // determine nang and pfct_angs from scat_data with finest za_grid
    for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
      for (Index i_se = scat_data[i_ss].nelem() - 1; i_se >= 0; i_se--)
        // considering scat elems within one species mostly sorted from small to
        // large sizes with large sizes corresponding to large za_grid. that is,
        // starting searching from back should trigger if-clause (and variable
        // update) less often.
        if (nang < scat_data[i_ss][i_se].za_grid.nelem()) {
          nang = scat_data[i_ss][i_se].za_grid.nelem();
          this_ss = i_ss;
          this_se = i_se;
        }
    pfct_angs = scat_data[this_ss][this_se].za_grid;
  } else if (Npfct < min_nang) {
    ostringstream os;
    os << "Number of requested angular grid points (Npfct=" << Npfct
       << ") is insufficient.\n"
       << "At least " << min_nang << " points required.\n";
    throw runtime_error(os.str());
  } else {
    nlinspace(pfct_angs, 0, 180, nang);
  }
}

void get_parZ(Tensor3& pha_bulk_par,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              ConstMatrixView pnd_profiles,
              ConstVectorView t_profile,
              ConstVectorView pfct_angs,
              const ArrayOfIndex& cloudbox_limits) {
  const Index Np_cloud = pnd_profiles.ncols();
  const Index nang = pfct_angs.nelem();

  // Initialization
  pha_bulk_par = 0.;

  // preparing input data
  Vector T_array = t_profile[Range(cloudbox_limits[0], Np_cloud)];
  Matrix idir_array(1, 2, 0.);  // we want pfct on sca ang grid, hence set
  // pdir(*,0) to sca ang, all other to 0.
  Matrix pdir_array(nang, 2, 0.);
  pdir_array(joker, 0) = pfct_angs;

  // making particle property output containers
  ArrayOfArrayOfTensor6 pha_mat_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 pha_mat_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 pha_mat_bulk;  //nf,nT,npdir,nidir,nst,nst
  Index ptype_bulk;

  // calculate phase matrix
  // FIXME: might be optimized by instead just executing pha_mat_1ScatElem where
  // ext_bulk_par or pnd_field are non-zero.
  pha_mat_NScatElems(pha_mat_Nse,
                     ptypes_Nse,
                     t_ok,
                     scat_data,
                     1,
                     T_array,
                     pdir_array,
                     idir_array,
                     -1);
  pha_mat_ScatSpecBulk(pha_mat_ssbulk,
                       ptype_ssbulk,
                       pha_mat_Nse,
                       ptypes_Nse,
                       pnd_profiles,
                       t_ok);
  pha_mat_Bulk(pha_mat_bulk, ptype_bulk, pha_mat_ssbulk, ptype_ssbulk);

  pha_bulk_par(joker, Range(cloudbox_limits[0], Np_cloud), joker) =
      pha_mat_bulk(joker, joker, joker, 0, 0, 0);
}

void get_pfct(Tensor3& pfct_bulk_par,
              ConstTensor3View& pha_bulk_par,
              ConstMatrixView ext_bulk_par,
              ConstMatrixView abs_bulk_par,
              const ArrayOfIndex& cloudbox_limits) {
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Np = pha_bulk_par.nrows();
  const Index nf = pha_bulk_par.npages();
  Index nang = pha_bulk_par.ncols();

  ARTS_ASSERT(pfct_bulk_par.npages() == nf);
  ARTS_ASSERT(pfct_bulk_par.nrows() == Np - 1);
  ARTS_ASSERT(pfct_bulk_par.ncols() == nang);

  // Initialization
  pfct_bulk_par = 0.;

  for (Index ip = cloudbox_limits[0]; ip < Np_cloud - 1; ip++)
    for (Index f_index = 0; f_index < nf; f_index++) {
      // Calculate layer averaged scattering (omitting factor 0.5 here as
      // omitting it for layer averaged Z below, too).
      Numeric sca =
          (ext_bulk_par(f_index, ip) + ext_bulk_par(f_index, ip + 1)) -
          (abs_bulk_par(f_index, ip) + abs_bulk_par(f_index, ip + 1));
      if (sca != 0.) {
        // Calculate layer averaged Z (omitting factor 0.5) and rescale from
        // Z (Csca) to P (4Pi)
        for (Index ia = 0; ia < nang; ia++)
          pfct_bulk_par(f_index, Np - 2 - ip, ia) +=
              pha_bulk_par(f_index, ip, ia) + pha_bulk_par(f_index, ip + 1, ia);
        pfct_bulk_par(f_index, Np - 2 - ip, joker) *= 4 * PI / sca;
      }
    }
}

void get_pmom(Tensor3View pmom,
              ConstTensor3View pfct_bulk_par,
              ConstVectorView pfct_angs,
              const Index& Nlegendre) {
  const Index nf = pfct_bulk_par.npages();
  const Index nlyr = pfct_bulk_par.nrows();
  const Index nang = pfct_bulk_par.ncols();

  ARTS_ASSERT(nang == pfct_angs.nelem());

  ARTS_ASSERT(pmom.npages() == nf);
  ARTS_ASSERT(pmom.nrows() == nlyr);
  ARTS_ASSERT(pmom.ncols() == Nlegendre);

  Numeric pfct_threshold = 0.1;

  // Initialization
  pmom = 0.;

  // we need the cosine of the pfct angles
  Vector u(nang), adu(nang - 1);
  Tensor3 px(nang - 1, Nlegendre, 2, 0.);
  u[0] = cos(pfct_angs[0] * PI / 180.);
  px(joker, 0, joker) = 1.;
  for (Index ia = 1; ia < nang; ia++) {
    u[ia] = cos(pfct_angs[ia] * PI / 180.);
    adu[ia - 1] = abs(u[ia] - u[ia - 1]);
    px(ia - 1, 1, 0) = u[ia - 1];
    px(ia - 1, 1, 1) = u[ia];
    for (Index l = 2; l < Nlegendre; l++) {
      Numeric dl = (double)l;
      px(ia - 1, l, 0) = (2 * dl - 1) / dl * u[ia - 1] * px(ia - 1, l - 1, 0) -
                         (dl - 1) / dl * px(ia - 1, l - 2, 0);
      px(ia - 1, l, 1) = (2 * dl - 1) / dl * u[ia] * px(ia - 1, l - 1, 1) -
                         (dl - 1) / dl * px(ia - 1, l - 2, 1);
    }
  }

  for (Index il = 0; il < nlyr; il++)
    if (pfct_bulk_par(joker, il, 0).sum() != 0.)
      for (Index f_index = 0; f_index < nf; f_index++) {
        if (pfct_bulk_par(f_index, il, 0) != 0) {
          Vector pfct = pfct_bulk_par(f_index, il, joker);

          // Check if phase function is properly normalized
          Numeric pint = 0.;
          for (Index ia = 0; ia < nang - 1; ia++)
            pint += 0.5 * adu[ia] * (pfct[ia] + pfct[ia + 1]);

          if (abs(pint / 2. - 1.) > pfct_threshold) {
            ostringstream os;
            os << "Phase function normalization deviates from expected value by\n"
               << 1e2 * pint / 2. - 1e2 << "(allowed: " << pfct_threshold * 1e2
               << "%).\n"
               << "Occurs at layer #" << il << " and frequency #" << f_index
               << ".\n"
               << "Something is wrong with your scattering data. Check!\n";
            throw runtime_error(os.str());
          }

          // for the rest, rescale pfct to norm 2
          pfct *= 2. / pint;

          pmom(f_index, il, 0) = 1.;
          for (Index ia = 0; ia < nang - 1; ia++) {
            //for (Index l=0; l<Nlegendre; l++)
            for (Index l = 1; l < Nlegendre; l++)
              pmom(f_index, il, l) +=
                  0.25 * adu[ia] *
                  (px(ia, l, 0) * pfct[ia] + px(ia, l, 1) * pfct[ia + 1]);
          }
        }
      }
}

// Use a thread_local variable to communicate the Verbosity to the
// Disort error and warning functions. Ugly workaround, to avoid
// passing a Verbosity argument throughout the whole cdisort code.
// We want to avoid changes to the original code to keep it maintainable
// in respect to upstream updates.
thread_local Verbosity disort_verbosity;

#define MAX_WARNINGS 100

/** Verbosity enabled replacement for the original cdisort function. */
void c_errmsg(const char* messag, int type) {
  Verbosity verbosity = disort_verbosity;
  static int warning_limit = FALSE, num_warnings = 0;

  if (type == DS_ERROR) {
    CREATE_OUT0;
    out0 << "  ******* ERROR >>>>>>  " << messag << "\n";
    arts_exit(1);
  }

  if (warning_limit) return;

  if (++num_warnings <= MAX_WARNINGS) {
    CREATE_OUT1;
    out1 << "  ******* WARNING >>>>>>  " << messag << "\n";
  } else {
    CREATE_OUT1;
    out1 << "  >>>>>>  TOO MANY WARNING MESSAGES --  They will no longer be "
            "printed  <<<<<<<\n\n";
    warning_limit = TRUE;
  }

  return;
}

#undef MAX_WARNINGS

/** Verbosity enabled replacement for the original cdisort function. */
int c_write_bad_var(int quiet, const char* varnam) {
  const int maxmsg = 50;
  static int nummsg = 0;

  nummsg++;
  if (quiet != QUIET) {
    Verbosity verbosity = disort_verbosity;
    CREATE_OUT1;
    out1 << "  ****  Input variable " << varnam << " in error  ****\n";
    if (nummsg == maxmsg) {
      c_errmsg("Too many input errors.  Aborting...", DS_ERROR);
    }
  }

  return TRUE;
}

/** Verbosity enabled replacement for the original cdisort function. */
int c_write_too_small_dim(int quiet, const char* dimnam, int minval) {
  if (quiet != QUIET) {
    Verbosity verbosity = disort_verbosity;
    CREATE_OUT1;
    out1 << "  ****  Symbolic dimension " << dimnam
         << " should be increased to at least " << minval << "  ****\n";
  }

  return TRUE;
}

void reduced_1datm(Vector& p,
                   Vector& z,
                   Vector& t,
                   Matrix& vmr,
                   Matrix& pnd,
                   ArrayOfIndex& cboxlims,
                   Index& ncboxremoved,
                   ConstVectorView p_grid,
                   ConstVectorView z_profile,
                   const Numeric& z_surface,
                   ConstVectorView t_profile,
                   ConstMatrixView vmr_profiles,
                   ConstMatrixView pnd_profiles,
                   const ArrayOfIndex& cloudbox_limits) {
  // Surface at p_grid[0] and we just need to copy the original data
  if (abs(z_surface - z_profile[0]) < 1e-3) {
    p = p_grid;
    z = z_profile;
    t = t_profile;
    vmr = vmr_profiles;
    pnd = pnd_profiles;
    cboxlims = cloudbox_limits;
    ncboxremoved = 0;
  }
  // Surface above p_grid[0]
  else {
    // Some counters
    Index np = p_grid.nelem(), ifirst = 0;
    // Determine where to start with respect to z_profile
    for (; z_surface >= z_profile[ifirst + 1]; ifirst++) {
    }
    np -= ifirst;
    // Start by copying from ifirst to end
    Range ind(ifirst, np);
    p = p_grid[ind];
    z = z_profile[ind];
    t = t_profile[ind];
    vmr = vmr_profiles(joker, ind);
    // Insert surface altitude
    z[0] = z_surface;
    // Prepare interpolation
    ArrayOfGridPos gp(1);
    gridpos(gp[0], z_profile, z_surface);
    Vector itw(2);
    interpweights(itw, gp[0]);
    // t and vmr
    t[0] = interp(itw, t, gp[0]);
    for (int i = 0; i < vmr.nrows(); i++) {
      vmr(i, 0) = interp(itw, vmr(i, joker), gp[0]);
    }
    // p (we need a matrix version of iwt to use the function *itw2p*)
    Matrix itw2(1, 2);
    itw2(0, 0) = itw[0];
    itw2(0, 1) = itw[1];
    itw2p(p[0], p, gp, itw2);
    // pnd_field and cloudbox limits need special treatment
    cboxlims = cloudbox_limits;
    if (ifirst < cloudbox_limits[0]) {  // Surface below cloudbox
      cboxlims[0] -= ifirst;
      cboxlims[1] -= ifirst;
      pnd = pnd_profiles;
      ncboxremoved = 0;
    } else {  // Surface inside cloudbox
      ncboxremoved = ifirst - cboxlims[0];
      cboxlims[0] = 0;
      cboxlims[1] = cloudbox_limits[1] - cloudbox_limits[0] - ncboxremoved;
      ind = Range(ncboxremoved, cboxlims[1] + 1);
      pnd = pnd_profiles(joker, ind);
      gp[0].idx -= cloudbox_limits[0] + ncboxremoved;
      for (int i = 0; i < pnd.nrows(); i++) {
        pnd(i, 0) = interp(itw, pnd(i, joker), gp[0]);
      }
    }
  }
}

void run_cdisort(Workspace& ws,
                 Tensor7& cloudbox_field,
                 ConstVectorView f_grid,
                 ConstVectorView p_grid,
                 ConstVectorView z_profile,
                 const Numeric& z_surface,
                 ConstVectorView t_profile,
                 ConstMatrixView vmr_profiles,
                 ConstMatrixView pnd_profiles,
                 const ArrayOfArrayOfSingleScatteringData& scat_data,
                 const Agenda& propmat_clearsky_agenda,
                 const ArrayOfIndex& cloudbox_limits,
                 const Numeric& surface_skin_t,
                 const Vector& surface_scalar_reflectivity,
                 ConstVectorView za_grid,
                 const Index& nstreams,
                 const Index& Npfct,
                 const Index& quiet,
                 const Verbosity& verbosity) {
  // Create an atmosphere starting at z_surface
  Vector p, z, t;
  Matrix vmr, pnd;
  ArrayOfIndex cboxlims;
  Index ncboxremoved;
  //
  reduced_1datm(p,
                z,
                t,
                vmr,
                pnd,
                cboxlims,
                ncboxremoved,
                p_grid,
                z_profile,
                z_surface,
                t_profile,
                vmr_profiles,
                pnd_profiles,
                cloudbox_limits);

  disort_state ds;
  disort_output out;

  if (quiet == 0)
    disort_verbosity = verbosity;
  else
    disort_verbosity = Verbosity(0, 0, 0);

  const Index nf = f_grid.nelem();

  ds.accur = 0.005;
  ds.flag.prnt[0] = FALSE;
  ds.flag.prnt[1] = FALSE;
  ds.flag.prnt[2] = FALSE;
  ds.flag.prnt[3] = FALSE;
  ds.flag.prnt[4] = TRUE;

  ds.flag.usrtau = FALSE;
  ds.flag.usrang = TRUE;
  ds.flag.spher = FALSE;
  ds.flag.general_source = FALSE;
  ds.flag.output_uum = FALSE;

  ds.nlyr = static_cast<int>(p.nelem() - 1);

  ds.flag.brdf_type = BRDF_NONE;

  ds.flag.ibcnd = GENERAL_BC;
  ds.flag.usrang = TRUE;
  ds.flag.planck = TRUE;
  ds.flag.onlyfl = FALSE;
  ds.flag.lamber = TRUE;
  ds.flag.quiet = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nstr = static_cast<int>(nstreams);
  ds.nphase = ds.nstr;
  ds.nmom = ds.nstr;
  //ds.ntau = ds.nlyr + 1;   // With ds.flag.usrtau = FALSE; set by cdisort
  ds.numu = static_cast<int>(za_grid.nelem());
  ds.nphi = 1;
  Index Nlegendre = nstreams + 1;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &out);

  // Properties of solar beam, set to zero as they are not needed
  ds.bc.fbeam = 0.;
  ds.bc.umu0 = 0.;
  ds.bc.phi0 = 0.;
  ds.bc.fluor = 0.;

  // Since we have no solar source there is no angular dependance
  ds.phi[0] = 0.;

  for (Index i = 0; i <= ds.nlyr; i++) ds.temper[i] = t[ds.nlyr - i];

  Matrix ext_bulk_gas(nf, ds.nlyr + 1);
  get_gasoptprop(ws, ext_bulk_gas, propmat_clearsky_agenda, t, vmr, p, f_grid);
  Matrix ext_bulk_par(nf, ds.nlyr + 1), abs_bulk_par(nf, ds.nlyr + 1);
  get_paroptprop(
      ext_bulk_par, abs_bulk_par, scat_data, pnd, t, p, cboxlims, f_grid);

  // Optical depth of layers
  Matrix dtauc(nf, ds.nlyr);
  // Single scattering albedo of layers
  Matrix ssalb(nf, ds.nlyr);
  get_dtauc_ssalb(dtauc, ssalb, ext_bulk_gas, ext_bulk_par, abs_bulk_par, z);

  // Transform to mu, starting with negative values
  for (Index i = 0; i < ds.numu; i++) ds.umu[i] = -cos(za_grid[i] * PI / 180);

  //upper boundary conditions:
  // DISORT offers isotropic incoming radiance or emissivity-scaled planck
  // emission. Both are applied additively.
  // We want to have cosmic background radiation, for which ttemp=COSMIC_BG_TEMP
  // and temis=1 should give identical results to fisot(COSMIC_BG_TEMP). As they
  // are additive we should use either the one or the other.
  // Note: previous setup (using fisot) setting temis=0 should be avoided.
  // Generally, temis!=1 should be avoided since that technically implies a
  // reflective upper boundary (though it seems that this is not exploited in
  // DISORT1.2, which we so far use).

  // Cosmic background
  // we use temis*ttemp as upper boundary specification, hence CBR set to 0.
  ds.bc.fisot = 0;

  // Top of the atmosphere temperature and emissivity
  ds.bc.ttemp = COSMIC_BG_TEMP;
  ds.bc.btemp = surface_skin_t;
  ds.bc.temis = 1.;

  Vector pfct_angs;
  get_angs(pfct_angs, scat_data, Npfct);
  Index nang = pfct_angs.nelem();

  Index nf_ssd = scat_data[0][0].f_grid.nelem();
  Tensor3 pha_bulk_par(nf_ssd, ds.nlyr + 1, nang);
  get_parZ(pha_bulk_par, scat_data, pnd, t, pfct_angs, cboxlims);
  Tensor3 pfct_bulk_par(nf_ssd, ds.nlyr, nang);
  get_pfct(pfct_bulk_par, pha_bulk_par, ext_bulk_par, abs_bulk_par, cboxlims);

  // Legendre polynomials of phase function
  Tensor3 pmom(nf_ssd, ds.nlyr, Nlegendre, 0.);
  get_pmom(pmom, pfct_bulk_par, pfct_angs, Nlegendre);

  for (Index f_index = 0; f_index < f_grid.nelem(); f_index++) {
    sprintf(ds.header, "ARTS Calc f_index = %ld", f_index);

    std::memcpy(ds.dtauc,
                dtauc(f_index, joker).get_c_array(),
                sizeof(Numeric) * ds.nlyr);
    std::memcpy(ds.ssalb,
                ssalb(f_index, joker).get_c_array(),
                sizeof(Numeric) * ds.nlyr);

    // Wavenumber in [1/cm]
    ds.wvnmhi = ds.wvnmlo = (f_grid[f_index]) / (100. * SPEED_OF_LIGHT);
    ds.wvnmhi += ds.wvnmhi * 1e-7;
    ds.wvnmlo -= ds.wvnmlo * 1e-7;

    ds.bc.albedo = surface_scalar_reflectivity[f_index];

    std::memcpy(ds.pmom,
                pmom(f_index, joker, joker).get_c_array(),
                sizeof(Numeric) * pmom.nrows() * pmom.ncols());

    c_disort(&ds, &out);

    for (Index j = 0; j < ds.numu; j++) {
      for (Index k = cboxlims[1] - cboxlims[0]; k >= 0; k--) {
        cloudbox_field(f_index, k + ncboxremoved, 0, 0, j, 0, 0) =
            out.uu[ds.numu * (ds.nlyr - k - cboxlims[0]) + j] /
            (ds.wvnmhi - ds.wvnmlo) / (100 * SPEED_OF_LIGHT);
      }
      // To avoid potential numerical problems at interpolation of the field,
      // we copy the surface field to underground altitudes
      for (Index k = ncboxremoved - 1; k >= 0; k--) {
        cloudbox_field(f_index, k, 0, 0, j, 0, 0) =
            cloudbox_field(f_index, k + 1, 0, 0, j, 0, 0);
      }
    }
  }

  /* Free allocated memory */
  c_disort_out_free(&ds, &out);
  c_disort_state_free(&ds);
}

void surf_albedoCalc(Workspace& ws,
                     //Output
                     VectorView albedo,
                     Numeric& btemp,
                     //Input
                     const Agenda& surface_rtprop_agenda,
                     ConstVectorView f_grid,
                     ConstVectorView scat_za_grid,
                     const Numeric& surf_alt,
                     const Verbosity& verbosity) {
  // Here, we derive an average surface albedo of the setup as given by ARTS'
  // surface_rtprop_agenda to use with Disorts's proprietary Lambertian surface.
  // In this way, ARTS-Disort can approximately mimick all surface reflection
  // types that ARTS itself can handle.
  // Surface temperature as derived from surface_rtprop_agenda is also returned.
  //
  // We derive the reflection matrices over all incident and reflected polar
  // angle directions and derive their integrated value (or weighted average).
  // The surface_rtprop_agenda handles one reflected direction at a time
  // (rtp_los) and for the reflected directions we loop over all (upwelling)
  // angles as given by scat_za_grid. The derived reflection matrices already
  // represent the reflectivity (or power reflection coefficient) for the given
  // reflected direction including proper angle weighting. For integrating/
  // averaging over the reflected directions, we use the same approach, i.e.
  // weight each angle by its associated range as given by the half distances to
  // the neighboring grid angles (using ARTS' scat_za_grid means 0/180deg are
  // grid points, 90deg shouldn't be among them (resulting from even number
  // requirement for Disort and the (implicitly assumed?) requirement of a
  // horizon-symmetric za_grid)).
  //
  // We do all frequencies here at once (assuming this is the faster variant as
  // the agenda anyway (can) provide output for full f_grid at once and as we
  // have to apply the same inter/extrapolation to all the frequencies).

  CREATE_OUT2;

  chk_not_empty("surface_rtprop_agenda", surface_rtprop_agenda);

  const Index nf = f_grid.nelem();
  Index frza = 0;
  while (frza < scat_za_grid.nelem() && scat_za_grid[frza] < 90.) frza++;
  if (frza == scat_za_grid.nelem()) {
    ostringstream os;
    os << "No upwelling direction found in scat_za_grid.\n";
    throw runtime_error(os.str());
  }
  const Index nrza = scat_za_grid.nelem() - frza;
  //cout << nrza << " upwelling directions found, starting from element #"
  //     << frza << " of scat_za_grid.\n";
  Matrix dir_refl_coeff(nrza, nf, 0.);

  // Local input of surface_rtprop_agenda.
  Vector rtp_pos(1, surf_alt);  //atmosphere_dim is 1

  // first derive the (reflected-)direction dependent power reflection
  // coefficient
  for (Index rza = 0; rza < nrza; rza++) {
    // Local output of surface_rtprop_agenda.
    Numeric surface_skin_t;
    Matrix surface_los;
    Tensor4 surface_rmatrix;
    Matrix surface_emission;

    Vector rtp_los(1, scat_za_grid[rza + frza]);
    out2 << "Doing reflected dir #" << rza << " at " << rtp_los[0] << " degs\n";

    surface_rtprop_agendaExecute(ws,
                                 surface_skin_t,
                                 surface_emission,
                                 surface_los,
                                 surface_rmatrix,
                                 f_grid,
                                 rtp_pos,
                                 rtp_los,
                                 surface_rtprop_agenda);
    //cout << "surf_los has " << surface_los.ncols() << " columns and "
    //     << surface_los.nrows() << " rows.\n";
    ARTS_ASSERT(surface_los.ncols() == 1 || surface_los.nrows() == 0);
    if (rza == 0)
      btemp = surface_skin_t;
    else if (surface_skin_t != btemp) {
      ostringstream os;
      os << "Something went wrong.\n"
         << "  *surface_rtprop_agenda* returned different surface_skin_t\n"
         << "  for different LOS.\n";
      throw runtime_error(os.str());
    }
    if (surface_los.nrows() > 0) {
      for (Index f_index = 0; f_index < nf; f_index++)
        dir_refl_coeff(rza, f_index) =
            surface_rmatrix(joker, f_index, 0, 0).sum();
    }
    out2 << "  directional albedos[f_grid] = " << dir_refl_coeff(rza, joker)
         << "\n";
  }

  if (btemp < 0. || btemp > 1000.) {
    ostringstream os;
    os << "Surface temperature has been derived as " << btemp << " K,\n"
       << "which is not considered a meaningful value.\n";
    throw runtime_error(os.str());
  }

  // now integrate/average the (reflected-)direction dependent power reflection
  // coefficients
  //
  // starting with deriving the angles defining the angle ranges
  Vector surf_int_grid(nrza + 1);
  // the first angle grid point should be around (but above) 90deg and should
  // cover the angle range between the 90deg and half-way point towards the next
  // angle grid point. we probably also want to check, that we don't
  // 'extrapolate' too much.
  if (is_same_within_epsilon(scat_za_grid[frza], 90., 1e-6)) {
    ostringstream os;
    os << "Looks like scat_za_grid contains the 90deg direction,\n"
       << "which it shouldn't for running Disort.\n";
    throw runtime_error(os.str());
  }
  Numeric za_extrapol = (scat_za_grid[frza] - 90.) /
                        (scat_za_grid[frza + 1] - scat_za_grid[frza]);
  const Numeric ok_extrapol = 0.5;
  if ((za_extrapol - ok_extrapol) > 1e-6) {
    ostringstream os;
    os << "Extrapolation range from shallowest scat_za_grid point\n"
       << "to horizon is too big.\n"
       << "  Allowed extrapolation factor is 0.5.\n  Yours is " << za_extrapol
       << ", which is " << za_extrapol - 0.5 << " too big.\n";
    throw runtime_error(os.str());
  }
  if (!is_same_within_epsilon(
          scat_za_grid[scat_za_grid.nelem() - 1], 180., 1e-6)) {
    ostringstream os;
    os << "Looks like last point in scat_za_grid is not 180deg.\n";
    throw runtime_error(os.str());
  }

  surf_int_grid[0] = 90.;
  surf_int_grid[nrza] = 180.;
  for (Index rza = 1; rza < nrza; rza++)
    surf_int_grid[rza] =
        0.5 * (scat_za_grid[frza + rza - 1] + scat_za_grid[frza + rza]);
  surf_int_grid *= DEG2RAD;

  // now calculating the actual weights and apply them
  for (Index rza = 0; rza < nrza; rza++) {
    //Numeric coslow = cos(2.*surf_int_grid[rza]);
    //Numeric coshigh = cos(2.*surf_int_grid[rza+1]);
    //Numeric w = 0.5*(coshigh-coslow);
    Numeric w =
        0.5 * (cos(2. * surf_int_grid[rza + 1]) - cos(2. * surf_int_grid[rza]));
    //cout << "at reflLOS[" << rza << "]=" << scat_za_grid[frza+rza] << ":\n";
    //cout << "  angle weight derives as w = 0.5*(" << coshigh << "-"
    //     << coslow << ") = " << w << "\n";
    //cout << "  weighting directional reflection coefficient from "
    //     <<  dir_refl_coeff(rza,joker);
    dir_refl_coeff(rza, joker) *= w;
    //cout << " to " <<  dir_refl_coeff(rza,joker) << "\n";
  }

  // eventually sum up the weighted directional power reflection coefficients
  for (Index f_index = 0; f_index < nf; f_index++) {
    albedo[f_index] = dir_refl_coeff(joker, f_index).sum();
    out2 << "at f=" << f_grid[f_index] * 1e-9
         << " GHz, ending up with albedo=" << albedo[f_index] << "\n";
    if (albedo[f_index] < 0 || albedo[f_index] > 1.) {
      ostringstream os;
      os << "Something went wrong: Albedo must be inside [0,1],\n"
         << "  but is not at freq #" << f_index << " , where it is "
         << albedo[f_index] << ".\n";
      throw runtime_error(os.str());
    }
  }
}
