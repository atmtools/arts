/* Copyright (C) 2003-2012 Cory Davis <cory@met.ed.ac.uk>
                            
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
 * @file   m_montecarlo.cc
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-19
 *
 * @brief  Workspace functions for the solution of cloud-box radiative transfer
 * by Monte Carlo methods.  All of these functions refer to 3D calculations
 *
 */
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include "arts.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "mc_interp.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath_OLD.h"
#include "refraction.h"
#include "rng.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric BOLTZMAN_CONST=Constant::boltzmann_constant;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussian(MCAntenna& mc_antenna,
                           //keyword arguments
                           const Numeric& za_sigma,
                           const Numeric& aa_sigma,
                           const Verbosity&) {
  mc_antenna.set_gaussian(za_sigma, aa_sigma);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetGaussianByFWHM(MCAntenna& mc_antenna,
                                 //keyword arguments
                                 const Numeric& za_fwhm,
                                 const Numeric& aa_fwhm,
                                 const Verbosity&) {
  mc_antenna.set_gaussian_fwhm(za_fwhm, aa_fwhm);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void mc_antennaSetPencilBeam(MCAntenna& mc_antenna, const Verbosity&) {
  mc_antenna.set_pencil_beam();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MCGeneral(Workspace& ws,
               Vector& y,
               Index& mc_iteration_count,
               Vector& mc_error,
               Tensor3& mc_points,
               ArrayOfIndex& mc_source_domain,
               ArrayOfIndex& mc_scat_order,
               const MCAntenna& mc_antenna,
               const Vector& f_grid,
               const Index& f_index,
               const Matrix& sensor_pos,
               const Matrix& sensor_los,
               const Index& stokes_dim,
               const Index& atmosphere_dim,
               const Agenda& ppath_step_agenda,
               const Numeric& ppath_lmax,
               const Numeric& ppath_lraytrace,
               const Agenda& iy_space_agenda,
               const Agenda& surface_rtprop_agenda,
               const Agenda& propmat_clearsky_agenda,
               const Vector& p_grid,
               const Vector& lat_grid,
               const Vector& lon_grid,
               const Tensor3& z_field,
               const Vector& refellipsoid,
               const Matrix& z_surface,
               const Tensor3& t_field,
               const Tensor4& vmr_field,
               const Index& cloudbox_on,
               const ArrayOfIndex& cloudbox_limits,
               const Tensor4& pnd_field,
               const ArrayOfArrayOfSingleScatteringData& scat_data,
               const Index& atmfields_checked,
               const Index& atmgeom_checked,
               const Index& scat_data_checked,
               const Index& cloudbox_checked,
               const String& iy_unit,
               const Index& mc_seed,
               const Numeric& std_err,
               const Index& max_time,
               const Index& max_iter,
               const Index& min_iter,
               const Numeric& taustep_limit,
               const Index& l_mc_scat_order,
               const Index& t_interp_order,
               const Verbosity& verbosity) {
  // Checks of input
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  if (atmgeom_checked != 1)
    throw runtime_error(
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");

  if (!cloudbox_on)
    throw runtime_error("The cloudbox  must be activated (cloudbox_on=1).");
  if (cloudbox_checked != 1)
    throw runtime_error(
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  if (scat_data_checked != 1)
    throw runtime_error(
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  if (min_iter < 100) throw runtime_error("*mc_min_iter* must be >= 100.");

  if (max_time < 0 && max_iter < 0 && std_err < 0)
    throw runtime_error(
        "At least one of std_err, max_time, and max_iter "
        "must be positive.");

  if (l_mc_scat_order <= 0)
    throw runtime_error("*l_max_scat_order* must be > 0.");

  if (f_index < 0)
    throw runtime_error(
        "The option of f_index < 0 is not handled by this "
        "method.");
  if (f_index >= f_grid.nelem())
    throw runtime_error("*f_index* is outside the range of *f_grid*.");

  if (atmosphere_dim != 3)
    throw runtime_error("Only 3D atmospheres are handled. ");

  if (sensor_pos.ncols() != 3) {
    ostringstream os;
    os << "Expected number of columns in sensor_pos: 3.\n";
    os << "Found: " << sensor_pos.ncols();
    throw runtime_error(os.str());
  }
  if (sensor_pos.nrows() != 1) {
    ostringstream os;
    os << "Expected number of rows in sensor_pos: 1.\n";
    os << "Found: " << sensor_pos.nrows();
    throw runtime_error(os.str());
  }

  if (!(sensor_los.ncols() == 2)) {
    ostringstream os;
    os << "Expected number of columns in sensor_los: 2.\n";
    os << "Found: " << sensor_los.ncols();
    throw runtime_error(os.str());
  }
  if (!(sensor_los.nrows() == 1)) {
    ostringstream os;
    os << "Expected number of rows in sensor_los: 1.\n";
    os << "Found: " << sensor_los.nrows();
    throw runtime_error(os.str());
  }

  Ppath ppath_step;
  Rng rng;  //Random Number generator
  time_t start_time = time(NULL);
  Index N_se = pnd_field.nbooks();  //Number of scattering elements
  Vector pnd_vec(
      N_se);  //Vector of particle number densities used at each point
  Vector Z11maxvector(
      N_se);  //Vector holding the maximum phase function for each

  // finding maximum phase function for each scat element
  Index i_total = -1;
  Index this_f_index = min(f_index, scat_data[0][0].f_grid.nelem());
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
      i_total++;
      Z11maxvector[i_total] = max(scat_data[i_ss][i_se].pha_mat_data(
          this_f_index, joker, joker, joker, joker, joker, 0));
    }
  }

  rng.seed(mc_seed, verbosity);
  Numeric g, temperature, albedo, g_los_csc_theta;
  Matrix A(stokes_dim, stokes_dim), Q(stokes_dim, stokes_dim);
  Matrix evol_op(stokes_dim, stokes_dim), ext_mat_mono(stokes_dim, stokes_dim);
  Matrix q(stokes_dim, stokes_dim), newQ(stokes_dim, stokes_dim);
  Matrix Z(stokes_dim, stokes_dim);
  Matrix R_ant2enu(3, 3),
      R_stokes(stokes_dim, stokes_dim);  // Needed for antenna rotations
  q = 0.0;
  newQ = 0.0;
  Vector vector1(stokes_dim), abs_vec_mono(stokes_dim), I_i(stokes_dim);
  Vector Isum(stokes_dim), Isquaredsum(stokes_dim);
  Index termination_flag = 0;
  const Numeric f_mono = f_grid[f_index];
  const Numeric prop_dir =
      -1.0;  // propagation direction opposite of los angles

  CREATE_OUT0;

  y.resize(stokes_dim);
  y = 0;

  mc_iteration_count = 0;
  mc_error.resize(stokes_dim);
  mc_points.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
  mc_points = 0;
  mc_scat_order.resize(l_mc_scat_order);
  mc_scat_order = 0;
  mc_source_domain.resize(4);
  mc_source_domain = 0;

  //local versions of workspace
  Numeric local_surface_skin_t;
  Matrix local_iy(1, stokes_dim), local_surface_emission(1, stokes_dim);
  Matrix local_surface_los;
  Tensor4 local_surface_rmatrix;
  Vector local_rte_pos(3);  // Fixed this (changed from 2 to 3)
  Vector local_rte_los(2);
  Vector new_rte_los(2);
  Index np;
  Isum = 0.0;
  Isquaredsum = 0.0;
  Numeric std_err_i;
  bool convert_to_rjbt = false;
  if (iy_unit == "RJBT") {
    std_err_i = f_mono * f_mono * 2 * BOLTZMAN_CONST / SPEED_OF_LIGHT /
                SPEED_OF_LIGHT * std_err;
    convert_to_rjbt = true;
  } else if (iy_unit == "1") {
    std_err_i = std_err;
  } else {
    ostringstream os;
    os << "Invalid value for *iy_unit*:" << iy_unit << ".\n"
       << "This method allows only the options \"RJBT\" and \"1\".";
    throw runtime_error(os.str());
  }

  // Calculate rotation matrix for boresight
  rotmat_enu(R_ant2enu, sensor_los(0, joker));

  //Begin Main Loop
  //
  bool keepgoing, oksampling;
  Index nfails = 0;
  //
  while (true) {
    // Complete content of while inside try/catch to handle occasional
    // failures in the ppath calculations
    try {
      bool inside_cloud;

      mc_iteration_count += 1;
      Index scattering_order = 0;

      keepgoing = true;   // indicating whether to continue tracing a photon
      oksampling = true;  // gets false if g becomes zero

      //Sample a FOV direction
      Matrix R_prop(3, 3);
      mc_antenna.draw_los(
          local_rte_los, R_prop, rng, R_ant2enu, sensor_los(0, joker));

      // Get stokes rotation matrix for rotating polarization
      rotmat_stokes(
          R_stokes, stokes_dim, prop_dir, prop_dir, R_prop, R_ant2enu);
      id_mat(Q);
      local_rte_pos = sensor_pos(0, joker);
      I_i = 0.0;

      while (keepgoing) {
        mcPathTraceGeneral(ws,
                           evol_op,
                           abs_vec_mono,
                           temperature,
                           ext_mat_mono,
                           rng,
                           local_rte_pos,
                           local_rte_los,
                           pnd_vec,
                           g,
                           ppath_step,
                           termination_flag,
                           inside_cloud,
                           ppath_step_agenda,
                           ppath_lmax,
                           ppath_lraytrace,
                           taustep_limit,
                           propmat_clearsky_agenda,
                           stokes_dim,
                           f_index,
                           f_grid,
                           p_grid,
                           lat_grid,
                           lon_grid,
                           z_field,
                           refellipsoid,
                           z_surface,
                           t_field,
                           vmr_field,
                           cloudbox_limits,
                           pnd_field,
                           scat_data,
                           verbosity);

        // GH 2011-09-08: if the lowest layer has large
        // extent and a thick cloud, g may be 0 due to
        // underflow, but then I_i should be 0 as well.
        // Don't turn it into nan for no reason.
        // If reaching underflow, no point in going on;
        // hence new photon.
        // GH 2011-09-14: moved this check to outside the different
        // scenarios, as this goes wrong regardless of the scenario.
        if (g == 0) {
          keepgoing = false;
          oksampling = false;
          mc_iteration_count -= 1;
          out0 << "WARNING: A rejected path sampling (g=0)!\n(if this"
               << "happens repeatedly, try to decrease *ppath_lmax*)";
        } else if (termination_flag == 1) {
          iy_space_agendaExecute(ws,
                                 local_iy,
                                 Vector(1, f_mono),
                                 local_rte_pos,
                                 local_rte_los,
                                 iy_space_agenda);
          mult(vector1, evol_op, local_iy(0, joker));
          mult(I_i, Q, vector1);
          I_i /= g;
          keepgoing = false;  //stop here. New photon.
          mc_source_domain[0] += 1;
        } else if (termination_flag == 2) {
          //Calculate surface properties
          surface_rtprop_agendaExecute(ws,
                                       local_surface_skin_t,
                                       local_surface_emission,
                                       local_surface_los,
                                       local_surface_rmatrix,
                                       Vector(1, f_mono),
                                       local_rte_pos,
                                       local_rte_los,
                                       surface_rtprop_agenda);

          //if( local_surface_los.nrows() > 1 )
          // throw runtime_error(
          //                "The method handles only specular reflections." );

          //deal with blackbody case
          if (local_surface_los.empty()) {
            mult(vector1, evol_op, local_surface_emission(0, joker));
            mult(I_i, Q, vector1);
            I_i /= g;
            keepgoing = false;
            mc_source_domain[1] += 1;
          } else
          //decide between reflection and emission
          {
            const Numeric rnd = rng.draw();

            Numeric R11 = 0;
            for (Index i = 0; i < local_surface_rmatrix.nbooks(); i++) {
              R11 += local_surface_rmatrix(i, 0, 0, 0);
            }

            if (rnd > R11) {
              //then we have emission
              mult(vector1, evol_op, local_surface_emission(0, joker));
              mult(I_i, Q, vector1);
              I_i /= g * (1 - R11);
              keepgoing = false;
              mc_source_domain[1] += 1;
            } else {
              //we have reflection
              // determine which reflection los to use
              Index i = 0;
              Numeric rsum = local_surface_rmatrix(i, 0, 0, 0);
              while (rsum < rnd) {
                i++;
                rsum += local_surface_rmatrix(i, 0, 0, 0);
              }

              local_rte_los = local_surface_los(i, joker);

              mult(q, evol_op, local_surface_rmatrix(i, 0, joker, joker));
              mult(newQ, Q, q);
              Q = newQ;
              Q /= g * local_surface_rmatrix(i, 0, 0, 0);
            }
          }
        } else if (inside_cloud) {
          //we have another scattering/emission point
          //Estimate single scattering albedo
          albedo = 1 - abs_vec_mono[0] / ext_mat_mono(0, 0);

          //determine whether photon is emitted or scattered
          if (rng.draw() > albedo) {
            //Calculate emission
            Numeric planck_value = planck(f_mono, temperature);
            Vector emission = abs_vec_mono;
            emission *= planck_value;
            Vector emissioncontri(stokes_dim);
            mult(emissioncontri, evol_op, emission);
            emissioncontri /= (g * (1 - albedo));  //yuck!
            mult(I_i, Q, emissioncontri);
            keepgoing = false;
            mc_source_domain[3] += 1;
          } else {
            //we have a scattering event
            Sample_los(new_rte_los,
                       g_los_csc_theta,
                       Z,
                       rng,
                       local_rte_los,
                       scat_data,
                       f_index,
                       stokes_dim,
                       pnd_vec,
                       Z11maxvector,
                       ext_mat_mono(0, 0) - abs_vec_mono[0],
                       temperature,
                       t_interp_order);

            Z /= g * g_los_csc_theta * albedo;

            mult(q, evol_op, Z);
            mult(newQ, Q, q);
            Q = newQ;
            scattering_order += 1;
            local_rte_los = new_rte_los;
          }
        } else {
          //Must be clear sky emission point
          //Calculate emission
          Numeric planck_value = planck(f_mono, temperature);
          Vector emission = abs_vec_mono;
          emission *= planck_value;
          Vector emissioncontri(stokes_dim);
          mult(emissioncontri, evol_op, emission);
          emissioncontri /= g;
          mult(I_i, Q, emissioncontri);
          keepgoing = false;
          mc_source_domain[2] += 1;
        }
      }  // keepgoing

      if (oksampling) {
        // Set spome of the bookkeeping variables
        np = ppath_step.np;
        mc_points(ppath_step.gp_p[np - 1].idx,
                  ppath_step.gp_lat[np - 1].idx,
                  ppath_step.gp_lon[np - 1].idx) += 1;
        if (scattering_order < l_mc_scat_order) {
          mc_scat_order[scattering_order] += 1;
        }

        // Rotate into antenna polarization frame
        Vector I_hold(stokes_dim);
        mult(I_hold, R_stokes, I_i);
        Isum += I_i;

        for (Index j = 0; j < stokes_dim; j++) {
          ARTS_ASSERT(!std::isnan(I_i[j]));
          Isquaredsum[j] += I_i[j] * I_i[j];
        }
        y = Isum;
        y /= (Numeric)mc_iteration_count;
        for (Index j = 0; j < stokes_dim; j++) {
          mc_error[j] = sqrt(
              (Isquaredsum[j] / (Numeric)mc_iteration_count - y[j] * y[j]) /
              (Numeric)mc_iteration_count);
        }
        if (std_err > 0 && mc_iteration_count >= min_iter &&
            mc_error[0] < std_err_i) {
          break;
        }
        if (max_time > 0 && (Index)(time(NULL) - start_time) >= max_time) {
          break;
        }
        if (max_iter > 0 && mc_iteration_count >= max_iter) {
          break;
        }
      }
    }  // Try

    catch (const std::runtime_error& e) {
      mc_iteration_count += 1;
      nfails += 1;
      out0 << "WARNING: A MC path sampling failed! Error was:\n";
      cout << e.what() << endl;
      if (nfails >= 5) {
        throw runtime_error(
            "The MC path sampling has failed five times. A few failures "
            "should be OK, but this number is suspiciously high and the "
            "reason to these failures should be tracked down.");
      }
    }
  }  // while

  if (convert_to_rjbt) {
    for (Index j = 0; j < stokes_dim; j++) {
      y[j] = invrayjean(y[j], f_mono);
      mc_error[j] = invrayjean(mc_error[j], f_mono);
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void MCRadar(  // Workspace reference:
    Workspace& ws,
    // WS Output:
    Vector& y,
    Vector& mc_error,
    // WS Input:
    const MCAntenna& mc_antenna,
    const Vector& f_grid,
    const Index& f_index,
    const Matrix& sensor_pos,
    const Matrix& sensor_los,
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Numeric& ppath_lmax,
    const Agenda& ppath_step_agenda,
    const Numeric& ppath_lraytrace,
    const Agenda& propmat_clearsky_agenda,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    const Vector& refellipsoid,
    const Matrix& z_surface,
    const Tensor3& t_field,
    const Tensor4& vmr_field,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor4& pnd_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& mc_y_tx,
    const Vector& range_bins,
    const Index& atmfields_checked,
    const Index& atmgeom_checked,
    const Index& scat_data_checked,
    const Index& cloudbox_checked,
    const String& iy_unit_radar,
    const Index& mc_max_scatorder,
    const Index& mc_seed,
    const Index& mc_max_iter,
    const Numeric& ze_tref,
    const Numeric& k2,
    const Index& t_interp_order,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUT0;

  // Important constants
  const Index nbins = range_bins.nelem() - 1;
  const Numeric r_min = min(range_bins);
  const Numeric r_max = max(range_bins);

  // Basics
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  if (stokes_dim < 2)
    throw runtime_error("This method requires that stokes_dim >= 2");
  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  if (atmgeom_checked != 1)
    throw runtime_error(
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");

  if (!cloudbox_on)
    throw runtime_error("The cloudbox  must be activated (cloudbox_on=1).");
  if (cloudbox_checked != 1)
    throw runtime_error(
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  if (scat_data_checked != 1)
    throw runtime_error(
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  if (f_index < 0)
    throw runtime_error(
        "The option of f_index < 0 is not handled by this "
        "method.");
  if (f_index >= f_grid.nelem())
    throw runtime_error("*f_index* is outside the range of *f_grid*.");

  if (atmosphere_dim != 3)
    throw runtime_error("Only 3D atmospheres are handled.");

  if (stokes_dim != mc_y_tx.nelem())
    throw runtime_error("*mc_y_tx* must have size of *stokes_dim*.");

  if (sensor_pos.ncols() != 3) {
    ostringstream os;
    os << "Expected number of columns in sensor_pos: 3.\n";
    os << "Found: " << sensor_pos.ncols();
    throw runtime_error(os.str());
  }

  if (sensor_los.ncols() != 2) {
    ostringstream os;
    os << "Expected number of columns in sensor_los: 2.\n";
    os << "Found: " << sensor_los.ncols();
    throw runtime_error(os.str());
  }

  if (mc_max_iter < 0)
    throw runtime_error(
        "mc_max_iter must be positive, "
        "as it is the only limiter.");

  if (!is_increasing(range_bins))
    throw runtime_error(
        "The vector *range_bins* must contain strictly "
        "increasing values.");

  if (r_min < 0)
    throw runtime_error(
        "The vector *range_bins* is not allowed to contain "
        "negative distance or round-trip time.");

  if (mc_antenna.atype != ANTENNA_TYPE_GAUSSIAN) {
    throw runtime_error(
        "MCRadar only works with "
        "Gaussian antenna patterns.");
  }

  Ppath ppath_step;
  Rng rng;                          //Random Number generator
  Index N_se = pnd_field.nbooks();  //Number of scattering elements
  Vector pnd_vec(
      N_se);  //Vector of particle number densities used at each point
  bool anyptype_nonTotRan = is_anyptype_nonTotRan(scat_data);
  bool is_dist = max(range_bins) > 1;  // Is it round trip time or distance
  rng.seed(mc_seed, verbosity);
  Numeric ppath_lraytrace_var;
  //Numeric temperature, albedo;
  Numeric albedo;
  Numeric Csca, Cext;
  Numeric antenna_wgt;
  Matrix evol_op(stokes_dim, stokes_dim), ext_mat_mono(stokes_dim, stokes_dim);
  Matrix trans_mat(stokes_dim, stokes_dim);
  Matrix Z(stokes_dim, stokes_dim);
  Matrix R_ant2enu(3, 3), R_enu2ant(3, 3), R_stokes(stokes_dim, stokes_dim);
  Vector abs_vec_mono(stokes_dim), I_i(stokes_dim), I_i_rot(stokes_dim);
  Vector Isum(nbins * stokes_dim), Isquaredsum(nbins * stokes_dim);
  Index termination_flag = 0;
  Vector bin_height(nbins);
  Vector range_bin_count(nbins);
  Index mc_iter;
  Index scat_order;

  // allocating variables needed for pha_mat extraction (don't want to do this
  // in every loop step again).
  ArrayOfArrayOfTensor6 pha_mat_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor6 pha_mat_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor6 pha_mat_bulk;
  Index ptype_bulk;
  Matrix pdir_array(1, 2), idir_array(1, 2);
  Vector t_array(1);
  Matrix pnds(N_se, 1);

  // for pha_mat handling, at the moment we still need scat_data_mono. Hence,
  // extract that here (but in its local container, not into the WSV
  // scat_data_mono).
  //ArrayOfArrayOfSingleScatteringData this_scat_data_mono;
  //scat_data_monoExtract( this_scat_data_mono, scat_data, f_index, verbosity );

  const Numeric f_mono = f_grid[f_index];
  const Numeric tx_dir = 1.0;
  const Numeric rx_dir = -1.0;

  for (Index ibin = 0; ibin < nbins; ibin++) {
    bin_height[ibin] = range_bins[ibin + 1] - range_bins[ibin];
  }
  if (!is_dist) {
    bin_height *= 0.5 * SPEED_OF_LIGHT;
  }

  y.resize(nbins * stokes_dim);
  y = 0;

  range_bin_count = 0;

  mc_iter = 0;
  // this will need to be reshaped differently for range gates
  mc_error.resize(stokes_dim * nbins);

  //local versions of workspace
  Matrix local_iy(1, stokes_dim), local_surface_emission(1, stokes_dim);
  Matrix local_surface_los;
  Tensor4 local_surface_rmatrix;
  Vector local_rte_pos(3);
  Vector local_rte_los(2);
  Vector new_rte_los(2);
  Vector Ipath(stokes_dim), Ihold(stokes_dim), Ipath_norm(stokes_dim);
  Isum = 0.0;
  Isquaredsum = 0.0;
  Numeric s_tot, s_return;  // photon distance traveled
  Numeric t_tot, t_return;  // photon time traveled
  Numeric r_trav, r_bin;  // range traveled (1-way distance) or round-trip time

  Numeric fac;
  if (iy_unit_radar == "1") {
    fac = 1.0;
  }
  // Conversion from intensity to reflectivity
  else if (iy_unit_radar == "Ze") {
    Vector cfac(1);
    ze_cfac(cfac, Vector(1, f_mono), ze_tref, k2);
    // Due to different definitions, the factor shall here be scaled with 1/(2pi)
    fac = cfac[0] / (2 * PI);
  } else {
    ostringstream os;
    os << "Invalid value for *iy_unit_radar*:" << iy_unit_radar << ".\n"
       << "This method allows only the options \"Ze\" and \"1\".";
    throw runtime_error(os.str());
  }

  // Calculate rotation matrix and polarization bases for boresight
  rotmat_enu(R_ant2enu, sensor_los(0, joker));
  R_enu2ant = transpose(R_ant2enu);

  //Begin Main Loop
  bool keepgoing, firstpass, integrity;
  while (mc_iter < mc_max_iter) {
    bool inside_cloud;

    mc_iter += 1;

    integrity = true;  // intensity is not nan or below threshold
    keepgoing = true;  // indicating whether to continue tracing a photon
    firstpass = true;  // ensure backscatter is properly calculated

    //Sample a FOV direction
    Matrix R_tx(3, 3);
    mc_antenna.draw_los(
        local_rte_los, R_tx, rng, R_ant2enu, sensor_los(0, joker));
    rotmat_stokes(R_stokes, stokes_dim, tx_dir, tx_dir, R_ant2enu, R_tx);
    mult(Ihold, R_stokes, mc_y_tx);

    // Initialize other variables
    local_rte_pos = sensor_pos(0, joker);
    s_tot = 0.0;
    t_tot = 0.0;
    scat_order = 0;
    while (keepgoing) {
      Numeric s_path, t_path;

      mcPathTraceRadar(ws,
                       evol_op,
                       abs_vec_mono,
                       t_array[0],
                       ext_mat_mono,
                       rng,
                       local_rte_pos,
                       local_rte_los,
                       pnd_vec,  //pnds(joker,0),
                       s_path,
                       t_path,
                       ppath_step,
                       termination_flag,
                       inside_cloud,
                       ppath_step_agenda,
                       ppath_lmax,
                       ppath_lraytrace,
                       propmat_clearsky_agenda,
                       anyptype_nonTotRan,
                       stokes_dim,
                       f_index,
                       f_grid,
                       Ihold,
                       p_grid,
                       lat_grid,
                       lon_grid,
                       z_field,
                       refellipsoid,
                       z_surface,
                       t_field,
                       vmr_field,
                       cloudbox_limits,
                       pnd_field,
                       scat_data,
                       verbosity);
      pnds(joker, 0) = pnd_vec;
      if (!inside_cloud || termination_flag != 0) {
        keepgoing = false;
      } else {
        s_tot += s_path;
        t_tot += t_path;

        //
        Csca = ext_mat_mono(0, 0) - abs_vec_mono[0];
        Cext = ext_mat_mono(0, 0);
        if (anyptype_nonTotRan) {
          const Numeric Irat = Ihold[1] / Ihold[0];
          Csca += Irat * (ext_mat_mono(1, 0) - abs_vec_mono[1]);
          Cext += Irat * ext_mat_mono(0, 1);
        }
        albedo = Csca / Cext;

        // Terminate if absorption event, outside cloud, or surface
        Numeric rn = rng.draw();
        if (rn > albedo) {
          keepgoing = false;
          continue;
        }

        Vector rte_los_geom(2);

        // Compute reflectivity contribution based on local-to-sensor
        // geometry, path attenuation
        // Get los angles at atmospheric locale to determine
        // scattering angle
        if (firstpass) {
          // Use this to ensure that the difference in azimuth angle
          // between incident and scattered lines-of-sight is 180
          // degrees
          mirror_los(rte_los_geom, local_rte_los, atmosphere_dim);
          firstpass = false;
        } else {
          // Replace with ppath_agendaExecute??
          rte_losGeometricFromRtePosToRtePos2(rte_los_geom,
                                              atmosphere_dim,
                                              lat_grid,
                                              lon_grid,
                                              refellipsoid,
                                              local_rte_pos,
                                              Vector{sensor_pos(0, joker)},
                                              verbosity);
        }

        // Get los angles at sensor to determine antenna pattern
        // weighting of return signal and ppath to determine
        // propagation path back to sensor
        // Replace with ppath_agendaExecute??
        Ppath ppath;
        Vector rte_los_antenna(2);
        ppath_lraytrace_var = ppath_lraytrace;
        Numeric za_accuracy = 2e-5;
        Numeric pplrt_factor = 5;
        Numeric pplrt_lowest = 0.5;

        rte_losGeometricFromRtePosToRtePos2(rte_los_antenna,
                                            atmosphere_dim,
                                            lat_grid,
                                            lon_grid,
                                            refellipsoid,
                                            Vector{sensor_pos(0, joker)},
                                            local_rte_pos,
                                            verbosity);

        ppathFromRtePos2(ws,
                         ppath,
                         rte_los_antenna,
                         ppath_lraytrace_var,
                         ppath_step_agenda,
                         atmosphere_dim,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         z_field,
                         f_grid,
                         refellipsoid,
                         z_surface,
                         Vector{sensor_pos(0, joker)},
                         local_rte_pos,
                         ppath_lmax,
                         za_accuracy,
                         pplrt_factor,
                         pplrt_lowest,
                         verbosity);

        // Return distance
        const Index np2 = ppath.np;
        s_return = ppath.end_lstep;
        t_return = s_return / SPEED_OF_LIGHT;
        for (Index ip = 1; ip < np2; ip++) {
          s_return += ppath.lstep[ip - 1];
          t_return += ppath.lstep[ip - 1] * 0.5 *
                      (ppath.ngroup[ip - 1] + ppath.ngroup[ip]) /
                      SPEED_OF_LIGHT;
        }

        // One-way distance
        if (is_dist) {
          r_trav = 0.5 * (s_tot + s_return);
        }

        // Round trip travel time
        else {
          r_trav = t_tot + t_return;
        }

        // Still within max range of radar?
        if (r_trav <= r_max) {
          // Compute path extinction as with radio link
          get_ppath_transmat(ws,
                             trans_mat,
                             ppath,
                             propmat_clearsky_agenda,
                             stokes_dim,
                             f_index,
                             f_grid,
                             p_grid,
                             t_field,
                             vmr_field,
                             cloudbox_limits,
                             pnd_field,
                             scat_data,
                             verbosity);

          // Obtain scattering matrix given incident and scattered angles
          Matrix P(stokes_dim, stokes_dim);

          pdir_array(0, joker) = rte_los_geom;
          idir_array(0, joker) = local_rte_los;
          pha_mat_NScatElems(pha_mat_Nse,
                             ptypes_Nse,
                             t_ok,
                             scat_data,
                             stokes_dim,
                             t_array,
                             pdir_array,
                             idir_array,
                             f_index,
                             t_interp_order);
          pha_mat_ScatSpecBulk(pha_mat_ssbulk,
                               ptype_ssbulk,
                               pha_mat_Nse,
                               ptypes_Nse,
                               pnds,
                               t_ok);
          pha_mat_Bulk(pha_mat_bulk, ptype_bulk, pha_mat_ssbulk, ptype_ssbulk);
          P = pha_mat_bulk(0, 0, 0, 0, joker, joker);

          P *= 4 * PI;
          P /= Csca;

          // Compute reflectivity contribution here
          mult(Ipath, evol_op, Ihold);
          Ipath /= Ipath[0];
          Ipath *= Ihold[0];
          mult(Ihold, P, Ipath);
          mult(I_i, trans_mat, Ihold);
          Ihold = Ipath;
          if (Ihold[0] < 1e-40 || std::isnan(Ihold[0]) ||
              std::isnan(Ihold[1]) ||
              (stokes_dim > 2 && std::isnan(Ihold[2])) ||
              (stokes_dim > 3 && std::isnan(Ihold[3]))) {
            integrity = false;
          }

          if (r_trav > r_min && integrity) {
            // Add reflectivity to proper range bin
            Index ibin = 0;
            r_bin = 0.0;
            while (r_bin < r_trav && ibin <= nbins + 1) {
              ibin++;
              r_bin = range_bins[ibin];
            }
            ibin -= 1;

            // Calculate rx antenna weight and polarization rotation
            Matrix R_rx(3, 3);
            rotmat_enu(R_rx, rte_los_antenna);
            mc_antenna.return_los(antenna_wgt, R_rx, R_enu2ant);
            rotmat_stokes(
                R_stokes, stokes_dim, rx_dir, tx_dir, R_rx, R_ant2enu);
            mult(I_i_rot, R_stokes, I_i);

            for (Index istokes = 0; istokes < stokes_dim; istokes++) {
              Index ibiny = ibin * stokes_dim + istokes;
              ARTS_ASSERT(!std::isnan(I_i_rot[istokes]));
              Isum[ibiny] += antenna_wgt * I_i_rot[istokes];
              Isquaredsum[ibiny] += antenna_wgt * antenna_wgt *
                                    I_i_rot[istokes] * I_i_rot[istokes];
            }
            range_bin_count[ibin] += 1;
          }

          scat_order++;
            
          Sample_los_uniform(new_rte_los, rng);
          pdir_array(0, joker) = new_rte_los;
          // alt:
          // Sample_los_uniform( pdir_array(0,joker), rng );
          pha_mat_NScatElems(pha_mat_Nse,
                             ptypes_Nse,
                             t_ok,
                             scat_data,
                             stokes_dim,
                             t_array,
                             pdir_array,
                             idir_array,
                             f_index,
                             t_interp_order);
          pha_mat_ScatSpecBulk(pha_mat_ssbulk,
                               ptype_ssbulk,
                               pha_mat_Nse,
                               ptypes_Nse,
                               pnds,
                               t_ok);
          pha_mat_Bulk(pha_mat_bulk, ptype_bulk, pha_mat_ssbulk, ptype_ssbulk);
          Z = pha_mat_bulk(0, 0, 0, 0, joker, joker);

          Z *= 4 * PI;
          Z /= Csca;
          mult(Ipath, Z, Ihold);
          Ihold = Ipath;
          local_rte_los = new_rte_los;
          // alt:
          //local_rte_los = pdir_array(0,joker);
          // or even (but also requires replacements of local_rte_los
          // with idir_array throughout the whole loop):
          //idir_array = pdir_array;
        } else {
          // Past farthest range
          keepgoing = false;
        }
      }

      // Some checks
      if (scat_order >= mc_max_scatorder) keepgoing = false;
      if (!integrity) keepgoing = false;
    }  // while (inner: keepgoing)

  }  // while (outer)

  // Normalize range bins and apply sensor response (polarization)
  for (Index ibin = 0; ibin < nbins; ibin++) {
    for (Index istokes = 0; istokes < stokes_dim; istokes++) {
      Index ibiny = ibin * stokes_dim + istokes;
      if (range_bin_count[ibin] > 0) {
        y[ibiny] = Isum[ibiny] / ((Numeric)mc_iter) / bin_height[ibin];
        mc_error[ibiny] = sqrt((Isquaredsum[ibiny] / (Numeric)mc_iter /
                                    bin_height[ibin] / bin_height[ibin] -
                                y[ibiny] * y[ibiny]) /
                               (Numeric)mc_iter);
      }
    }
  }

  y *= fac;
  mc_error *= fac;
}  // end MCRadar

/* Workspace method: Doxygen documentation will be auto-generated */
void MCSetSeedFromTime(Index& mc_seed, const Verbosity&) {
  mc_seed = (Index)time(NULL);
}
