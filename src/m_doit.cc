/*!
  \file   m_doit.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \author Sreerekha T.R. <rekha@uni-bremen.de>
  \date   Wed Jun 19 11:03:57 2002
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox using the DOIT method.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <workspace.h>
#include "array.h"
#include "arts.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "check_input.h"
#include "debug.h"
#include "doit.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "m_general.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "species_tags.h"
#include "surf.h"
#include "xml_io.h"
#include "arts_omp.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void DOAngularGridsSet(  // WS Output:
    Index& doit_za_grid_size,
    Vector& aa_grid,
    Vector& za_grid,
    // Keywords:
    const Index& N_za_grid,
    const Index& N_aa_grid,
    const String& za_grid_opt_file) {
  // Azimuth angle grid (the same is used for the scattering integral and
  // for the radiative transfer.
  if (N_aa_grid > 1)
    nlinspace(aa_grid, 0, 360, N_aa_grid);
  else {
    ARTS_USER_ERROR_IF (N_aa_grid < 1,
      "N_aa_grid must be > 0 (even for 1D / DISORT cases).")
    aa_grid.resize(1);
    aa_grid[0] = 0.;
  }

  // Zenith angle grid:
  // Number of zenith angle grid points (only for scattering integral):
  ARTS_USER_ERROR_IF (N_za_grid < 0,
    "N_za_grid must be >= 0.")
  doit_za_grid_size = N_za_grid;

  if (za_grid_opt_file == "")
    if (N_za_grid == 0)
      za_grid.resize(0);
    else {
      ARTS_USER_ERROR_IF (N_za_grid == 1,
        "N_za_grid must be >1 or =0 (the latter only allowed for RT4).")
      nlinspace(za_grid, 0, 180, N_za_grid);
    }
  else
    xml_read_from_file(za_grid_opt_file, za_grid);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagAbs(  //WS Input and Output:
    Index& doit_conv_flag,
    Index& doit_iteration_counter,
    Tensor6& cloudbox_field_mono,
    // WS Input:
    const Tensor6& cloudbox_field_mono_old,
    // Keyword:
    const Vector& epsilon,
    const Index& max_iterations,
    const Index& throw_nonconv_error) {
  //------------Check the input---------------------------------------
  ARTS_USER_ERROR_IF (doit_conv_flag != 0,
        "Convergence flag is non-zero, which means that this\n"
        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
        "be used only in *doit_conv_test_agenda*\n");

  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_lat = cloudbox_field_mono.nshelves();
  const Index N_lon = cloudbox_field_mono.nbooks();
  const Index N_za = cloudbox_field_mono.npages();
  const Index N_aa = cloudbox_field_mono.nrows();

  // Check keyword "epsilon":
  ARTS_USER_ERROR_IF (epsilon.nelem() != 4,
        "You have to specify limiting values for the "
        "convergence test for each Stokes component "
        "separately. That means that *epsilon* must "
        "have *4* elements!");

  // Check if cloudbox_field and cloudbox_field_old have the same dimensions:
  ARTS_USER_ERROR_IF (!is_size(
          cloudbox_field_mono_old, N_p, N_lat, N_lon, N_za, N_aa, 4),
        "The fields (Tensor6) *cloudbox_field* and \n"
        "*cloudbox_field_old* which are compared in the \n"
        "convergence test do not have the same size.\n");

  //-----------End of checks-------------------------------------------------

  doit_iteration_counter += 1;

  if (doit_iteration_counter > max_iterations) {
    ostringstream out;
    out << "Method does not converge (number of iterations \n"
        << "is > " << max_iterations << "). Either the cloud "
        << "particle number density \n"
        << "is too large or the numerical setup for the DOIT \n"
        << "calculation is not correct. In case of limb \n"
        << "simulations please make sure that you use an \n"
        << "optimized zenith angle grid. \n"
        << "*cloudbox_field* might be wrong.\n";
    if (throw_nonconv_error != 0) {
      // FIXME: OLE: Remove this later
      //          ostringstream os;
      //          os << "Error in DOIT calculation:\n"
      //             << out.str();
      //          throw runtime_error( os.str() );
      cloudbox_field_mono = NAN;
      doit_conv_flag = 1;
    } else {
      doit_conv_flag = 1;
    }
  } else {
    for (Index p_index = 0; p_index < N_p; p_index++) {
      for (Index lat_index = 0; lat_index < N_lat; lat_index++) {
        for (Index lon_index = 0; lon_index < N_lon; lon_index++) {
          for (Index za_index = 0; za_index < N_za; za_index++) {
            for (Index aa_index = 0; aa_index < N_aa; aa_index++) {
              for (Index stokes_index = 0; stokes_index < 4;
                   stokes_index++) {
                Numeric diff = (cloudbox_field_mono(p_index,
                                                    lat_index,
                                                    lon_index,
                                                    za_index,
                                                    aa_index,
                                                    stokes_index) -
                                cloudbox_field_mono_old(p_index,
                                                        lat_index,
                                                        lon_index,
                                                        za_index,
                                                        aa_index,
                                                        stokes_index));

                // If the absolute difference of the components
                // is larger than the pre-defined values, return
                // to *cloudbox_fieldIterarte* and do next iteration

                if (abs(diff) > epsilon[stokes_index]) {
                  return;
                }

              }  // End loop stokes_dom.
            }    // End loop aa_grid.
          }      // End loop za_grid.
        }        // End loop lon_grid.
      }          // End loop lat_grid.
    }            // End p_grid.

    // Convergence test has been successful, doit_conv_flag can be set to 1.
    doit_conv_flag = 1;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagAbsBT(  //WS Input and Output:
    Index& doit_conv_flag,
    Index& doit_iteration_counter,
    Tensor6& cloudbox_field_mono,
    // WS Input:
    const Tensor6& cloudbox_field_mono_old,
    const Vector& f_grid,
    const Index& f_index,
    // Keyword:
    const Vector& epsilon,
    const Index& max_iterations,
    const Index& throw_nonconv_error) {
  //------------Check the input---------------------------------------

  ARTS_USER_ERROR_IF (doit_conv_flag != 0,
        "Convergence flag is non-zero, which means that this \n"
        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
        "be used only in *doit_conv_test_agenda*\n");

  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_lat = cloudbox_field_mono.nshelves();
  const Index N_lon = cloudbox_field_mono.nbooks();
  const Index N_za = cloudbox_field_mono.npages();
  const Index N_aa = cloudbox_field_mono.nrows();

  // Check keyword "epsilon":
  ARTS_USER_ERROR_IF (epsilon.nelem() != 4,
        "You have to specify limiting values for the "
        "convergence test for each Stokes component "
        "separately. That means that *epsilon* must "
        "have *4* elements!");

  // Check if cloudbox_field and cloudbox_field_old have the same dimensions:
  ARTS_USER_ERROR_IF (!is_size(
          cloudbox_field_mono_old, N_p, N_lat, N_lon, N_za, N_aa, 4),
        "The fields (Tensor6) *cloudbox_field* and \n"
        "*cloudbox_field_old* which are compared in the \n"
        "convergence test do not have the same size.\n");

  // Frequency grid
  //
  ARTS_USER_ERROR_IF (f_grid.empty(), "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  // Is the frequency index valid?
  ARTS_USER_ERROR_IF (f_index >= f_grid.nelem(),
        "*f_index* is greater than number of elements in the\n"
        "frequency grid.\n");

  //-----------End of checks--------------------------------

  doit_iteration_counter += 1;

  //Numeric max_diff_bt = 0.;
  if (doit_iteration_counter > max_iterations) {
    ostringstream out;
    out << "At frequency " << f_grid[f_index] << " GHz \n"
        << "method does not converge (number of iterations \n"
        << "is > " << max_iterations << "). Either the particle"
        << " number density \n"
        << "is too large or the numerical setup for the DOIT \n"
        << "calculation is not correct. In case of limb \n"
        << "simulations please make sure that you use an \n"
        << "optimized zenith angle grid. \n"
        << "*cloudbox_field* might be wrong.\n";
    if (throw_nonconv_error != 0) {
      // FIXME: OLE: Remove this later
      //          ostringstream os;
      //          os << "Error in DOIT calculation:\n"
      //             << out.str();
      //          throw runtime_error( os.str() );
      cloudbox_field_mono = NAN;
      doit_conv_flag = 1;
    } else {
      doit_conv_flag = 1;
    }
  } else {
    for (Index p_index = 0; p_index < N_p; p_index++) {
      for (Index lat_index = 0; lat_index < N_lat; lat_index++) {
        for (Index lon_index = 0; lon_index < N_lon; lon_index++) {
          for (Index za_index = 0; za_index < N_za; za_index++) {
            for (Index aa_index = 0; aa_index < N_aa; aa_index++) {
              for (Index stokes_index = 0; stokes_index < 4;
                   stokes_index++) {
                Numeric diff = cloudbox_field_mono(p_index,
                                                   lat_index,
                                                   lon_index,
                                                   za_index,
                                                   aa_index,
                                                   stokes_index) -
                               cloudbox_field_mono_old(p_index,
                                                       lat_index,
                                                       lon_index,
                                                       za_index,
                                                       aa_index,
                                                       stokes_index);

                // If the absolute difference of the components
                // is larger than the pre-defined values, return
                // to *cloudbox_fieldIterate* and do next iteration
                Numeric diff_bt = invrayjean(diff, f_grid[f_index]);
                //                        if( abs(diff_bt) > max_diff_bt )
                //                          max_diff_bt = abs(diff_bt);
                if (abs(diff_bt) > epsilon[stokes_index]) {
                  //                            cout << "max BT difference in iteration #"
                  //                                 << doit_iteration_counter << ": "
                  //                                 << max_diff_bt << "\n";
                  return;
                }
              }  // End loop stokes_dom.
            }    // End loop aa_grid.
          }      // End loop za_grid.
        }        // End loop lon_grid.
      }          // End loop lat_grid.
    }            // End p_grid.

    //cout << "max BT difference in iteration #" << doit_iteration_counter
    //     << ": " << max_diff_bt << "\n";
    // Convergence test has been successful, doit_conv_flag can be set to 1.
    doit_conv_flag = 1;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagLsq(  //WS Output:
    Index& doit_conv_flag,
    Index& doit_iteration_counter,
    Tensor6& cloudbox_field_mono,
    // WS Input:
    const Tensor6& cloudbox_field_mono_old,
    const Vector& f_grid,
    const Index& f_index,
    // Keyword:
    const Vector& epsilon,
    const Index& max_iterations,
    const Index& throw_nonconv_error) {
  //------------Check the input---------------------------------------

  ARTS_USER_ERROR_IF (doit_conv_flag != 0,
        "Convergence flag is non-zero, which means that this \n"
        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
        "be used only in *doit_conv_test_agenda*\n");

  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_lat = cloudbox_field_mono.nshelves();
  const Index N_lon = cloudbox_field_mono.nbooks();
  const Index N_za = cloudbox_field_mono.npages();
  const Index N_aa = cloudbox_field_mono.nrows();

  // Check keyword "epsilon":
  ARTS_USER_ERROR_IF (epsilon.nelem() != 4,
        "You have to specify limiting values for the "
        "convergence test for each Stokes component "
        "separately. That means that *epsilon* must "
        "have *4* elements!");

  // Check if cloudbox_field and cloudbox_field_old have the same dimensions:
  ARTS_USER_ERROR_IF (!is_size(
          cloudbox_field_mono_old, N_p, N_lat, N_lon, N_za, N_aa, 4),
        "The fields (Tensor6) *cloudbox_field* and \n"
        "*cloudbox_field_old* which are compared in the \n"
        "convergence test do not have the same size.\n");

  // Frequency grid
  //
  ARTS_USER_ERROR_IF (f_grid.empty(), "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  // Is the frequency index valid?
  ARTS_USER_ERROR_IF (f_index >= f_grid.nelem(),
        "*f_index* is greater than number of elements in the\n"
        "frequency grid.\n");

  //-----------End of checks--------------------------------

  doit_iteration_counter += 1;

  if (doit_iteration_counter > max_iterations) {
    ostringstream out;
    out << "Method does not converge (number of iterations \n"
        << "is > " << max_iterations << "). Either the"
        << " particle number density \n"
        << "is too large or the numerical setup for the DOIT \n"
        << "calculation is not correct. In case of limb \n"
        << "simulations please make sure that you use an \n"
        << "optimized zenith angle grid. \n";
    if (throw_nonconv_error != 0) {
      // FIXME: OLE: Remove this later
      //          ostringstream os;
      //          os << "Error in DOIT calculation:\n"
      //             << out.str();
      //          throw runtime_error( os.str() );
      cloudbox_field_mono = NAN;
      doit_conv_flag = 1;
    } else {
      doit_conv_flag = 1;
    }
  } else {
    Vector lqs(4, 0.);

    // Will be set to zero if convergence not fullfilled
    doit_conv_flag = 1;
    for (Index i = 0; i < epsilon.nelem(); i++) {
      for (Index p_index = 0; p_index < N_p; p_index++) {
        for (Index lat_index = 0; lat_index < N_lat; lat_index++) {
          for (Index lon_index = 0; lon_index < N_lon; lon_index++) {
            for (Index za_index = 0; za_index < N_za; za_index++) {
              for (Index aa_index = 0; aa_index < N_aa; aa_index++) {
                lqs[i] += pow(
                    cloudbox_field_mono(
                        p_index, lat_index, lon_index, za_index, aa_index, i) -
                        cloudbox_field_mono_old(p_index,
                                                lat_index,
                                                lon_index,
                                                za_index,
                                                aa_index,
                                                i),
                    2);
              }  // End loop aa_grid.
            }    // End loop za_grid.
          }      // End loop lon_grid.
        }        // End loop lat_grid.
      }          // End p_grid.

      lqs[i] = sqrt(lqs[i]);
      lqs[i] /= (Numeric)(N_p * N_lat * N_lon * N_za * N_aa);

      // Convert difference to Rayleigh Jeans BT
      lqs[i] = invrayjean(lqs[i], f_grid[f_index]);

      if (lqs[i] >= epsilon[i]) doit_conv_flag = 0;
    }
    // end loop stokes_index
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_field_monoIterate(const Workspace& ws,
                                // WS Input and Output:
                                Tensor6& cloudbox_field_mono,

                                // WS Input:
                                const Agenda& doit_scat_field_agenda,
                                const Agenda& doit_rte_agenda,
                                const Agenda& doit_conv_test_agenda,
                                const Index& accelerated)

{
  for (Index v = 0; v < cloudbox_field_mono.nvitrines(); v++)
    for (Index s = 0; s < cloudbox_field_mono.nshelves(); s++)
      for (Index b = 0; b < cloudbox_field_mono.nbooks(); b++)
        for (Index p = 0; p < cloudbox_field_mono.npages(); p++)
          for (Index r = 0; r < cloudbox_field_mono.nrows(); r++)
            for (Index c = 0; c < cloudbox_field_mono.ncols(); c++)
              ARTS_USER_ERROR_IF (std::isnan(cloudbox_field_mono(v, s, b, p, r, c)),
                    "*cloudbox_field_mono* contains at least one NaN value.\n"
                    "This indicates an improper initialization of *cloudbox_field*.");

  //cloudbox_field_mono can not be further checked here, because there is no way
  //to find out the size without including a lot more interface
  //variables
  //-----------End of checks--------------------------------------

  Tensor6 cloudbox_field_mono_old_local;
  Index doit_conv_flag_local;
  Index doit_iteration_counter_local;

  // Resize and initialize doit_scat_field,
  // which  has the same dimensions as cloudbox_field
  Tensor6 doit_scat_field_local(cloudbox_field_mono.nvitrines(),
                                cloudbox_field_mono.nshelves(),
                                cloudbox_field_mono.nbooks(),
                                cloudbox_field_mono.npages(),
                                cloudbox_field_mono.nrows(),
                                cloudbox_field_mono.ncols(),
                                0.);

  doit_conv_flag_local = 0;
  doit_iteration_counter_local = 0;
  // Array to save the last iteration steps
  ArrayOfTensor6 acceleration_input;
  if (accelerated) {
    acceleration_input.resize(4);
  }
  while (doit_conv_flag_local == 0) {
    // 1. Copy cloudbox_field to cloudbox_field_old.
    cloudbox_field_mono_old_local = cloudbox_field_mono;

    // 2.Calculate scattered field vector for all points in the cloudbox.

    // Calculate the scattered field.
    doit_scat_field_agendaExecute(
        ws, doit_scat_field_local, cloudbox_field_mono, doit_scat_field_agenda);

    // Update cloudbox_field.
    doit_rte_agendaExecute(
        ws, cloudbox_field_mono, doit_scat_field_local, doit_rte_agenda);

    //Convergence test.
    doit_conv_test_agendaExecute(ws,
                                 doit_conv_flag_local,
                                 doit_iteration_counter_local,
                                 cloudbox_field_mono,
                                 cloudbox_field_mono_old_local,
                                 doit_conv_test_agenda);

    // Convergence Acceleration, if wished.
    if (accelerated > 0 && doit_conv_flag_local == 0) {
      acceleration_input[(doit_iteration_counter_local - 1) % 4] =
          cloudbox_field_mono;
      // NG - Acceleration
      if (doit_iteration_counter_local % 4 == 0) {
        cloudbox_field_ngAcceleration(
            cloudbox_field_mono, acceleration_input, accelerated);
      }
    }
  }  //end of while loop, convergence is reached.
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DoitInit(  //WS Output
    Tensor6& doit_scat_field,
    Tensor7& cloudbox_field,
    Index& doit_is_initialized,
    // WS Input
    
    const Vector& f_grid,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& doit_za_grid_size,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits) {
  if (!cloudbox_on) {
    return;
    //throw runtime_error( "Cloudbox is off, no scattering calculations to be"
    //                     "performed." );
  }

  // -------------- Check the input ------------------------------

  // Number of zenith angles.
  const Index N_scat_za = za_grid.nelem();

  // The recommended values were found by testing the accuracy and the speed of
  // 1D DOIT calculations for different grid sizes. For 3D calculations it can
  // be necessary to use more grid points.
  ARTS_USER_ERROR_IF (N_scat_za < 16,
        "For accurate results, za_grid must have "
        "more than 15 elements.");
  if (N_scat_za > 100) {
  }

  ARTS_USER_ERROR_IF (za_grid[0] != 0. || za_grid[N_scat_za - 1] != 180.,
                      "The range of *za_grid* must [0 180].");

  ARTS_USER_ERROR_IF (!is_increasing(za_grid),
                      "*za_grid* must be increasing.");

  // Number of azimuth angles.
  const Index N_scat_aa = aa_grid.nelem();

  ARTS_USER_ERROR_IF (N_scat_aa < 6,
        "For accurate results, aa_grid must have "
        "more than 5 elements.");
  if (N_scat_aa > 100) {
  }

  ARTS_USER_ERROR_IF (aa_grid[0] != 0. || aa_grid[N_scat_aa - 1] != 360.,
                      "The range of *aa_grid* must [0 360].");

  ARTS_USER_ERROR_IF (doit_za_grid_size < 16,
        "*doit_za_grid_size* must be greater than 15 for accurate results");
  if (doit_za_grid_size > 100) {
  }

  ARTS_USER_ERROR_IF (cloudbox_limits.nelem() != 6,
        "*cloudbox_limits* is a vector which contains the"
        "upper and lower limit of the cloud for all "
        "atmospheric dimensions. So its dimension must"
        "be 2 x *3*");

  //------------- end of checks ---------------------------------------

  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index Nza = za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  const Index Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  const Index Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  const Index Naa = aa_grid.nelem();

  cloudbox_field.resize(Nf, Np_cloud, Nlat_cloud, Nlon_cloud, Nza, Naa, 4);
  doit_scat_field.resize(Np_cloud, Nlat_cloud, Nlon_cloud, Nza, Naa, 4);

  cloudbox_field = NAN;
  doit_scat_field = NAN;
  doit_is_initialized = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_field_monoOptimizeReverse(  //WS input
    Tensor6& cloudbox_field_mono,
    const Vector& p_grid_orig,
    const Vector& p_grid,
    const ArrayOfIndex& cloudbox_limits) {
  Tensor6 cloudbox_field_mono_opt(
      p_grid_orig.nelem(), 1, 1, cloudbox_field_mono.npages(), 1, 1);
  ArrayOfGridPos Z_gp(p_grid_orig.nelem());
  Matrix itw_z(Z_gp.nelem(), 2);
  // We only need the p_grid inside the cloudbox as cloudbox_field_mono is only defined in the
  // cloudbox
  Vector p_grid_cloudbox {p_grid[Range(
      cloudbox_limits[0], cloudbox_limits[1] - cloudbox_limits[0] + 1)]};
  ostringstream os;
  os << "There is a problem with the pressure grid interpolation";
  chk_interpolation_grids(os.str(), p_grid, p_grid_orig);

  // Gridpositions:
  gridpos(Z_gp, p_grid_cloudbox, p_grid_orig);
  // Interpolation weights:
  interpweights(itw_z, Z_gp);
  // Interpolate cloudbox_field_mono
  for (Index i = 0; i < cloudbox_field_mono.npages(); i++) {
    interp(cloudbox_field_mono_opt(joker, 0, 0, i, 0, 0),
           itw_z,
           cloudbox_field_mono(joker, 0, 0, i, 0, 0),
           Z_gp);
  }
  cloudbox_field_mono = cloudbox_field_mono_opt;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void OptimizeDoitPressureGrid(
    const Workspace& ws,
    //WS input
    Vector& p_grid,
    Tensor4& pnd_field,
    Tensor3& t_field,
    ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    Tensor3& z_field,
    ArrayOfIndex& cloudbox_limits,
    Tensor6& cloudbox_field_mono,
    Tensor7& pha_mat_doit,
    Tensor4& vmr_field,
    Vector& p_grid_orig,
    const Vector& f_grid,
    const Index& f_index,
    const Agenda& propmat_clearsky_agenda,
    const Numeric& tau_scat_max,
    const Numeric& sgl_alb_max,
    const Index& cloudbox_size_max) {
  // Make sure that 4 = 1 and that ScatSpeciesMerge has been applied:
  ARTS_USER_ERROR_IF (cloudbox_field_mono.ncols() != 1,
        " This method currently only works for unpolarized radiation "
        "( 4 = 1)");
  // If ScatSpeciesMerged has been applied, then scat_data_mono should have only one element, and
  // pnd_field should have the number of pressure grid points in the cloudbox as first dimension
  ARTS_USER_ERROR_IF (scat_data_mono.nelem() > 1 ||
      pnd_field.nbooks() != cloudbox_limits[1] - cloudbox_limits[0] + 1,
        " ScatSpeciesMerge has to be applied in order to use this method");

  bool was_too_much = false;
  Numeric tau_scat_max_internal = tau_scat_max;
  ArrayOfSingleScatteringData& scat_data_local = scat_data_mono[0];
  p_grid_orig = p_grid;
  vector<Numeric> z_grid_new;
  Vector ext_mat(cloudbox_limits[1] - cloudbox_limits[0] + 1);
  Vector abs_vec(cloudbox_limits[1] - cloudbox_limits[0] + 1);
  Vector scat_vec(cloudbox_limits[1] - cloudbox_limits[0] + 1);
  ArrayOfIndex cloudbox_limits_opt(2);
  z_grid_new.reserve(1000);
  for (Index i = 0; i < cloudbox_limits[0]; i++)
    z_grid_new.push_back(z_field(i, 0, 0));
  //-----------------------------------------------
  // Calculate optical thicknesses of the layers
  //------------------------------------------------

  // Fields for scalar gas absorption
  const Vector rtp_mag_dummy(3, 0);
  const Vector ppath_los_dummy;
  StokvecVector nlte_dummy;
  PropmatMatrix partial_dummy;
  StokvecMatrix partial_nlte_dummy;
  PropmatVector cur_propmat_clearsky;

  Index scat_data_insert_offset = 0;
  Vector single_scattering_albedo(cloudbox_limits[1] - cloudbox_limits[0] + 1,
                                  0.);
  for (Index k = cloudbox_limits[0]; k < cloudbox_limits[1] + 1; ++k) {
    Index cloudbox_index = k - cloudbox_limits[0];
    Numeric abs_coeff = 0;
    // Scattering particles
    ext_mat[cloudbox_index] =
        scat_data_local[cloudbox_index].ext_mat_data(0, 0, 0, 0, 0);
    abs_vec[cloudbox_index] =
        scat_data_local[cloudbox_index].abs_vec_data(0, 0, 0, 0, 0);
    scat_vec[cloudbox_index] =
        ext_mat[cloudbox_index] - abs_vec[cloudbox_index];
    // Calculate scalar gas absorption
    propmat_clearsky_agendaExecute(ws,
                                   cur_propmat_clearsky,
                                   nlte_dummy,
                                   partial_dummy,
                                   partial_nlte_dummy,
                                   ArrayOfRetrievalQuantity(0),
                                   {},
                                   Vector{f_grid[Range(f_index, 1)]},
                                   ppath_los_dummy,
                                   AtmPoint{},  // FIXME: DUMMY VALUE,
                                   propmat_clearsky_agenda);
    abs_coeff += cur_propmat_clearsky[0].A();
    abs_coeff /= (Numeric)vmr_field.nbooks();  // FIXME: Is this really as intended???
    single_scattering_albedo[cloudbox_index] =
        scat_vec[cloudbox_index] / (ext_mat[cloudbox_index] + abs_coeff);
  }
  //
  // Try out the current settings of tau_scat_max and sgl_alb_max
  //
  do {
    scat_data_insert_offset = 0;
    for (Index k = cloudbox_limits[0]; k < cloudbox_limits[1]; ++k) {
      Index cloudbox_index = k - cloudbox_limits[0];
      const Numeric sgl_alb = (single_scattering_albedo[cloudbox_index] +
                               single_scattering_albedo[cloudbox_index + 1]) /
                              2;
      const Numeric scat_opt_thk =
          (z_field(k + 1, 0, 0) - z_field(k, 0, 0)) *
          ((scat_vec[cloudbox_index] + scat_vec[cloudbox_index + 1]) / 2);
      if (scat_opt_thk > tau_scat_max_internal && sgl_alb > sgl_alb_max) {
        Index factor = (Index)ceil(scat_opt_thk / tau_scat_max_internal);
        for (Index j = 1; j < factor; j++) {
          scat_data_insert_offset++;
        }
      }
    }
    // If enhancement is too large, change tau_scat_max to a higher value:
    if (scat_data_insert_offset + cloudbox_limits[1] - cloudbox_limits[0] + 1 >
        cloudbox_size_max) {
      tau_scat_max_internal += 0.01;
      was_too_much = true;
    }
  } while (scat_data_insert_offset + cloudbox_limits[1] - cloudbox_limits[0] +
               1 >
           cloudbox_size_max);
  scat_data_insert_offset = 0;
  // Give warning if enhancement was too much and threshold had to be changed:
  if (was_too_much) {
  }

  //--------------------------------------
  //Optimize the altitude grid z_grid
  //---------------------------------------
  for (Index k = cloudbox_limits[0]; k < cloudbox_limits[1]; ++k) {
    Index cloudbox_index = k - cloudbox_limits[0];
    const Numeric sgl_alb = (single_scattering_albedo[cloudbox_index] +
                             single_scattering_albedo[cloudbox_index + 1]) /
                            2;
    const Numeric scat_opt_thk =
        (z_field(k + 1, 0, 0) - z_field(k, 0, 0)) *
        ((scat_vec[cloudbox_index] + scat_vec[cloudbox_index + 1]) / 2);
    z_grid_new.push_back(z_field(k, 0, 0));

    if (scat_opt_thk > tau_scat_max_internal && sgl_alb > sgl_alb_max) {
      Index factor = (Index)ceil(scat_opt_thk / tau_scat_max_internal);
      Numeric step =
          (z_field(k + 1, 0, 0) - z_field(k, 0, 0)) / (Numeric)factor;
      const SingleScatteringData nextLayer =
          scat_data_local[cloudbox_index + scat_data_insert_offset + 1];
      const SingleScatteringData currentLayer =
          scat_data_local[cloudbox_index + scat_data_insert_offset];

      for (Index j = 1; j < factor; j++) {
        z_grid_new.push_back(z_field(k, 0, 0) + (Numeric)j * step);
        // Perform manual interpolation of scat_data
        const Numeric weight = (Numeric)j / (Numeric)factor;
        SingleScatteringData newLayer = currentLayer;
        Tensor7 weightednextLayerPhamat = nextLayer.pha_mat_data;
        Tensor5 weightednextLayerExtmat = nextLayer.ext_mat_data;
        Tensor5 weightednextLayerAbsvec = nextLayer.abs_vec_data;

        weightednextLayerPhamat *= weight;
        weightednextLayerExtmat *= weight;
        weightednextLayerAbsvec *= weight;

        newLayer.pha_mat_data *= 1. - weight;
        newLayer.ext_mat_data *= 1. - weight;
        newLayer.abs_vec_data *= 1. - weight;

        newLayer.pha_mat_data += weightednextLayerPhamat;
        newLayer.ext_mat_data += weightednextLayerExtmat;
        newLayer.abs_vec_data += weightednextLayerAbsvec;

        // Optimize scat_data
        scat_data_local.insert(scat_data_local.begin() + cloudbox_index +
                                   scat_data_insert_offset + 1,
                               std::move(newLayer));

        scat_data_insert_offset++;
      }
    }
  }
  // New cloudbox limits
  cloudbox_limits_opt[0] = cloudbox_limits[0];
  cloudbox_limits_opt[1] = scat_data_insert_offset + cloudbox_limits[1];
  const Index cloudbox_opt_size =
      cloudbox_limits_opt[1] - cloudbox_limits_opt[0] + 1;

  for (Index i = cloudbox_limits[1]; i < z_field.npages(); i++)
    z_grid_new.push_back(z_field(i, 0, 0));

  Vector z_grid(z_grid_new.size());
  for (Index i = 0; i < z_grid.nelem(); i++) z_grid[i] = z_grid_new[i];
  p_grid_orig = p_grid[Range(cloudbox_limits[0],
                             cloudbox_limits[1] - cloudbox_limits[0] + 1)];
  // ---------------------------------------
  // Interpolate fields to new z_grid
  // ----------------------------------------
  ArrayOfArrayOfSingleScatteringData scat_data_new;
  Tensor3 t_field_new(z_grid.nelem(), 1, 1);
  Vector p_grid_opt(z_grid.nelem());
  Tensor6 cloudbox_field_mono_opt(
      cloudbox_opt_size, 1, 1, cloudbox_field_mono.npages(), 1, 1);
  Tensor7 pha_mat_doit_opt(cloudbox_opt_size,
                           pha_mat_doit.nvitrines(),
                           1,
                           pha_mat_doit.nbooks(),
                           pha_mat_doit.npages(),
                           1,
                           1);
  ArrayOfGridPos Z_gp(z_grid.nelem());
  Matrix itw_z(Z_gp.nelem(), 2);
  ostringstream os;
  os << "At the current frequency " << f_grid[f_index]
     << " there was an error while interpolating the fields to the new z_field";
  chk_interpolation_grids(os.str(), z_field(joker, 0, 0), z_grid);

  // Gridpositions of interpolation:
  gridpos(Z_gp, z_field(joker, 0, 0), z_grid);
  // Interpolation weights:
  interpweights(itw_z, Z_gp);

  // Interpolate Temperature
  interp(t_field_new(joker, 0, 0), itw_z, t_field(joker, 0, 0), Z_gp);
  // Write new Temperature to scat_data
  for (Index k = cloudbox_limits_opt[0]; k < cloudbox_limits_opt[1]; k++) {
    Index i = k - cloudbox_limits[0];
    scat_data_local[i].T_grid = t_field_new(i, 0, 0);
  }

  // Interpolate p_grid
  interp(p_grid_opt, itw_z, p_grid, Z_gp);

  // Interpolate vmr_field
  Tensor4 vmr_field_opt(vmr_field.nbooks(), p_grid_opt.nelem(), 1, 1);
  for (Index i = 0; i < vmr_field.nbooks(); i++)
    interp(
        vmr_field_opt(i, joker, 0, 0), itw_z, vmr_field(i, joker, 0, 0), Z_gp);

  // Interpolate cloudbox_field_mono and pha_mat_doit
  ArrayOfGridPos Z_gp_2(cloudbox_opt_size);
  Matrix itw_z_2(Z_gp_2.nelem(), 2);
  Range r1 =
      Range(cloudbox_limits[0], cloudbox_limits[1] - cloudbox_limits[0] + 1);
  Range r2 = Range(cloudbox_limits_opt[0], cloudbox_opt_size);
  chk_interpolation_grids(os.str(), z_field(r1, 0, 0), z_grid[r2]);
  gridpos(Z_gp_2,
          z_field(Range(cloudbox_limits[0],
                        cloudbox_limits[1] - cloudbox_limits[0] + 1),
                  0,
                  0),
          z_grid[Range(cloudbox_limits_opt[0], cloudbox_opt_size)]);
  interpweights(itw_z_2, Z_gp_2);

  for (Index i = 0; i < cloudbox_field_mono.npages(); i++) {
    interp(cloudbox_field_mono_opt(joker, 0, 0, i, 0, 0),
           itw_z_2,
           cloudbox_field_mono(joker, 0, 0, i, 0, 0),
           Z_gp_2);
  }
  for (Index i = 0; i < pha_mat_doit.nvitrines(); i++) {
    for (Index j = 0; j < pha_mat_doit.nbooks(); j++) {
      for (Index k = 0; k < pha_mat_doit.npages(); k++) {
        interp(pha_mat_doit_opt(joker, i, 0, j, k, 0, 0),
               itw_z_2,
               pha_mat_doit(joker, i, 0, j, k, 0, 0),
               Z_gp_2);
      }
    }
  }

  // Interpolate pnd-field
  pnd_field.resize(cloudbox_opt_size, cloudbox_opt_size, 1, 1);
  pnd_field = 0.;
  for (Index i = 0; i < cloudbox_opt_size; i++) pnd_field(i, i, 0, 0) = 1.;

  //Write new fields
  p_grid = p_grid_opt;
  t_field = t_field_new;
  cloudbox_limits = cloudbox_limits_opt;
  cloudbox_field_mono = cloudbox_field_mono_opt;
  pha_mat_doit = pha_mat_doit_opt;
  z_field.resize(z_grid.nelem(), 1, 1);
  z_field(joker, 0, 0) = z_grid;
  vmr_field = vmr_field_opt;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DoitWriteIterationFields(  //WS input
    const Index& doit_iteration_counter,
    const Tensor6& cloudbox_field_mono,
    const Index& f_index,
    //Keyword:
    const ArrayOfIndex& iterations,
    const ArrayOfIndex& frequencies) {
  if (!frequencies.nelem() || !iterations.nelem()) return;

  // If the number of iterations is less than a number specified in the
  // keyword *iterations*, this number will be ignored.

  ostringstream os;
  os << "doit_iteration_f" << f_index << "_i" << doit_iteration_counter;

  // All iterations for all frequencies are written to files
  if (frequencies[0] == -1 && iterations[0] == -1) {
    xml_write_to_file(
        os.str() + ".xml", cloudbox_field_mono, FILE_TYPE_ASCII, 0);
  }

  for (Index j = 0; j < frequencies.nelem(); j++) {
    if (f_index == frequencies[j] || (!j && frequencies[j] == -1)) {
      // All iterations are written to files
      if (iterations[0] == -1) {
        xml_write_to_file(os.str() + ".xml",
                          cloudbox_field_mono,
                          FILE_TYPE_ASCII,
                          0);
      }

      // Only the iterations given by the keyword are written to a file
      else {
        for (Index i = 0; i < iterations.nelem(); i++) {
          if (doit_iteration_counter == iterations[i])
            xml_write_to_file(os.str() + ".xml",
                              cloudbox_field_mono,
                              FILE_TYPE_ASCII,
                              0);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_za_grid_optCalc(  //WS Output
    Vector& doit_za_grid_opt,
    // WS Input:
    const Tensor6& cloudbox_field_mono,
    const Vector& za_grid,
    const Index& doit_za_interp,
    //Keywords:
    const Numeric& acc) {
  //-------- Check the input ---------------------------------

  // Here it is checked whether cloudbox_field is 1D and whether it is
  // consistent with za_grid. The number of pressure levels and the
  // number of stokes components does not matter.
  chk_size("cloudbox_field",
           cloudbox_field_mono,
           cloudbox_field_mono.nvitrines(),
           1,
           1,
           za_grid.nelem(),
           1,
           cloudbox_field_mono.ncols());

  ARTS_USER_ERROR_IF (cloudbox_field_mono.ncols() < 1 || cloudbox_field_mono.ncols() > 4,
        "The last dimension of *cloudbox_field* corresponds\n"
        "to the Stokes dimension, therefore the number of\n"
        "columns in *cloudbox_field* must be a number between\n"
        "1 and 4, but it is not!");

  ARTS_USER_ERROR_IF (!(doit_za_interp == 0 || doit_za_interp == 1),
        "Interpolation method is not defined. Use \n"
        "*doit_za_interpSet*.\n");

  if (za_grid.nelem() < 500) {
    /*    throw runtime_error("The fine grid (*za_grid*) has less than \n"
                        "500 grid points which is not sufficient for \n"
                        "grid_optimization");
*/
  }
  // ------------- end of checks -------------------------------------

  // Here only used as dummy variable.
  Matrix cloudbox_field_opt_mat;
  cloudbox_field_opt_mat = 0.;

  // Optimize zenith angle grid.
  za_gridOpt(doit_za_grid_opt,
             cloudbox_field_opt_mat,
             za_grid,
             cloudbox_field_mono,
             acc,
             doit_za_interp);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_za_interpSet(Index& doit_za_interp,
                       //Keyword
                       const String& method) {
  ARTS_USER_ERROR_IF (3 != 1 && method == "polynomial",
        "Polynomial interpolation is only implemented for\n"
        "1D DOIT calculations as \n"
        "in 3D there can be numerical problems.\n"
        "Please use 'linear' interpolation method.");

  if (method == "linear")
    doit_za_interp = 0;
  else if (method == "polynomial")
    doit_za_interp = 1;
  else {
    ARTS_USER_ERROR (
        "Possible interpolation methods are 'linear' "
        "and 'polynomial'.\n");
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DoitCalc(const Workspace& ws,
              Tensor7& cloudbox_field,
              const Index& atmfields_checked,
              const Index& atmgeom_checked,
              const Index& cloudbox_checked,
              const Index& scat_data_checked,
              const Index& cloudbox_on,
              const Vector& f_grid,
              const Agenda& doit_mono_agenda,
              const Index& doit_is_initialized)

{
  if (!cloudbox_on) {
    return;
    //throw runtime_error( "Cloudbox is off, no scattering calculations to be"
    //                     "performed." );
  }

  //-------- Check input -------------------------------------------

  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;

  ARTS_USER_ERROR_IF (scat_data_checked != 1,
        "The scattering data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  // Frequency grid
  //
  ARTS_USER_ERROR_IF (f_grid.empty(), "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  // Check whether DoitInit was executed
  ARTS_USER_ERROR_IF (!doit_is_initialized,
    "Initialization method *DoitInit* has to be called before "
    "*DoitCalc*")

  //-------- end of checks ----------------------------------------

  // OMP likes simple loop end conditions, so we make a local copy here:
  const Index nf = f_grid.nelem();

  if (nf) {
    String fail_msg;
    bool failed = false;

#pragma omp parallel for if (!arts_omp_in_parallel() && nf > 1)
    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) {
        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        continue;
      }

      try {
        ostringstream os;
        os << "Frequency: " << f_grid[f_index] / 1e9 << " GHz \n";

        Tensor6 cloudbox_field_mono_local{
            cloudbox_field(f_index, joker, joker, joker, joker, joker, joker)};
        doit_mono_agendaExecute(ws,
                                cloudbox_field_mono_local,
                                f_grid,
                                f_index,
                                doit_mono_agenda);
        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) =
            cloudbox_field_mono_local;
      } catch (const std::exception& e) {
        cloudbox_field(f_index, joker, joker, joker, joker, joker, joker) = NAN;
        ostringstream os;
        os << "Error for f_index = " << f_index << " (" << f_grid[f_index]
           << " Hz)" << endl
           << e.what();
#pragma omp critical(DoitCalc_fail)
        {
          failed = true;
          fail_msg = os.str();
        }
        continue;
      }
    }

    ARTS_USER_ERROR_IF (failed, fail_msg);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_fieldSetClearsky(Tensor7& cloudbox_field,
                               const Vector& p_grid,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const ArrayOfIndex& cloudbox_limits,
                               const Index& cloudbox_on,
                               const Index& doit_is_initialized,
                               const Index& all_frequencies) {
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) return;

  // Check whether DoitInit was executed
  ARTS_USER_ERROR_IF (!doit_is_initialized,
      "Initialization method *DoitInit* has to be called before "
      "*cloudbox_fieldSetClearsky*.")

  // Initial field only needs to be calculated from clearsky field for the
  // first frequency. For the next frequencies the solution field from the
  // previous frequencies is used.
  {
    ARTS_USER_ERROR_IF (all_frequencies == false,
          "Error in cloudbox_fieldSetClearsky: For 3D "
          "all_frequencies option is not implemented \n");

    for (Index f_index = 0; f_index < cloudbox_field.nvitrines(); f_index++) {
      Index N_p = cloudbox_field.nvitrines();
      Index N_lat = cloudbox_field.nshelves();
      Index N_lon = cloudbox_field.nbooks();
      Index N_za = cloudbox_field.npages();
      Index N_aa = cloudbox_field.nrows();
      Index N_i = cloudbox_field.ncols();

      Tensor6 scat_i_p(2, N_lat, N_lon, N_za, N_aa, N_i);
      scat_i_p(0, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index, 0, joker, joker, joker, joker, joker);
      scat_i_p(1, joker, joker, joker, joker, joker) =
          cloudbox_field(f_index, N_p - 1, joker, joker, joker, joker, joker);

      Tensor6 scat_i_lat(N_p, 2, N_lon, N_za, N_aa, N_i);
      scat_i_lat(joker, 0, joker, joker, joker, joker) =
          cloudbox_field(f_index, joker, 0, joker, joker, joker, joker);
      scat_i_lat(joker, 1, joker, joker, joker, joker) =
          cloudbox_field(f_index, joker, N_lat - 1, joker, joker, joker, joker);

      Tensor6 scat_i_lon(N_p, N_lat, 2, N_za, N_aa, N_i);
      scat_i_lon(joker, joker, 0, joker, joker, joker) =
          cloudbox_field(f_index, joker, joker, 0, joker, joker, joker);
      scat_i_lon(joker, joker, 1, joker, joker, joker) =
          cloudbox_field(f_index, joker, joker, N_lon - 1, joker, joker, joker);

      //1. interpolation - pressure grid, latitude grid and longitude grid

      ArrayOfGridPos p_gp((cloudbox_limits[1] - cloudbox_limits[0]) + 1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3] - cloudbox_limits[2]) + 1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5] - cloudbox_limits[4]) + 1);

      /*the old grid is having only two elements, corresponding to the
             cloudbox_limits and the new grid have elements corresponding to
             all grid points inside the cloudbox plus the cloud_box_limits*/

      p2gridpos(p_gp,
                p_grid[Range(cloudbox_limits[0],
                             2,
                             (cloudbox_limits[1] - cloudbox_limits[0]))],
                p_grid[Range(cloudbox_limits[0],
                             (cloudbox_limits[1] - cloudbox_limits[0]) + 1)]);
      gridpos(lat_gp,
              lat_grid[Range(cloudbox_limits[2],
                             2,
                             (cloudbox_limits[3] - cloudbox_limits[2]))],
              lat_grid[Range(cloudbox_limits[2],
                             (cloudbox_limits[3] - cloudbox_limits[2]) + 1)]);
      gridpos(lon_gp,
              lon_grid[Range(cloudbox_limits[4],
                             2,
                             (cloudbox_limits[5] - cloudbox_limits[4]))],
              lon_grid[Range(cloudbox_limits[4],
                             (cloudbox_limits[5] - cloudbox_limits[4]) + 1)]);

      //interpolation weights corresponding to pressure, latitude and
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1] - cloudbox_limits[0]) + 1, 2);
      Matrix itw_lat((cloudbox_limits[3] - cloudbox_limits[2]) + 1, 2);
      Matrix itw_lon((cloudbox_limits[5] - cloudbox_limits[4]) + 1, 2);

      interpweights(itw_p, p_gp);
      interpweights(itw_lat, lat_gp);
      interpweights(itw_lon, lon_gp);

      // interpolation - pressure grid
      for (Index lat_index = 0;
           lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
           ++lat_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                         Range(joker),
                                                         lat_index,
                                                         lon_index,
                                                         za_index,
                                                         aa_index,
                                                         i);

                ConstVectorView source_field = scat_i_p(
                    Range(joker), lat_index, lon_index, za_index, aa_index, i);

                interp(target_field, itw_p, source_field, p_gp);
              }
            }
          }
        }
      }
      //interpolation latitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lon_index = 0;
             lon_index <= (cloudbox_limits[5] - cloudbox_limits[4]);
             ++lon_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                         p_index,
                                                         Range(joker),
                                                         lon_index,
                                                         za_index,
                                                         aa_index,
                                                         i);

                ConstVectorView source_field = scat_i_lat(
                    p_index, Range(joker), lon_index, za_index, aa_index, i);

                interp(target_field, itw_lat, source_field, lat_gp);
              }
            }
          }
        }
      }
      //interpolation -longitude
      for (Index p_index = 0;
           p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
           ++p_index) {
        for (Index lat_index = 0;
             lat_index <= (cloudbox_limits[3] - cloudbox_limits[2]);
             ++lat_index) {
          for (Index za_index = 0; za_index < N_za; ++za_index) {
            for (Index aa_index = 0; aa_index < N_aa; ++aa_index) {
              for (Index i = 0; i < N_i; ++i) {
                VectorView target_field = cloudbox_field(f_index,
                                                         p_index,
                                                         lat_index,
                                                         Range(joker),
                                                         za_index,
                                                         aa_index,
                                                         i);

                ConstVectorView source_field = scat_i_lon(
                    p_index, lat_index, Range(joker), za_index, aa_index, i);

                interp(target_field, itw_lon, source_field, lon_gp);
              }
            }
          }
        }
      }  //end of interpolation
    }    // end of frequency loop
  }      //ends 3 = 3
}