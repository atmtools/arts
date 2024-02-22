/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_checked.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2013-08-20

  \brief  Workspace functions setting the checked WSVs

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

#include <workspace.h>

#include "arts_conversions.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "matpack_data.h"
#include "surf.h"
#include "wigner_functions.h"

inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);

/* Workspace method: Doxygen documentation will be auto-generated */
void atmfields_checkedCalc(Index& atmfields_checked,
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const AtmField& atm_field) {
  // Consistency between dim, grids and atmospheric fields/surfaces
  ARTS_USER_ERROR_IF(not atm_field.has(AtmKey::p), "No pressure field")
  ARTS_USER_ERROR_IF(not atm_field.has(AtmKey::t), "No temperature field")
  for (auto& spec : abs_species)
    ARTS_USER_ERROR_IF(not atm_field.has(spec.Species()),
                       "No ",
                       std::quoted(var_string(spec)),
                       " field")

    {
      const bool u = atm_field.has(AtmKey::mag_u),
                 v = atm_field.has(AtmKey::mag_v),
                 w = atm_field.has(AtmKey::mag_w);
      ARTS_USER_ERROR_IF(
          (u or v or w) and (not u or not v or not w),
          "If any magnetic field component exist, all three must.\n The "
          "following component exists (1) and does not exist (0):\n",
          "U: ",
          static_cast<int>(u),
          "V: ",
          static_cast<int>(v),
          "W: ",
          static_cast<int>(w))
    }

  {
    const bool u = atm_field.has(AtmKey::wind_u),
               v = atm_field.has(AtmKey::wind_v),
               w = atm_field.has(AtmKey::wind_w);
    ARTS_USER_ERROR_IF(
        (u or v or w) and (not u or not v or not w),
        "If any wind field component exist, all three must.\n The "
        "following component exists (1) and does not exist (0):\n",
        "U: ",
        static_cast<int>(u),
        "V: ",
        static_cast<int>(v),
        "W: ",
        static_cast<int>(w))
  }

  // If here, all OK
  atmfields_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atmgeom_checkedCalc(
    Index& atmgeom_checked,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    const SurfaceField& surface_field,  // Adopt to new version ZZZ
    const Matrix& z_surface,
    const Vector& lat_true,
    const Vector& lon_true,
    const Numeric& max500hpa_gradient) {
  // A repetition from atmfields_checked, but we do this to make the two parts
  // independent (the other option would be to demand atmfields_checkec == 1)
  chk_atm_grids(3, p_grid, lat_grid, lon_grid);

  // *refellipsoid*
  chk_refellipsoid(surface_field.ellipsoid);

  chk_atm_field("z_field", z_field, 3, p_grid, lat_grid, lon_grid);
  chk_atm_surface("z_surface", z_surface, 3, lat_grid, lon_grid);

  // Check that z_field has strictly increasing pages.
  for (Index row = 0; row < z_field.nrows(); row++) {
    for (Index col = 0; col < z_field.ncols(); col++) {
      std::ostringstream os;
      os << "z_field (for latitude nr " << row << " and longitude nr " << col
         << ")";
      chk_if_increasing(os.str(), z_field(joker, row, col));
    }
  }

  // Check that there is no gap between the surface and lowest pressure
  // level
  // (A copy of this code piece is found in z_fieldFromHSE. Make this to an
  // internal function if used in more places.)
  for (Index row = 0; row < z_surface.nrows(); row++) {
    for (Index col = 0; col < z_surface.ncols(); col++) {
      ARTS_USER_ERROR_IF(
          z_surface(row, col) < z_field(0, row, col) ||
              z_surface(row, col) >= z_field(z_field.npages() - 1, row, col),
          "The surface altitude (*z_surface*) cannot be outside\n"
          "of the altitudes in *z_field*.\n"
          "z_surface: ",
          z_surface(row, col),
          "\n"
          "min of z_field: ",
          z_field(0, row, col),
          "\n"
          "max of z_field: ",
          z_field(z_field.npages() - 1, row, col),
          "\n",
          (3 > 1) ? var_string("\nThis was found to be the case for:\n",
                               "latitude ",
                               lat_grid[row])
                  : var_string(),
          (3 > 2) ? var_string("\nlongitude ", lon_grid[col]) : var_string())
    }
  }

  // A rough check of gradients in the 500 hPa level (or closest pressure level)
  // This check is just run in the latitude direction, along some longitudes to
  // save time
  if (3 > 1) {
    // Find 500 hPa
    Index ip = -1;
    Numeric dpmin = 99e99;
    for (Index i = 0; i < p_grid.nelem(); i++) {
      const Numeric dp = abs(p_grid[i] - 500e2);
      if (dp < dpmin) {
        ip = i;
        dpmin = dp;
      }
    }
    Numeric maxgrad = -1;
    for (Index ilon = 0; ilon < max(1, lon_grid.nelem()); ilon += 10) {
      for (Index ilat = 1; ilat < lat_grid.nelem(); ilat++) {
        const Numeric grad =
            abs((z_field(ip, ilat, ilon) - z_field(ip, ilat - 1, ilon)) /
                (lat_grid[ilat] - lat_grid[ilat - 1]));
        if (grad > maxgrad) maxgrad = grad;
      }
    }
    // Scale to change over 100 km
    maxgrad *= 100.0 / 111.0;
    if (maxgrad > max500hpa_gradient) {
      std::ostringstream os;
      os << "A check of the altitude of the " << p_grid[ip] / 100
         << " hPa level has been made.\nThe maximum gradient found matches "
         << maxgrad << " m/100km, that exceeds\nthe set limit of "
         << max500hpa_gradient << " m/100km (by GIN *max500hpa_gradient*).\n"
         << "Please check the smoothness of *z_field*.";
      throw std::runtime_error(os.str());
    }
  }

  // lat/lon true
  if (3 < 3 && (lat_true.nelem() || lon_true.nelem())) {
    if (3 == 1) {
      ARTS_USER_ERROR_IF(lat_true.nelem() != 1,
                         "For 1D, *lat_true* must have length 1.");
      ARTS_USER_ERROR_IF(lon_true.nelem() != 1,
                         "For 1D, *lon_true* must have length 1.");
    } else if (3 == 2) {
      ARTS_USER_ERROR_IF(
          lat_true.nelem() != lat_grid.nelem(),
          "For 2D, *lat_true* must have same length as *lat_grid*.");
      ARTS_USER_ERROR_IF(
          lon_true.nelem() != lat_grid.nelem(),
          "For 2D, *lon_true* must have same length as *lat_grid*.");
    }
    ARTS_USER_ERROR_IF(lon_true.nelem() != lat_true.nelem(),
                       "If *lat_true* is set, also *lon_true* must be "
                       "set (and have the same length).");
    ARTS_USER_ERROR_IF(min(lat_true) < -90 || max(lat_true) > 90,
                       "Values in *lat_true* must be inside [-90,90].");
    ARTS_USER_ERROR_IF(min(lon_true) < -180 || max(lon_true) > 360,
                       "Values in *lon_true* must be inside [-180,360].");
  }

  // If here, all OK
  atmgeom_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_checkedCalc(Index& scat_data_checked,
                           const ArrayOfArrayOfSingleScatteringData& scat_data,
                           const Vector& f_grid,
                           const Numeric& dfrel_threshold,
                           const String& check_level,
                           const Numeric& )//sca_mat_threshold)
// FIXME: when we allow K, a, Z to be on different f and T grids, their use in
// the scatt solvers needs to be reviewed again and adapted to this!
{
  // Prevent the user from producing nonsense by using too much freedom in
  // setting the single freq validity threshold...
  ARTS_USER_ERROR_IF(dfrel_threshold > 0.5,
                     "*dfrel_threshold* too large (max. allowed: 0.5, your's: ",
                     dfrel_threshold,
                     ").")

  // freq range of calc covered?
  ARTS_USER_ERROR_IF(f_grid.empty(), "The frequency grid is empty.");
  if (f_grid.nelem() > 1) chk_if_increasing("f_grid", f_grid);

  Index nf = f_grid.nelem();
  Index N_ss = scat_data.size();
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    Index N_se = scat_data[i_ss].size();
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // For each scattering element (se) check that se's f_grid is either
      // identical to f_grid or contains a single entry only. In the latter
      // case, the f/f_grid-ratio (equv. to a size parameter ratio) not allowed
      // to exceed the dfrel_threshold.
      Index nf_se = scat_data[i_ss][i_se].f_grid.nelem();
      if (nf_se != 1) {
        ARTS_USER_ERROR_IF(nf_se != nf,
                           "*scat_data* must have either one or *f_grid* (=",
                           nf,
                           ") frequency entries,\n"
                           "but scattering element #",
                           i_se,
                           " in scattering species #",
                           i_ss,
                           " has ",
                           nf_se,
                           ".")
        for (Index f = 0; f < nf_se; f++) {
          ARTS_USER_ERROR_IF(
              !is_same_within_epsilon(
                  scat_data[i_ss][i_se].f_grid[f], f_grid[f], 0.5e-9),
              "*scat_data* frequency grid has to be identical to *f_grid*\n"
              "(or contain only a single entry),\n"
              "but scattering element #",
              i_se,
              " in scattering species #",
              i_ss,
              " deviates for f_index ",
              f,
              ".")
        }
      } else {
        ARTS_USER_ERROR_IF((abs(1. - scat_data[i_ss][i_se].f_grid[0] /
                                         f_grid[0]) > dfrel_threshold) or
                               (abs(1. - scat_data[i_ss][i_se].f_grid[0] /
                                             f_grid[nf - 1]) > dfrel_threshold),
                           "Frequency entry (f=",
                           scat_data[i_ss][i_se].f_grid[0],
                           "Hz) of scattering element #",
                           i_se,
                           "\n"
                           "in scattering species #",
                           i_ss,
                           " is too far (>",
                           dfrel_threshold * 1e2,
                           "%) from one or more\n"
                           "of the f_grid limits (fmin=",
                           f_grid[0],
                           "Hz, fmax=",
                           f_grid[nf - 1],
                           "Hz).")
      }

      // check that the freq dimension of sca_mat, ext_mat, and abs_vec is
      // either ssd.f_grid.nelem() or 1 (FIXME: so far, being freq dim !=
      // ssd.f_grid.nelem() switched off, as usage in scatt solvers so far
      // doesn't allow this. see FIXME at start.).
      {
        std::ostringstream bs1, bs2;
        bs1 << "Frequency dimension of ";
        //bs2 << " must be either one or ssd.f_grid.nelem() (=" << nf_se << "),\n"
        bs2 << " must be ssd.f_grid.nelem() (=" << nf_se << "),\n"
            << "but scattering element #" << i_se << " in scattering species #"
            << i_ss << " is ";
        Index nf_sd = scat_data[i_ss][i_se].pha_mat_data.nlibraries();
        ARTS_USER_ERROR_IF(
            nf_sd != nf_se, bs1.str(), "pha_mat_data", bs2.str(), nf_se, ".")
        nf_sd = scat_data[i_ss][i_se].ext_mat_data.nshelves();
        ARTS_USER_ERROR_IF(
            nf_sd != nf_se, bs1.str(), "ext_mat_data", bs2.str(), nf_se, ".")
        nf_sd = scat_data[i_ss][i_se].abs_vec_data.nshelves();
        ARTS_USER_ERROR_IF(
            nf_sd != nf_se, bs1.str(), "abs_vec_data", bs2.str(), nf_se, ".")
      }

      // check that the temp dimension of K and a is ssd.T_grid.nelem(). For Z
      // it might be ssd.T_grid.nelem() or 1.
      {
        std::ostringstream bs1, bs2;
        Index nt_se = scat_data[i_ss][i_se].T_grid.nelem();
        bs1 << "Temperature dimension of ";
        //bs2 << " must be either one or ssd.T_grid.nelem(),\n"
        bs2 << " must be ssd.T_grid.nelem() (=" << nt_se << "),\n"
            << "but for scattering element #" << i_se
            << " in scattering species #" << i_ss << " it is ";
        Index nt_sd = scat_data[i_ss][i_se].pha_mat_data.nvitrines();
        ARTS_USER_ERROR_IF(nt_sd != nt_se and nt_sd != 1,
                           bs1.str(),
                           "pha_mat_data",
                           bs2.str(),
                           nt_sd,
                           ".")
        nt_sd = scat_data[i_ss][i_se].ext_mat_data.nbooks();
        ARTS_USER_ERROR_IF(
            nt_sd != nt_se, bs1.str(), "ext_mat_data", bs2.str(), nt_se, ".")
        nt_sd = scat_data[i_ss][i_se].abs_vec_data.nbooks();
        ARTS_USER_ERROR_IF(
            nt_sd != nt_se, bs1.str(), "abs_vec_data", bs2.str(), nt_se, ".")
      }
    }
  }

  if (toupper(check_level) != "NONE") {
    // handing over to scat_dataCheck which checks whether
    // 1) scat_data containing any NaN?
    // 2) any negative values in Z11, K11, or a1?
    // 3) sca_mat norm sufficiently good (int(Z11)~=K11-a1?)
    // 1) & 2) always done
    // 3) only done if scat_data_check_level is "all"
    assert(false);
    // scat_dataCheck(scat_data, check_level, sca_mat_threshold);
  }

  // If here, all OK
  scat_data_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_checkedCalc(Index& sensor_checked,
                        const AscendingGrid& f_grid,
                        const Matrix& sensor_pos,
                        const Matrix& sensor_los,
                        const Matrix& transmitter_pos,
                        const Matrix& mblock_dlos,
                        const Sparse& sensor_response,
                        const Vector& sensor_response_f,
                        const ArrayOfIndex& sensor_response_pol,
                        const Matrix& sensor_response_dlos) {
  // Some sizes
  const Index nf = f_grid.nelem();
  const Index nlos = mblock_dlos.nrows();
  const Index n1y = sensor_response.nrows();
  const Index nmblock = sensor_pos.nrows();
  const Index niyb = nf * nlos * 4;

  // Sensor position and LOS.
  //
  ARTS_USER_ERROR_IF(sensor_pos.empty(),
                     "*sensor_pos* is empty. This is not allowed.");
  ARTS_USER_ERROR_IF(sensor_los.empty(),
                     "*sensor_los* is empty. This is not allowed.");
  //
  ARTS_USER_ERROR_IF(sensor_pos.ncols() != 3,
                     "The number of columns of sensor_pos must be "
                     "equal to the atmospheric dimensionality.");
  ARTS_USER_ERROR_IF(3 <= 2 && sensor_los.ncols() != 1,
                     "For 1D and 2D, sensor_los shall have one column.");
  ARTS_USER_ERROR_IF(3 == 3 && sensor_los.ncols() != 2,
                     "For 3D, sensor_los shall have two columns.");
  ARTS_USER_ERROR_IF(sensor_los.nrows() != nmblock,
                     "The number of rows of sensor_pos and sensor_los must be "
                     "identical, but sensor_pos has ",
                     nmblock,
                     " rows,\n"
                     "while sensor_los has ",
                     sensor_los.nrows(),
                     " rows.")
  ARTS_USER_ERROR_IF(
      max(sensor_los(joker, 0)) > 180,
      "First column of *sensor_los* is not allowed to have values above 180.");
  if (3 == 2) {
    ARTS_USER_ERROR_IF(
        min(sensor_los(joker, 0)) < -180,
        "For 3 = 2, first column of "
        "*sensor_los* is not allowed to have values below -180.");
  } else {
    ARTS_USER_ERROR_IF(min(sensor_los(joker, 0)) < 0,
                       "For 3 != 2, first column of "
                       "*sensor_los* is not allowed to have values below 0.");
  }
  if (3 == 3) {
    ARTS_USER_ERROR_IF(
        max(sensor_los(joker, 1)) > 180,
        "Second column of *sensor_los* is not allowed to have values above 180.");
    ARTS_USER_ERROR_IF(
        min(sensor_los(joker, 1)) < -180,
        "Second column of *sensor_los* is not allowed to have values below -180.");
  }

  // Transmission position.
  if (!transmitter_pos.empty()) {
    ARTS_USER_ERROR_IF(transmitter_pos.nrows() != sensor_pos.nrows(),
                       "*transmitter_pos* must either be empty or have "
                       "the same number of rows as *sensor_pos*.");
    ARTS_USER_ERROR_IF(transmitter_pos.ncols() != max(Index(2), 3),
                       "*transmitter_pos* must either be empty, have "
                       "2 for 1D/2D or 3 columns for 3D.");
  }

  // mblock_dlos
  //
  ARTS_USER_ERROR_IF(mblock_dlos.empty(), "*mblock_dlos* is empty.");
  ARTS_USER_ERROR_IF(mblock_dlos.ncols() > 2,
                     "The maximum number of columns in *mblock_dlos* is two.");
  if (3 < 3) {
    ARTS_USER_ERROR_IF(
        mblock_dlos.ncols() != 1,
        "For 1D and 2D *mblock_dlos* must have exactly one column.");
  }

  // Sensor
  //
  ARTS_USER_ERROR_IF(
      sensor_response.ncols() != niyb,
      "The *sensor_response* matrix does not have the right size,\n"
      "either the method *sensor_responseInit* has not been run or some\n"
      "of the other sensor response methods has not been correctly\n"
      "configured.")

  // Sensor aux variables
  //
  ARTS_USER_ERROR_IF(
      n1y != sensor_response_f.nelem() ||
          static_cast<Size>(n1y) != sensor_response_pol.size() ||
          n1y != sensor_response_dlos.nrows(),
      "Sensor auxiliary variables do not have the correct size.\n"
      "The following variables should all have same size:\n"
      "length of y for one block     : ",
      n1y,
      "\n"
      "sensor_response_f.nelem()     : ",
      sensor_response_f.nelem(),
      "\nsensor_response_pol.nelem() : ",
      sensor_response_pol.size(),
      "\nsensor_response_dlos.nrows(): ",
      sensor_response_dlos.nrows(),
      "\n")

  // If here, all OK
  sensor_checked = 1;
}
