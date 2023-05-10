/* Copyright (C) 2013
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler(at)ltu.se>

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

#include "arts.h"
#include "arts_conversions.h"
#include "atm.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "matpack_data.h"
#include "wigner_functions.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
using  Cloudbox::LAT_LON_MIN;

/* Workspace method: Doxygen documentation will be auto-generated */
void atmfields_checkedCalc(Index &atmfields_checked,
                           const ArrayOfArrayOfSpeciesTag &abs_species,
                           const AtmField &atm_field, const Verbosity &) {
  // Consistency between dim, grids and atmospheric fields/surfaces
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::p), "No pressure field")
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t), "No temperature field")
  for (auto &spec : abs_species)
    ARTS_USER_ERROR_IF(not atm_field.has(spec), "No ",
                       std::quoted(var_string(spec)), " field")

    {
      const bool u = atm_field.has(Atm::Key::mag_u),
                 v = atm_field.has(Atm::Key::mag_v),
                 w = atm_field.has(Atm::Key::mag_w);
      ARTS_USER_ERROR_IF(
          (u or v or w) and (not u or not v or not w),
          "If any magnetic field component exist, all three must.\n The "
          "following component exists (1) and does not exist (0):\n",
          "U: ", static_cast<int>(u), "V: ", static_cast<int>(v),
          "W: ", static_cast<int>(w))
    }

  {
    const bool u = atm_field.has(Atm::Key::wind_u),
               v = atm_field.has(Atm::Key::wind_v),
               w = atm_field.has(Atm::Key::wind_w);
    ARTS_USER_ERROR_IF(
        (u or v or w) and (not u or not v or not w),
        "If any wind field component exist, all three must.\n The "
        "following component exists (1) and does not exist (0):\n",
        "U: ", static_cast<int>(u), "V: ", static_cast<int>(v),
        "W: ", static_cast<int>(w))
  }

  // If here, all OK
  atmfields_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atmgeom_checkedCalc(Index& atmgeom_checked,
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& z_field,
                         const Vector& refellipsoid,  // Adopt to new version ZZZ
                         const Matrix& z_surface,
                         const Vector& lat_true,
                         const Vector& lon_true,
                         const Numeric& max500hpa_gradient,
                         const Verbosity&) {
  // A repetition from atmfields_checked, but we do this to make the two parts
  // independent (the other option would be to demand atmfields_checkec == 1)
  chk_atm_grids(3, p_grid, lat_grid, lon_grid);

  // *refellipsoid*
  ARTS_USER_ERROR_IF (refellipsoid.nelem() != 2,
        "The WSV *refellispoid* must be a vector of "
        "length 2*.");
  ARTS_USER_ERROR_IF (refellipsoid[0] <= 0,
        "The first element of *refellipsoid* must "
        "be > 0.");
  ARTS_USER_ERROR_IF (refellipsoid[1] < 0 || refellipsoid[1] > 1,
        "The second element of *refellipsoid* must be "
        "inside [0,1].");
  ARTS_USER_ERROR_IF (3 == 1 && refellipsoid[1] != 0,
        "For 1D, the second element of *refellipsoid* "
        "(the eccentricity) must be 0.");

  chk_atm_field("z_field", z_field, 3, p_grid, lat_grid, lon_grid);
  chk_atm_surface("z_surface", z_surface, 3, lat_grid, lon_grid);

  // Check that z_field has strictly increasing pages.
  for (Index row = 0; row < z_field.nrows(); row++) {
    for (Index col = 0; col < z_field.ncols(); col++) {
      ostringstream os;
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
      ARTS_USER_ERROR_IF (z_surface(row, col) < z_field(0, row, col) ||
                          z_surface(row, col) >= z_field(z_field.npages() - 1, row, col),
        "The surface altitude (*z_surface*) cannot be outside\n"
        "of the altitudes in *z_field*.\n"
        "z_surface: ", z_surface(row, col), "\n"
        "min of z_field: ", z_field(0, row, col), "\n"
        "max of z_field: ", z_field(z_field.npages() - 1, row, col),
        "\n",
        (3 > 1) ? var_string("\nThis was found to be the case for:\n", "latitude ", lat_grid[row]) : var_string(),
        (3 > 2) ? var_string("\nlongitude ", lon_grid[col]) : var_string())
    }
  }

  // A rough check of gradients in the 500 hPa level (or closest pressure level)
  // This check is just run in the latitude direction, along some longitudes to
  // save time
  if (3 > 1) {
    // Find 500 hPa
    Index ip = -1; Numeric dpmin = 99e99;
    for (Index i=0; i<p_grid.nelem(); i++) {
      const Numeric dp = abs(p_grid[i] - 500e2);
      if (dp < dpmin) {
        ip = i;
        dpmin = dp; 
      }
    }
    Numeric maxgrad = -1;
    for (Index ilon=0; ilon<max(1,lon_grid.nelem()); ilon += 10) {
      for (Index ilat=1; ilat<lat_grid.nelem(); ilat++) {
        const Numeric grad = abs((z_field(ip,ilat,ilon)-z_field(ip,ilat-1,ilon)) /
                                 (lat_grid[ilat]-lat_grid[ilat-1]));
        if (grad > maxgrad)
          maxgrad = grad;
      }
    }
    // Scale to change over 100 km
    maxgrad *= 100.0/111.0;
    if (maxgrad > max500hpa_gradient) {
      ostringstream os;
      os << "A check of the altitude of the " << p_grid[ip]/100
         << " hPa level has been made.\nThe maximum gradient found matches "
         << maxgrad << " m/100km, that exceeds\nthe set limit of "
         << max500hpa_gradient << " m/100km (by GIN *max500hpa_gradient*).\n"
         << "Please check the smoothness of *z_field*.";
      throw runtime_error(os.str());
    }
  }
  
  // lat/lon true
  if (3 < 3 && (lat_true.nelem() || lon_true.nelem())) {
    if (3 == 1) {
      ARTS_USER_ERROR_IF (lat_true.nelem() != 1,
                          "For 1D, *lat_true* must have length 1.");
      ARTS_USER_ERROR_IF (lon_true.nelem() != 1,
                          "For 1D, *lon_true* must have length 1.");
    } else if (3 == 2) {
      ARTS_USER_ERROR_IF (lat_true.nelem() != lat_grid.nelem(),
            "For 2D, *lat_true* must have same length as *lat_grid*.");
      ARTS_USER_ERROR_IF (lon_true.nelem() != lat_grid.nelem(),
            "For 2D, *lon_true* must have same length as *lat_grid*.");
    }
    ARTS_USER_ERROR_IF (lon_true.nelem() != lat_true.nelem(),
          "If *lat_true* is set, also *lon_true* must be "
          "set (and have the same length).");
    ARTS_USER_ERROR_IF (min(lat_true) < -90 || max(lat_true) > 90,
                        "Values in *lat_true* must be inside [-90,90].");
    ARTS_USER_ERROR_IF (min(lon_true) < -180 || max(lon_true) > 360,
                        "Values in *lon_true* must be inside [-180,360].");
  }

  // If here, all OK
  atmgeom_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_checkedCalc(Index& cloudbox_checked,
                          const Index& atmfields_checked,
                          const AtmField& atm_field,
                          const Matrix& z_surface,
                          const Index& cloudbox_on,
                          const ArrayOfIndex& cloudbox_limits,
                          const Tensor4& pnd_field,
                          const ArrayOfTensor4& dpnd_field_dx,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const ArrayOfString& scat_species,
                          const Matrix& particle_masses,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const Index& demand_latlon_margin,
                          const Index& negative_pnd_ok,
                          const Verbosity&) {

  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");

  chk_if_bool("cloudbox_on", cloudbox_on);

  if (cloudbox_on) {
  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field, wind_u_field;
  ARTS_USER_ERROR("ERROR")
  //const auto& [z_grid, lat_grid, lon_grid] = atm_field.grid;

    // Winds, must be empty variables (i.e. no winds allowed)
    ARTS_USER_ERROR_IF(
        atm_field.has(Atm::Key::wind_u) or atm_field.has(Atm::Key::wind_v) or
            atm_field.has(Atm::Key::wind_w),
        "The scattering methods are not (yet?) handling winds. For this\n"
        "reason, the WSVs for wind fields must all be empty with an\n."
        "active cloudbox.");

    // no "particles" in abs_species if cloudbox is on (they act on the same
    // scat_data! and there is no good reason to have some particles as
    // abs-only, if we anyway do a scattering calculation.).
    Index has_absparticles = 0;
    for (Index sp = 0; sp < abs_species.nelem() && has_absparticles < 1; sp++) {
      if (abs_species[sp].Particles()) {
        has_absparticles = 1;
      }
    }
    ARTS_USER_ERROR_IF (has_absparticles,
          "For scattering calculations (cloudbox is on),"
          "abs_species is not allowed to contain\n"
          "'particles' (absorbing-only particles)!");

    // Cloudbox limits
    ARTS_USER_ERROR_IF (cloudbox_limits.nelem() != 3 * 2,
        "The array *cloudbox_limits* has incorrect length.\n"
        "For atmospheric dim. = ", 3,
        " the length shall be ", 3 * 2, " but it is ",
        cloudbox_limits.nelem(), ".")
    ARTS_USER_ERROR_IF (cloudbox_limits[1] <= cloudbox_limits[0] || cloudbox_limits[0] < 0 ||
                        cloudbox_limits[1] >= z_grid.nelem(),
        "Incorrect value(s) for cloud box pressure limit(s) found."
        "\nValues are either out of range or upper limit is not "
        "greater than lower limit.\nWith present length of "
        "*p_grid*, OK values are 0 - ", z_grid.nelem() - 1,
        ".\nThe pressure index limits are set to ", cloudbox_limits[0],
        " - ", cloudbox_limits[1], ".")

    Index nlat = 1, nlon = 1;

    if (3 > 1) {
      nlat = lat_grid.nelem();
      if (demand_latlon_margin) {
        ARTS_USER_ERROR_IF (cloudbox_limits[3] <= cloudbox_limits[2] ||
                            cloudbox_limits[2] < 1 ||
                            cloudbox_limits[3] >= nlat - 1,
          "Incorrect value(s) for cloud box latitude limit(s) found."
          "\nValues are either out of range or upper limit is not "
          "greater than lower limit.\nWith present length of "
          "*lat_grid* and demand_latlon_margin set to true, "
          "OK values are 1 - ", nlat - 2,
          ".\nThe latitude index limits are set to ", cloudbox_limits[2],
          " - ", cloudbox_limits[3], ".")
        ARTS_USER_ERROR_IF ((lat_grid[cloudbox_limits[2]] -
                             lat_grid[0] < LAT_LON_MIN) &&
                            (3 == 2 ||
                             (3 == 3 && lat_grid[0] > -90)),
          "Too small distance between cloudbox and lower end of "
          "latitude grid.\n"
          "This distance must be ", LAT_LON_MIN, " degrees.\n"
          "Cloudbox ends at ", lat_grid[cloudbox_limits[2]],
          " and latitude grid starts at ", lat_grid[0], ".")
        ARTS_USER_ERROR_IF ((lat_grid[nlat - 1] -
                             lat_grid[cloudbox_limits[3]] < LAT_LON_MIN) &&
                            (3 == 2 ||
                             (3 == 3 && lat_grid[nlat - 1] < 90)),
          "Too small distance between cloudbox and upper end of "
          "latitude grid.\n"
          "This distance must be ", LAT_LON_MIN, " degrees.\n"
          "Cloudbox ends at ", lat_grid[cloudbox_limits[3]],
          " and latitude grid ends at ", lat_grid[nlat - 1], ".")
      } else {
        ARTS_USER_ERROR_IF (cloudbox_limits[3] <= cloudbox_limits[2] ||
                            cloudbox_limits[2] < 0 ||
                            cloudbox_limits[3] >= nlat,
          "Incorrect value(s) for cloud box latitude limit(s) found."
          "\nValues are either out of range or upper limit is not "
          "greater than lower limit.\nWith present length of "
          "*lat_grid* and demand_latlon_margin set to false, "
          "OK values are 0 - ", nlat - 1,
          ".\nThe latitude index limits are set to ", cloudbox_limits[2],
          " - ", cloudbox_limits[3], ".")
      }  
    }

    if (3 > 2) {
      nlon = lon_grid.nelem();
      if (demand_latlon_margin) {
        ARTS_USER_ERROR_IF (cloudbox_limits[5] <= cloudbox_limits[4] ||
                            cloudbox_limits[4] < 1 ||
                            cloudbox_limits[5] >= nlon - 1,
          "Incorrect value(s) for cloud box longitude limit(s) found"
          ".\nValues are either out of range or upper limit is not "
          "greater than lower limit.\nWith present length of "
          "*lon_grid* and demand_latlon_margin set to true,"
          "OK values are 1 - ", nlon - 2,
          ".\nThe longitude limits are set to ", cloudbox_limits[4],
          " - ", cloudbox_limits[5], ".")
        if (lon_grid[nlon - 1] - lon_grid[0] < 360) {
          const Numeric latmax = max(abs(lat_grid[cloudbox_limits[2]]),
                                     abs(lat_grid[cloudbox_limits[3]]));
          const Numeric lfac = 1 / cos(DEG2RAD * latmax);
          ARTS_USER_ERROR_IF (lon_grid[cloudbox_limits[4]] - lon_grid[0] <
                              LAT_LON_MIN / lfac,
            "Too small distance between cloudbox and lower end of"
            "the longitude\ngrid. This distance must here be ",
            LAT_LON_MIN / lfac, " degrees.")
          ARTS_USER_ERROR_IF (lon_grid[nlon - 1] - lon_grid[cloudbox_limits[5]] <
                              LAT_LON_MIN / lfac,
            "Too small distance between cloudbox and upper end of"
            "the longitude\ngrid. This distance must here be ",
            LAT_LON_MIN / lfac, " degrees.")
        }
      } else {
        ARTS_USER_ERROR_IF (cloudbox_limits[5] <= cloudbox_limits[4] ||
                            cloudbox_limits[4] < 0 ||
                            cloudbox_limits[5] >= nlon,
          "Incorrect value(s) for cloud box longitude limit(s) found"
          ".\nValues are either out of range or upper limit is not "
          "greater than lower limit.\nWith present length of "
          "*lon_grid* and demand_latlon_margin set to false,"
          "OK values are 0 - ", nlon - 1,
          ".\nThe longitude limits are set to ", cloudbox_limits[4],
          " - ", cloudbox_limits[5], ".")
      }
    }

    // Check with respect to z_surface
    for (Index o = 0; o < nlon; o++) {
      for (Index a = 0; a < nlat; a++) {
        ARTS_USER_ERROR_IF (z_grid[cloudbox_limits[1]] <= z_surface(a, o),
              "The upper vertical limit of the cloudbox must be above "
              "the surface altitude (for all latitudes and longitudes).");
      }
    }

    // Check pnd_field
    //
    const Index np = TotalNumberOfElements(scat_data);
    // Dummy variables to mimic grids of correct size
    Vector g1(cloudbox_limits[1] - cloudbox_limits[0] + 1), g2(0), g3(0);
    if (3 > 1) {
      g2.resize(cloudbox_limits[3] - cloudbox_limits[2] + 1);
    }
    if (3 > 2) {
      g3.resize(cloudbox_limits[5] - cloudbox_limits[4] + 1);
    }

    chk_atm_field("pnd_field", pnd_field, 3, np, g1, g2, g3);

    ARTS_USER_ERROR_IF (!negative_pnd_ok && min(pnd_field) < 0,
                        "Negative values in *pnd_field* not allowed.");

    // No non-zero pnd at lower boundary unless lower boundary is at or below
    // surface
    for (Index a = 0; a < g2.nelem(); a++) {
      for (Index o = 0; o < g3.nelem(); o++) {
        ARTS_USER_ERROR_IF (max(pnd_field(joker, 0, a, o)) > 0 &&
            z_grid[cloudbox_limits[0]] > z_surface(a, o),
              "A non-zero value found in *pnd_field* at the"
              " lower altitude limit of the cloudbox (but the "
              "position is not at or below the surface altitude).");
      }
    }

    // No non-zero pnd at upper boundary unless upper boundary is top of
    // atmosphere
    if (cloudbox_limits[1] != z_grid.nelem() - 1)
      ARTS_USER_ERROR_IF (max(pnd_field(joker, g1.nelem() - 1, joker, joker)) > 0,
            "A non-zero value found in *pnd_field* at "
            "upper altitude limit of the cloudbox.");
    if (3 >= 2) {
      ARTS_USER_ERROR_IF (max(pnd_field(joker, joker, 0, joker)) > 0,
            "A non-zero value found in *pnd_field* at "
            "lower latitude limit of the cloudbox.");
      ARTS_USER_ERROR_IF (max(pnd_field(joker, joker, g2.nelem() - 1, joker)) > 0,
            "A non-zero value found in *pnd_field* at "
            "upper latitude limit of the cloudbox.");
    }
    if (3 == 3) {
      ARTS_USER_ERROR_IF (max(pnd_field(joker, joker, joker, 0)) > 0,
            "A non-zero value found in *pnd_field* at "
            "lower longitude limit of the cloudbox.");
      ARTS_USER_ERROR_IF (max(pnd_field(joker, joker, joker, g3.nelem() - 1)) > 0,
            "A non-zero value found in *pnd_field* at "
            "upper longitude limit of the cloudbox.");
    }

    // And dpnd_field_dx
    ARTS_USER_ERROR_IF (dpnd_field_dx.nelem() != jacobian_quantities.nelem(),
          "Size of *dpnd_field_dx* inconsistent with number "
          "of *jacobian_quantities*.");

    // Check semi-mandatory variables, that are allowed to be empty
    //

    // scat_species:
    if (scat_species.nelem() > 0)
      ARTS_USER_ERROR_IF (scat_species.nelem() != scat_data.nelem(),
          "Number of scattering species specified by scat_species does\n"
          "not agree with number of scattering species in scat_data:\n"
          "scat_species has ", scat_species.nelem(),
          " entries, while scat_data has ", scat_data.nelem(), ".")

    // particle_masses:
    if (!particle_masses.empty()) {
      ARTS_USER_ERROR_IF (particle_masses.nrows() != np,
            "The WSV *particle_masses* must either be "
            "empty or have a row size matching the "
            "length of *scat_data*.");
      ARTS_USER_ERROR_IF (min(particle_masses) < 0,
                          "All values in *particles_masses* must be >= 0.");
    }
  }

  // If here, all OK
  cloudbox_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_checkedCalc(Index& scat_data_checked,
                           const ArrayOfArrayOfSingleScatteringData& scat_data,
                           const Vector& f_grid,
                           const Numeric& dfrel_threshold,
                           const String& check_level,
                           const Numeric& sca_mat_threshold,
                           const Verbosity& verbosity)
// FIXME: when we allow K, a, Z to be on different f and T grids, their use in
// the scatt solvers needs to be reviewed again and adapted to this!
{
  // Prevent the user from producing nonsense by using too much freedom in
  // setting the single freq validity threshold...
  ARTS_USER_ERROR_IF (dfrel_threshold > 0.5,
      "*dfrel_threshold* too large (max. allowed: 0.5, your's: ",
      dfrel_threshold, ").")

  // freq range of calc covered?
  ARTS_USER_ERROR_IF (f_grid.empty(), "The frequency grid is empty.");
  if (f_grid.nelem() > 1) chk_if_increasing("f_grid", f_grid);

  Index nf = f_grid.nelem();
  Index N_ss = scat_data.nelem();
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    Index N_se = scat_data[i_ss].nelem();
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // For each scattering element (se) check that se's f_grid is either
      // identical to f_grid or contains a single entry only. In the latter
      // case, the f/f_grid-ratio (equv. to a size parameter ratio) not allowed
      // to exceed the dfrel_threshold.
      Index nf_se = scat_data[i_ss][i_se].f_grid.nelem();
      if (nf_se != 1) {
        ARTS_USER_ERROR_IF (nf_se != nf,
            "*scat_data* must have either one or *f_grid* (=", nf,
            ") frequency entries,\n"
            "but scattering element #", i_se, " in scattering species #",
            i_ss, " has ", nf_se, ".")
        for (Index f = 0; f < nf_se; f++) {
          ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                  scat_data[i_ss][i_se].f_grid[f], f_grid[f], 0.5e-9),
              "*scat_data* frequency grid has to be identical to *f_grid*\n"
              "(or contain only a single entry),\n"
              "but scattering element #", i_se,
              " in scattering species #", i_ss,
              " deviates for f_index ", f, ".")
        }
      } else {
        ARTS_USER_ERROR_IF ((abs(1. - scat_data[i_ss][i_se].f_grid[0] / f_grid[0]) >
             dfrel_threshold) or
            (abs(1. - scat_data[i_ss][i_se].f_grid[0] / f_grid[nf - 1]) >
             dfrel_threshold),
            "Frequency entry (f=", scat_data[i_ss][i_se].f_grid[0],
            "Hz) of scattering element #", i_se, "\n"
            "in scattering species #", i_ss, " is too far (>",
            dfrel_threshold * 1e2, "%) from one or more\n"
            "of the f_grid limits (fmin=", f_grid[0],
            "Hz, fmax=", f_grid[nf - 1], "Hz).")
      }

      // check that the freq dimension of sca_mat, ext_mat, and abs_vec is
      // either ssd.f_grid.nelem() or 1 (FIXME: so far, being freq dim !=
      // ssd.f_grid.nelem() switched off, as usage in scatt solvers so far
      // doesn't allow this. see FIXME at start.).
      {
        ostringstream bs1, bs2;
        bs1 << "Frequency dimension of ";
        //bs2 << " must be either one or ssd.f_grid.nelem() (=" << nf_se << "),\n"
        bs2 << " must be ssd.f_grid.nelem() (=" << nf_se << "),\n"
            << "but scattering element #" << i_se << " in scattering species #"
            << i_ss << " is ";
        Index nf_sd = scat_data[i_ss][i_se].pha_mat_data.nlibraries();
        ARTS_USER_ERROR_IF (nf_sd != nf_se,
          bs1.str(), "pha_mat_data", bs2.str(), nf_se, ".")
        nf_sd = scat_data[i_ss][i_se].ext_mat_data.nshelves();
        ARTS_USER_ERROR_IF (nf_sd != nf_se,
          bs1.str(), "ext_mat_data", bs2.str(), nf_se, ".")
        nf_sd = scat_data[i_ss][i_se].abs_vec_data.nshelves();
        ARTS_USER_ERROR_IF (nf_sd != nf_se,
          bs1.str(), "abs_vec_data", bs2.str(), nf_se, ".")
      }

      // check that the temp dimension of K and a is ssd.T_grid.nelem(). For Z
      // it might be ssd.T_grid.nelem() or 1.
      {
        ostringstream bs1, bs2;
        Index nt_se = scat_data[i_ss][i_se].T_grid.nelem();
        bs1 << "Temperature dimension of ";
        //bs2 << " must be either one or ssd.T_grid.nelem(),\n"
        bs2 << " must be ssd.T_grid.nelem() (=" << nt_se << "),\n"
            << "but for scattering element #" << i_se
            << " in scattering species #" << i_ss << " it is ";
        Index nt_sd = scat_data[i_ss][i_se].pha_mat_data.nvitrines();
        ARTS_USER_ERROR_IF (nt_sd != nt_se and nt_sd != 1,
          bs1.str(), "pha_mat_data", bs2.str(), nt_sd, ".")
        nt_sd = scat_data[i_ss][i_se].ext_mat_data.nbooks();
        ARTS_USER_ERROR_IF (nt_sd != nt_se,
          bs1.str(), "ext_mat_data", bs2.str(), nt_se, ".")
        nt_sd = scat_data[i_ss][i_se].abs_vec_data.nbooks();
        ARTS_USER_ERROR_IF (nt_sd != nt_se,
          bs1.str(), "abs_vec_data", bs2.str(), nt_se, ".")
      }
    }
  }

  if (check_level.toupper() != "NONE") {
    // handing over to scat_dataCheck which checks whether
    // 1) scat_data containing any NaN?
    // 2) any negative values in Z11, K11, or a1?
    // 3) sca_mat norm sufficiently good (int(Z11)~=K11-a1?)
    // 1) & 2) always done
    // 3) only done if scat_data_check_level is "all"
    scat_dataCheck(scat_data, check_level, sca_mat_threshold, verbosity);
  }

  // If here, all OK
  scat_data_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void lbl_checkedCalc(Index& lbl_checked,
                     const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const SpeciesIsotopologueRatios& isotopologue_ratios,
                     const Verbosity&)
{
  checkIsotopologueRatios(abs_lines_per_species, isotopologue_ratios);
  checkPartitionFunctions(abs_lines_per_species);
  
  lbl_checked = false;
  
  ARTS_USER_ERROR_IF (abs_lines_per_species.nelem() not_eq abs_species.nelem(),
      "abs_lines_per_species and abs_species must have same length.\n"
      "Instead len(abs_lines_per_species) = ",
      abs_lines_per_species.nelem(),
      " and len(abs_species) = ",
      abs_species.nelem(),
      '\n')
  
  for (Index i=0; i<abs_species.nelem(); i++) {
    auto& specs = abs_species[i];
    auto& lines = abs_lines_per_species[i];
    
    if (not specs.nelem()) {
      if (not lines.nelem()) {
        continue;
      }
      ARTS_USER_ERROR ( "Lines for non-existent species discovered!\n");
    }
    
    const bool any_zeeman = std::any_of(specs.cbegin(), specs.cend(), [](auto& x){return x.Type() == Species::TagType::Zeeman;});
    ARTS_USER_ERROR_IF (any_zeeman and (not std::all_of(specs.cbegin(), specs.cend(), [](auto& x){return x.Type() == Species::TagType::Zeeman;})),
      "Zeeman species found but not all sub-species tags support Zeeman effect.\n"
      "Offending tag: ", specs, '\n')
    
    if (any_zeeman) {
      for (auto& band: lines) {
        ARTS_USER_ERROR_IF (band.cutoff not_eq Absorption::CutoffType::None,
                            "Zeeman effects are not symmetric, you cannot use cutoff.\n");
        for (Index k=0; k<band.NumLines(); k++) {
          bool hasJ = band.lines[k].localquanta.val.has(QuantumNumberType::J);
          bool hasF = band.lines[k].localquanta.val.has(QuantumNumberType::F);
          ARTS_USER_ERROR_IF (not(hasF or hasJ),
                              "No J(s) or F(s) yet declared Zeeman splitting.\n");

          auto& qn = hasF ? band.lines[k].localquanta.val[QuantumNumberType::F] : band.lines[k].localquanta.val[QuantumNumberType::J];
          ARTS_USER_ERROR_IF (not is_wigner3_ready(qn.upp()),
                              "Bad Wigner numbers for upper state F or J.  Try increasing the Wigner memory allocation.\n");
          ARTS_USER_ERROR_IF (not is_wigner3_ready(qn.low()),
                              "Bad Wigner numbers for lower state F or J.  Try increasing the Wigner memory allocation.\n");

          auto Ze = band.lines[k].zeeman;
          ARTS_USER_ERROR_IF (Ze.gu() == 0 ? false : not std::isnormal(Ze.gu()),
                              "Bad value(s) in the upper Zeeman data not allowed when modeling Zeeman effect.\n");
          ARTS_USER_ERROR_IF (Ze.gl() == 0 ? false : not std::isnormal(Ze.gl()),
                              "Bad value(s) in the lower Zeeman data not allowed when modeling Zeeman effect.\n");
        }
      }
    } else /*if (not any any_zeeman)*/ {
    }
    
    // Checks per band
    for (auto& band: lines) {
      ARTS_USER_ERROR_IF (band.mirroring not_eq Absorption::MirroringType::Manual and
                          std::any_of(band.lines.cbegin(), band.lines.cend(),
                                      [](auto& x){return x.F0 <= 0;}),
                          "Negative or zero frequency in non-Manual mirrored band.\n");
    }
  }
  
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      
      // Mirroring checks
      for (auto& line: band.lines) {
        if (band.mirroring == Absorption::MirroringType::Manual) {
          ARTS_USER_ERROR_IF(line.F0 >= 0,
                             "Must have negative frequency, finds " , line.F0)
        } else {
          ARTS_USER_ERROR_IF(line.F0 <= 0,
                             "Must have positive frequency, finds " , line.F0)
        }
      }
      
      // Cutoff checks
      switch (band.cutoff) {
        case Absorption::CutoffType::None: break;
        case Absorption::CutoffType::ByLine: {
          ARTS_USER_ERROR_IF (not (band.mirroring == Absorption::MirroringType::None or
                                   band.mirroring == Absorption::MirroringType::Manual),
                              "Cutoff only possible with symmetric mirroring types")
          ARTS_USER_ERROR_IF(Absorption::relaxationtype_relmat(band.population),
                             "Cannot have relaxation matrix line mixing with cutoff calculations")
          ARTS_USER_ERROR_IF (not (band.lineshapetype == LineShape::Type::DP or
                                   band.lineshapetype == LineShape::Type::LP or
                                   band.lineshapetype == LineShape::Type::VP),
                              "Cutoff only possible with symmetric line shape types")
          for (auto& line: band.lines) {
            for (auto& single_data: line.lineshape.Data()) {
              // Skipping G() intentionally, since the calculations technically works with it
              ARTS_USER_ERROR_IF(
                single_data.Y().type not_eq LineShape::TemperatureModel::None or
                single_data.DV().type not_eq LineShape::TemperatureModel::None,
                                 "Cannot have Rosenkranz-style line mixing with cutoff calculations\n"
                                 "Note that abs_lines_per_speciesTurnOffLineMixing will make this error go away by modifying the data\n"
              )
            }
          }
        } break;
        case Absorption::CutoffType::FINAL: ARTS_USER_ERROR("You have a band with undefined cutoff type.")
      }
    }
  }
  
  lbl_checked = true;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearsky_agenda_checkedCalc(
    Workspace& ws _U_,
    Index& propmat_clearsky_agenda_checked,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Agenda& propmat_clearsky_agenda,
    const Verbosity&) {
  bool needs_lines = false;
  bool needs_zeeman = false;
  bool needs_predefined = false;
  bool needs_continua = false;
  bool needs_cia = false;
  bool needs_free_electrons = false;
  bool needs_particles = false;
  bool needs_hxsec = false;

  for (auto& tag_groups: abs_species) {
    for (auto& tag: tag_groups) {
      switch (tag.Type()) {
        case Species::TagType::Plain:
          needs_lines = true;
          break;
        case Species::TagType::Zeeman:
          needs_zeeman = true;
          break;
        case Species::TagType::Predefined:
          needs_predefined = true;
          break;
        case Species::TagType::Cia:
          needs_cia = true;
          break;
        case Species::TagType::FreeElectrons:
          needs_free_electrons = true;
          break;
        case Species::TagType::Particles:
          needs_particles = true;
          break;
        case Species::TagType::XsecFit:
          needs_hxsec = true;
          break;
        default:
          ARTS_ASSERT(false, "Unknown species type: ", tag.Type())
          break;
      }
    }
  }

  ARTS_USER_ERROR_IF(
      needs_lines and
          not(propmat_clearsky_agenda.has_method("propmat_clearskyAddLines") or
              propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddFromLookup")),
      "*abs_species* contains normal lines species but *propmat_clearsky_agenda*\n"
      "cannot support lines.");

  ARTS_USER_ERROR_IF(
      needs_continua and
          not(propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddXsecAgenda") or
              propmat_clearsky_agenda.has_method("propmat_clearskyAddConts") or
              propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddFromLookup")),
      "*abs_species* contains legacy continua but *propmat_clearsky_agenda*\n"
      "cannot support these.");

  ARTS_USER_ERROR_IF(
      needs_cia and
          not(propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddXsecAgenda") or
              propmat_clearsky_agenda.has_method("propmat_clearskyAddCIA") or
              propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddFromLookup")),
      "*abs_species* contains CIA models but *propmat_clearsky_agenda*\n"
      "cannot support these.");

  ARTS_USER_ERROR_IF(
      needs_hxsec and
          not(propmat_clearsky_agenda.has_method("propmat_clearskyAddXsecAgenda") or
              propmat_clearsky_agenda.has_method("propmat_clearskyAddXsecFit") or
              propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddFromLookup")),
      "*abs_species* contains Hitran XSEC models but *propmat_clearsky_agenda*\n"
      "cannot support these.");

  ARTS_USER_ERROR_IF(
      needs_predefined and
          not(propmat_clearsky_agenda.has_method("propmat_clearskyAddPredefined") or
              propmat_clearsky_agenda.has_method(
                  "propmat_clearskyAddFromLookup")),
      "*abs_species* contains modern continua models but *propmat_clearsky_agenda*\n"
      "cannot support these.");

  ARTS_USER_ERROR_IF (needs_zeeman and
      not propmat_clearsky_agenda.has_method("propmat_clearskyAddZeeman"),
        "*abs_species* contains Zeeman species but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddZeeman*.");

  ARTS_USER_ERROR_IF (needs_free_electrons and
      not propmat_clearsky_agenda.has_method("propmat_clearskyAddFaraday"),
        "*abs_species* contains Zeeman species but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddFaraday*.");
  
  ARTS_USER_ERROR_IF (needs_particles &&
      !(propmat_clearsky_agenda.has_method("propmat_clearskyAddParticles")),
        "*abs_species* contains particles but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddParticles*.");

  propmat_clearsky_agenda_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_checkedCalc(Index& sensor_checked,
                        const Index& stokes_dim,
                        const Vector& f_grid,
                        const Matrix& sensor_pos,
                        const Matrix& sensor_los,
                        const Matrix& transmitter_pos,
                        const Matrix& mblock_dlos,
                        const Sparse& sensor_response,
                        const Vector& sensor_response_f,
                        const ArrayOfIndex& sensor_response_pol,
                        const Matrix& sensor_response_dlos,
                        const Verbosity&) {
  // Some sizes
  const Index nf = f_grid.nelem();
  const Index nlos = mblock_dlos.nrows();
  const Index n1y = sensor_response.nrows();
  const Index nmblock = sensor_pos.nrows();
  const Index niyb = nf * nlos * stokes_dim;

  // Sensor position and LOS.
  //
  ARTS_USER_ERROR_IF (!is_increasing(f_grid),
                      "*f_grid* must be a strictly increasing vector.");

  // Sensor position and LOS.
  //
  ARTS_USER_ERROR_IF (sensor_pos.empty(),
                      "*sensor_pos* is empty. This is not allowed.");
  ARTS_USER_ERROR_IF (sensor_los.empty(),
                      "*sensor_los* is empty. This is not allowed.");
  //
  ARTS_USER_ERROR_IF (sensor_pos.ncols() != 3,
        "The number of columns of sensor_pos must be "
        "equal to the atmospheric dimensionality.");
  ARTS_USER_ERROR_IF (3 <= 2 && sensor_los.ncols() != 1,
                      "For 1D and 2D, sensor_los shall have one column.");
  ARTS_USER_ERROR_IF (3 == 3 && sensor_los.ncols() != 2,
                      "For 3D, sensor_los shall have two columns.");
  ARTS_USER_ERROR_IF (sensor_los.nrows() != nmblock,
      "The number of rows of sensor_pos and sensor_los must be "
      "identical, but sensor_pos has ", nmblock, " rows,\n"
      "while sensor_los has ", sensor_los.nrows(), " rows.")
  ARTS_USER_ERROR_IF (max(sensor_los(joker, 0)) > 180,
        "First column of *sensor_los* is not allowed to have values above 180.");
  if (3 == 2) {
    ARTS_USER_ERROR_IF (min(sensor_los(joker, 0)) < -180,
          "For 3 = 2, first column of "
          "*sensor_los* is not allowed to have values below -180.");
  } else {
    ARTS_USER_ERROR_IF (min(sensor_los(joker, 0)) < 0,
          "For 3 != 2, first column of "
          "*sensor_los* is not allowed to have values below 0.");
  }
  if (3 == 3) {
    ARTS_USER_ERROR_IF (max(sensor_los(joker, 1)) > 180,
          "Second column of *sensor_los* is not allowed to have values above 180.");
    ARTS_USER_ERROR_IF (min(sensor_los(joker, 1)) < -180,
          "Second column of *sensor_los* is not allowed to have values below -180.");
  }

  // Transmission position.
  if (!transmitter_pos.empty()) {
    ARTS_USER_ERROR_IF (transmitter_pos.nrows() != sensor_pos.nrows(),
          "*transmitter_pos* must either be empty or have "
          "the same number of rows as *sensor_pos*.");
    ARTS_USER_ERROR_IF (transmitter_pos.ncols() != max(Index(2), 3),
          "*transmitter_pos* must either be empty, have "
          "2 for 1D/2D or 3 columns for 3D.");
  }

  // mblock_dlos
  //
  ARTS_USER_ERROR_IF (mblock_dlos.empty(),
        "*mblock_dlos* is empty.");
  ARTS_USER_ERROR_IF (mblock_dlos.ncols() > 2,
        "The maximum number of columns in *mblock_dlos* is two.");
  if (3 < 3) {
    ARTS_USER_ERROR_IF (mblock_dlos.ncols() != 1,
          "For 1D and 2D *mblock_dlos* must have exactly one column.");
  }

  // Sensor
  //
  ARTS_USER_ERROR_IF (sensor_response.ncols() != niyb,
      "The *sensor_response* matrix does not have the right size,\n"
      "either the method *sensor_responseInit* has not been run or some\n"
      "of the other sensor response methods has not been correctly\n"
      "configured.")

  // Sensor aux variables
  //
  ARTS_USER_ERROR_IF (n1y != sensor_response_f.nelem() || n1y != sensor_response_pol.nelem() ||
      n1y != sensor_response_dlos.nrows(),
      "Sensor auxiliary variables do not have the correct size.\n"
      "The following variables should all have same size:\n"
      "length of y for one block     : ", n1y, "\n"
      "sensor_response_f.nelem()     : ", sensor_response_f.nelem(),
      "\nsensor_response_pol.nelem() : ", sensor_response_pol.nelem(),
      "\nsensor_response_dlos.nrows(): ", sensor_response_dlos.nrows(),
      "\n")

  // If here, all OK
  sensor_checked = 1;
}
