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
#include "auto_md.h"
#include "cloudbox.h"
#include "matpackI.h"

extern const Numeric DEG2RAD;
extern const Numeric LAT_LON_MIN;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_agenda_checkedCalc(Workspace& ws _U_,
                                 Index& abs_xsec_agenda_checked,
                                 // WS Input:
                                 const ArrayOfArrayOfSpeciesTag& abs_species,
                                 const Agenda& abs_xsec_agenda,
                                 const Verbosity&) {
  bool needs_lines = false;
  bool needs_continua = false;
  bool needs_cia = false;
  bool needs_hxsec = false;

  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    for (Index tgs = 0; tgs < abs_species[sp].nelem(); tgs++) {
      switch (abs_species[sp][tgs].Type()) {
        case SpeciesTag::TYPE_PLAIN:
          needs_lines = true;
          break;
        case SpeciesTag::TYPE_ZEEMAN:
          break;
        case SpeciesTag::TYPE_PREDEF:
          needs_continua = true;
          break;
        case SpeciesTag::TYPE_CIA:
          needs_cia = true;
          break;
        case SpeciesTag::TYPE_FREE_ELECTRONS:
          break;
        case SpeciesTag::TYPE_PARTICLES:
          break;
        case SpeciesTag::TYPE_HITRAN_XSEC:
          needs_hxsec = true;
          break;
        default:
          ostringstream os;
          os << "Unknown species type: " << abs_species[sp][tgs].Type();
          throw runtime_error(os.str());
          break;
      }
    }
  }

  if (needs_lines and
      not(abs_xsec_agenda.has_method("abs_xsec_per_speciesAddLines") or
          abs_xsec_agenda.has_method("abs_xsec_per_speciesAddLines2"))) {
    throw runtime_error(
        "*abs_species* contains line species but *abs_xsec_agenda*\n"
        "does not contain *abs_xsec_per_speciesAddLines*.");
  }

  if (needs_continua &&
      !abs_xsec_agenda.has_method("abs_xsec_per_speciesAddConts")) {
    throw runtime_error(
        "*abs_species* contains continuum species but *abs_xsec_agenda*\n"
        "does not contain *abs_xsec_per_speciesAddConts*.");
  }

  if (needs_cia && !abs_xsec_agenda.has_method("abs_xsec_per_speciesAddCIA")) {
    throw runtime_error(
        "*abs_species* contains CIA species but *abs_xsec_agenda*\n"
        "does not contain *abs_xsec_per_speciesAddCIA*.");
  }

  if (needs_hxsec &&
      !abs_xsec_agenda.has_method("abs_xsec_per_speciesAddHitranXsec")) {
    throw runtime_error(
        "*abs_species* contains HITRAN xsec species but *abs_xsec_agenda*\n"
        "does not contain *abs_xsec_per_speciesAddHitranXsec*.");
  }

  // If here, all OK
  abs_xsec_agenda_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atmfields_checkedCalc(Index& atmfields_checked,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Tensor3& t_field,
                           const Tensor4& vmr_field,
                           const Tensor3& wind_u_field,
                           const Tensor3& wind_v_field,
                           const Tensor3& wind_w_field,
                           const Tensor3& mag_u_field,
                           const Tensor3& mag_v_field,
                           const Tensor3& mag_w_field,
                           const SpeciesAuxData& partition_functions,
                           const Index& abs_f_interp_order,
                           const Index& negative_vmr_ok,
                           const Index& bad_partition_functions_ok,
                           const Verbosity&) {
  // Consistency between dim, grids and atmospheric fields/surfaces
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);
  chk_atm_field("t_field", t_field, atmosphere_dim, p_grid, lat_grid, lon_grid);
  chk_atm_field("vmr_field",
                vmr_field,
                atmosphere_dim,
                abs_species.nelem(),
                p_grid,
                lat_grid,
                lon_grid);

  // More for vmr_field.
  if (!negative_vmr_ok && abs_species.nelem() && min(vmr_field) < 0)
    throw runtime_error("All values in *vmr_field* must be >= 0.");

  // More for t_field.
  if (min(t_field) <= 0)
    throw runtime_error("All temperatures in *t_field* must be > 0.");

  // Winds
  if (wind_w_field.npages() > 0) {
    chk_atm_field("wind_w_field",
                  wind_w_field,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid);
  }
  if (atmosphere_dim < 3 && wind_v_field.npages() > 0) {
    chk_atm_field("wind_v_field",
                  wind_v_field,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid);
  }
  if (atmosphere_dim > 2) {
    if (wind_u_field.npages() > 0) {
      if (wind_v_field.npages() > 0) {
        bool chk_poles = false;
        chk_atm_field("wind_u_field",
                      wind_u_field,
                      atmosphere_dim,
                      p_grid,
                      lat_grid,
                      lon_grid,
                      chk_poles);
        chk_atm_field("wind_v_field",
                      wind_v_field,
                      atmosphere_dim,
                      p_grid,
                      lat_grid,
                      lon_grid,
                      chk_poles);
        chk_atm_vecfield_lat90("wind_v_field",
                               wind_v_field,
                               "wind_u_field",
                               wind_u_field,
                               atmosphere_dim,
                               lat_grid);
      } else {
        chk_atm_field("wind_u_field",
                      wind_u_field,
                      atmosphere_dim,
                      p_grid,
                      lat_grid,
                      lon_grid);
      }
    } else {
      if (wind_v_field.npages() > 0) {
        chk_atm_field("wind_v_field",
                      wind_v_field,
                      atmosphere_dim,
                      p_grid,
                      lat_grid,
                      lon_grid);
      }
    }
  }

  // If any of the wind fields exist, abs_f_interp_order must not be zero.
  if (wind_u_field.npages() > 0 || wind_v_field.npages() > 0 ||
      wind_w_field.npages() > 0) {
    if (abs_f_interp_order == 0) {
      ostringstream os;
      os << "You have a wind field set, but abs_f_interp_order zero.\n"
         << "This is not allowed. Though abs_f_interp_order only is\n"
         << "required and has an effect if absorption lookup tables\n"
         << "are used, for safety reasons you also have to set it >0\n"
         << "in case of on-the-fly absorption.";
      throw runtime_error(os.str());
    }
  }

  // Magnetic field
  if (mag_w_field.npages() > 0) {
    chk_atm_field("mag_w_field (vertical magfield component)",
                  mag_w_field,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid);
  }
  if (mag_u_field.npages() > 0) {
    if (mag_v_field.npages() > 0) {
      bool chk_poles = false;
      chk_atm_field("mag_v_field",
                    mag_v_field,
                    atmosphere_dim,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    chk_poles);
      chk_atm_field("mag_u_field",
                    mag_u_field,
                    atmosphere_dim,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    chk_poles);
      chk_atm_vecfield_lat90("mag_v_field",
                             mag_v_field,
                             "mag_u_field",
                             mag_u_field,
                             atmosphere_dim,
                             lat_grid);
    } else {
      chk_atm_field("mag_u_field",
                    mag_u_field,
                    atmosphere_dim,
                    p_grid,
                    lat_grid,
                    lon_grid);
    }
  } else {
    if (mag_v_field.npages() > 0) {
      chk_atm_field("mag_v_field",
                    mag_v_field,
                    atmosphere_dim,
                    p_grid,
                    lat_grid,
                    lon_grid);
    }
  }

  if (partition_functions.nspecies()) {
    // Partition functions have partition functions and finding their temperature limits?
    Numeric min_T = 0, max_T = 1e6;  //1 MK is max tested
    for (Index ii = 0; ii < partition_functions.nspecies(); ii++)
      for (Index jj = 0; jj < partition_functions.nisotopologues(ii); jj++) {
        // Test if species is important at all (maybe, since code is unclear; will at worst cause failures later on in the program)
        bool test_spec = false;
        for (auto& as : abs_species)
          for (auto& s : as)
            if (s.Species() == ii and s.Isotopologue() == jj) test_spec = true;
        if (not test_spec) continue;

        ArrayOfGriddedField1 part_fun;
        switch (partition_functions.getParamType(ii, jj)) {
          case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
            part_fun = partition_functions.getParam(ii, jj);
            if (part_fun.nelem() == 2) {
              if (part_fun[1].data.nelem() == 2 &&
                  part_fun[0].data.nelem() > 1) {
                if (part_fun[1].data[0] > min_T) min_T = part_fun[1].data[0];
                if (part_fun[1].data[1] < max_T) max_T = part_fun[1].data[1];
              } else
                throw std::runtime_error(
                    "Bad coefficient parameter in partition_function.\n");
            } else
              throw std::runtime_error(
                  "Bad coefficient parameter in partition_function.\n");
            break;

          case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF_VIBROT:
            part_fun = partition_functions.getParam(ii, jj);
            if (part_fun.nelem() == 3) {
              if (part_fun[2].data.nelem() == 2 &&
                  (part_fun[1].data.nelem() == part_fun[0].data.nelem())) {
                if (part_fun[2].data[0] > min_T) min_T = part_fun[2].data[0];
                if (part_fun[2].data[1] < max_T) max_T = part_fun[2].data[1];
              } else
                throw std::runtime_error(
                    "Bad coefficient parameter in partition_function.\n");
            } else
              throw std::runtime_error(
                  "Bad coefficient parameter in partition_function.\n");
            break;

          case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
            part_fun = partition_functions.getParam(ii, jj);
            if (part_fun.nelem() == 1) {
              if (part_fun[0].data.nelem() > 1) {
                if (part_fun[0].get_numeric_grid(0)[0] > min_T)
                  min_T = part_fun[0].get_numeric_grid(0)[0];
                if (part_fun[0].get_numeric_grid(
                        0)[part_fun[0].data.nelem() - 1] < max_T)
                  max_T = part_fun[0].get_numeric_grid(
                      0)[part_fun[0].data.nelem() - 1];
              } else
                throw std::runtime_error(
                    "Bad t_field parameter in partition_function.\n");
            } else
              throw std::runtime_error(
                  "Bad t_field parameter in partition_function.\n");
            break;

          default:
            throw std::runtime_error(
                "Bad parameter type in partition_functions.\n");
            break;
        }
      }

    // Check that partition functions are OK if not explicitly turned off
    if (!bad_partition_functions_ok) {
      for (Index ii = 0; ii < t_field.npages(); ii++)
        for (Index jj = 0; jj < t_field.nrows(); jj++)
          for (Index kk = 0; kk < t_field.ncols(); kk++) {
            if (min_T > t_field(ii, jj, kk) || max_T < t_field(ii, jj, kk)) {
              ostringstream os;
              os << "There are bad partition functions in your setup.\n"
                 << "Minimum temperature for defined partition functions is: "
                 << min_T
                 << " K.\nMaximum temperature for defined partition functions is: "
                 << max_T << " K\nThere is a t_field entry of "
                 << t_field(ii, jj, kk) << " K.\n";
              throw std::runtime_error(os.str());
            }
          }
    }
  }

  // If here, all OK
  atmfields_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atmgeom_checkedCalc(Index& atmgeom_checked,
                         const Index& atmosphere_dim,
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& z_field,
                         const Vector& refellipsoid,
                         const Matrix& z_surface,
                         const Vector& lat_true,
                         const Vector& lon_true,
                         const Verbosity&) {
  // A repetition from atmfields_checked, but we do this to make the two parts
  // independent (the other option would be to demand atmfields_checkec == 1)
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  // *refellipsoid*
  if (refellipsoid.nelem() != 2)
    throw runtime_error(
        "The WSV *refellispoid* must be a vector of "
        "length 2*.");
  if (refellipsoid[0] <= 0)
    throw runtime_error(
        "The first element of *refellipsoid* must "
        "be > 0.");
  if (refellipsoid[1] < 0 || refellipsoid[1] > 1)
    throw runtime_error(
        "The second element of *refellipsoid* must be "
        "inside [0,1].");
  if (atmosphere_dim == 1 && refellipsoid[1] != 0)
    throw runtime_error(
        "For 1D, the second element of *refellipsoid* "
        "(the eccentricity) must be 0.");

  chk_atm_field("z_field", z_field, atmosphere_dim, p_grid, lat_grid, lon_grid);
  chk_atm_surface("z_surface", z_surface, atmosphere_dim, lat_grid, lon_grid);

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
      if (z_surface(row, col) < z_field(0, row, col) ||
          z_surface(row, col) >= z_field(z_field.npages() - 1, row, col)) {
        ostringstream os;
        os << "The surface altitude (*z_surface*) cannot be outside\n"
           << "of the altitudes in *z_field*.\n"
           << "z_surface: " << z_surface(row, col) << "\n"
           << "min of z_field: " << z_field(0, row, col) << "\n"
           << "max of z_field: " << z_field(z_field.npages() - 1, row, col)
           << "\n";
        if (atmosphere_dim > 1)
          os << "\nThis was found to be the case for:\n"
             << "latitude " << lat_grid[row];
        if (atmosphere_dim > 2) os << "\nlongitude " << lon_grid[col];
        throw runtime_error(os.str());
      }
    }
  }

  // lat/lon true
  if (atmosphere_dim < 3 && (lat_true.nelem() || lon_true.nelem())) {
    if (atmosphere_dim == 1) {
      if (lat_true.nelem() != 1)
        throw runtime_error("For 1D, *lat_true* must have length 1.");
      if (lon_true.nelem() != 1)
        throw runtime_error("For 1D, *lon_true* must have length 1.");
    } else if (atmosphere_dim == 2) {
      if (lat_true.nelem() != lat_grid.nelem())
        throw runtime_error(
            "For 2D, *lat_true* must have same length as *lat_grid*.");
      if (lon_true.nelem() != lat_grid.nelem())
        throw runtime_error(
            "For 2D, *lon_true* must have same length as *lat_grid*.");
    }
    if (lon_true.nelem() != lat_true.nelem())
      throw runtime_error(
          "If *lat_true* is set, also *lon_true* must be "
          "set (and have the same length).");
    if (min(lat_true) < -90 || max(lat_true) > 90)
      throw runtime_error("Values in *lat_true* must be inside [-90,90].");
    if (min(lon_true) < -180 || max(lon_true) > 360)
      throw runtime_error("Values in *lon_true* must be inside [-180,360].");
  }

  // If here, all OK
  atmgeom_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_checkedCalc(Index& cloudbox_checked,
                          const Index& atmfields_checked,
                          const Index& atmosphere_dim,
                          const Vector& p_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Tensor3& z_field,
                          const Matrix& z_surface,
                          const Tensor3& wind_u_field,
                          const Tensor3& wind_v_field,
                          const Tensor3& wind_w_field,
                          const Index& cloudbox_on,
                          const ArrayOfIndex& cloudbox_limits,
                          const Tensor4& pnd_field,
                          const ArrayOfTensor4& dpnd_field_dx,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const ArrayOfString& scat_species,
                          const Matrix& particle_masses,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const Index& negative_pnd_ok,
                          const Verbosity&) {
  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");

  chk_if_bool("cloudbox_on", cloudbox_on);

  if (cloudbox_on) {
    // Winds, must be empty variables (i.e. no winds allowed)
    {
      ostringstream ow;
      ow << "The scattering methods are not (yet?) handling winds. For this\n"
         << "reason, the WSVs for wind fields must all be empty with an\n."
         << "active cloudbox.";
      if (!wind_w_field.empty() || !wind_v_field.empty() ||
          !wind_u_field.empty()) {
        throw runtime_error(ow.str());
      }
    }

    // no "particles" in abs_species if cloudbox is on (they act on the same
    // scat_data! and there is no good reason to have some particles as
    // abs-only, if we anyway do a scattering calculation.).
    Index has_absparticles = 0;
    for (Index sp = 0; sp < abs_species.nelem() && has_absparticles < 1; sp++) {
      if (abs_species[sp][0].Type() == SpeciesTag::TYPE_PARTICLES) {
        has_absparticles = 1;
      }
    }
    if (has_absparticles) {
      throw runtime_error(
          "For scattering calculations (cloudbox is on),"
          "abs_species is not allowed to contain\n"
          "'particles' (absorbing-only particles)!");
    }

    // Cloudbox limits
    if (cloudbox_limits.nelem() != atmosphere_dim * 2) {
      ostringstream os;
      os << "The array *cloudbox_limits* has incorrect length.\n"
         << "For atmospheric dim. = " << atmosphere_dim
         << " the length shall be " << atmosphere_dim * 2 << " but it is "
         << cloudbox_limits.nelem() << ".";
      throw runtime_error(os.str());
    }
    if (cloudbox_limits[1] <= cloudbox_limits[0] || cloudbox_limits[0] < 0 ||
        cloudbox_limits[1] >= p_grid.nelem()) {
      ostringstream os;
      os << "Incorrect value(s) for cloud box pressure limit(s) found."
         << "\nValues are either out of range or upper limit is not "
         << "greater than lower limit.\nWith present length of "
         << "*p_grid*, OK values are 0 - " << p_grid.nelem() - 1
         << ".\nThe pressure index limits are set to " << cloudbox_limits[0]
         << " - " << cloudbox_limits[1] << ".";
      throw runtime_error(os.str());
    }

    Index nlat = 1, nlon = 1;

    if (atmosphere_dim > 1) {
      nlat = lat_grid.nelem();
      if (cloudbox_limits[3] <= cloudbox_limits[2] || cloudbox_limits[2] < 1 ||
          cloudbox_limits[3] >= nlat - 1) {
        ostringstream os;
        os << "Incorrect value(s) for cloud box latitude limit(s) found."
           << "\nValues are either out of range or upper limit is not "
           << "greater than lower limit.\nWith present length of "
           << "*lat_grid*, OK values are 1 - " << nlat - 2
           << ".\nThe latitude index limits are set to " << cloudbox_limits[2]
           << " - " << cloudbox_limits[3] << ".";
        throw runtime_error(os.str());
      }
      if ((lat_grid[cloudbox_limits[2]] - lat_grid[0] < LAT_LON_MIN) &&
          (atmosphere_dim == 2 || (atmosphere_dim == 3 && lat_grid[0] > -90))) {
        ostringstream os;
        os << "Too small distance between cloudbox and lower end of "
           << "latitude grid.\n"
           << "This distance must be " << LAT_LON_MIN << " degrees.\n"
           << "Cloudbox ends at " << lat_grid[cloudbox_limits[2]]
           << " and latitude grid starts at " << lat_grid[0] << ".";
        throw runtime_error(os.str());
      }
      if ((lat_grid[nlat - 1] - lat_grid[cloudbox_limits[3]] < LAT_LON_MIN) &&
          (atmosphere_dim == 2 ||
           (atmosphere_dim == 3 && lat_grid[nlat - 1] < 90))) {
        ostringstream os;
        os << "Too small distance between cloudbox and upper end of "
           << "latitude grid.\n"
           << "This distance must be " << LAT_LON_MIN << " degrees.\n"
           << "Cloudbox ends at " << lat_grid[cloudbox_limits[3]]
           << " and latitude grid ends at " << lat_grid[nlat - 1] << ".";
        throw runtime_error(os.str());
      }
    }

    if (atmosphere_dim > 2) {
      nlon = lon_grid.nelem();
      if (cloudbox_limits[5] <= cloudbox_limits[4] || cloudbox_limits[4] < 1 ||
          cloudbox_limits[5] >= nlon - 1) {
        ostringstream os;
        os << "Incorrect value(s) for cloud box longitude limit(s) found"
           << ".\nValues are either out of range or upper limit is not "
           << "greater than lower limit.\nWith present length of "
           << "*lon_grid*, OK values are 1 - " << nlon - 2
           << ".\nThe longitude limits are set to " << cloudbox_limits[4]
           << " - " << cloudbox_limits[5] << ".";
        throw runtime_error(os.str());
      }
      if (lon_grid[nlon - 1] - lon_grid[0] < 360) {
        const Numeric latmax = max(abs(lat_grid[cloudbox_limits[2]]),
                                   abs(lat_grid[cloudbox_limits[3]]));
        const Numeric lfac = 1 / cos(DEG2RAD * latmax);
        if (lon_grid[cloudbox_limits[4]] - lon_grid[0] < LAT_LON_MIN / lfac) {
          ostringstream os;
          os << "Too small distance between cloudbox and lower end of"
             << "the longitude\ngrid. This distance must here be "
             << LAT_LON_MIN / lfac << " degrees.";
          throw runtime_error(os.str());
        }
        if (lon_grid[nlon - 1] - lon_grid[cloudbox_limits[5]] <
            LAT_LON_MIN / lfac) {
          ostringstream os;
          os << "Too small distance between cloudbox and upper end of"
             << "the longitude\ngrid. This distance must here be "
             << LAT_LON_MIN / lfac << " degrees.";
          throw runtime_error(os.str());
        }
      }
    }

    // Check with respect to z_surface
    for (Index o = 0; o < nlon; o++) {
      for (Index a = 0; a < nlat; a++) {
        if (z_field(cloudbox_limits[1], a, o) <= z_surface(a, o))
          throw runtime_error(
              "The upper vertical limit of the cloudbox must be above "
              "the surface altitude (for all latitudes and longitudes).");
      }
    }

    // Check pnd_field
    //
    const Index np = TotalNumberOfElements(scat_data);
    // Dummy variables to mimic grids of correct size
    Vector g1(cloudbox_limits[1] - cloudbox_limits[0] + 1), g2(0), g3(0);
    if (atmosphere_dim > 1) {
      g2.resize(cloudbox_limits[3] - cloudbox_limits[2] + 1);
    }
    if (atmosphere_dim > 2) {
      g3.resize(cloudbox_limits[5] - cloudbox_limits[4] + 1);
    }

    chk_atm_field("pnd_field", pnd_field, atmosphere_dim, np, g1, g2, g3);

    if (!negative_pnd_ok && min(pnd_field) < 0)
      throw runtime_error("Negative values in *pnd_field* not allowed.");

    // No non-zero pnd at lower boundary unless lower boundary is at or below
    // surface
    for (Index a = 0; a < g2.nelem(); a++) {
      for (Index o = 0; o < g3.nelem(); o++) {
        if (max(pnd_field(joker, 0, a, o)) > 0 &&
            z_field(cloudbox_limits[0], a, o) > z_surface(a, o))
          throw runtime_error(
              "A non-zero value found in *pnd_field* at the"
              " lower altitude limit of the cloudbox (but the "
              "position is not at or below the surface altitude).");
      }
    }

    // No non-zero pnd at upper boundary unless upper boundary is top of
    // atmosphere
    if (cloudbox_limits[1] != p_grid.nelem() - 1)
      if (max(pnd_field(joker, g1.nelem() - 1, joker, joker)) > 0)
        throw runtime_error(
            "A non-zero value found in *pnd_field* at "
            "upper altitude limit of the cloudbox.");
    if (atmosphere_dim >= 2) {
      if (max(pnd_field(joker, joker, 0, joker)) > 0)
        throw runtime_error(
            "A non-zero value found in *pnd_field* at "
            "lower latitude limit of the cloudbox.");
      if (max(pnd_field(joker, joker, g2.nelem() - 1, joker)) > 0)
        throw runtime_error(
            "A non-zero value found in *pnd_field* at "
            "upper latitude limit of the cloudbox.");
    }
    if (atmosphere_dim == 3) {
      if (max(pnd_field(joker, joker, joker, 0)) > 0)
        throw runtime_error(
            "A non-zero value found in *pnd_field* at "
            "lower longitude limit of the cloudbox.");
      if (max(pnd_field(joker, joker, joker, g3.nelem() - 1)) > 0)
        throw runtime_error(
            "A non-zero value found in *pnd_field* at "
            "upper longitude limit of the cloudbox.");
    }

    // And dpnd_field_dx
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "Size of *dpnd_field_dx* inconsistent with number "
          "of *jacobian_quantities*.");

    // Check semi-mandatory variables, that are allowed to be empty
    //

    // scat_species:
    if (scat_species.nelem() > 0)
      if (scat_species.nelem() != scat_data.nelem()) {
        ostringstream os;
        os << "Number of scattering species specified by scat_species does\n"
           << "not agree with number of scattering species in scat_data:\n"
           << "scat_species has " << scat_species.nelem()
           << " entries, while scat_data has " << scat_data.nelem() << ".";
        throw runtime_error(os.str());
      }

    // particle_masses:
    if (!particle_masses.empty()) {
      if (particle_masses.nrows() != np)
        throw runtime_error(
            "The WSV *particle_masses* must either be "
            "empty or have a row size matching the "
            "length of *scat_data*.");
      if (min(particle_masses) < 0)
        throw runtime_error("All values in *particles_masses* must be >= 0.");
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
  if (dfrel_threshold > 0.5) {
    ostringstream os;
    os << "*dfrel_threshold* too large (max. allowed: 0.5, your's: "
       << dfrel_threshold << ").";
    throw runtime_error(os.str());
  }

  // freq range of calc covered?
  if (f_grid.empty()) throw runtime_error("The frequency grid is empty.");
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
        if (nf_se != nf) {
          ostringstream os;
          os << "*scat_data* must have either one or *f_grid* (=" << nf
             << ") frequency entries,\n"
             << "but scattering element #" << i_se << " in scattering species #"
             << i_ss << " has " << nf_se << ".";
          throw runtime_error(os.str());
        } else {
          for (Index f = 0; f < nf_se; f++) {
            if (!is_same_within_epsilon(
                    scat_data[i_ss][i_se].f_grid[f], f_grid[f], 0.5e-9)) {
              ostringstream os;
              os << "*scat_data* frequency grid has to be identical to *f_grid*\n"
                 << "(or contain only a single entry),\n"
                 << "but scattering element #" << i_se
                 << " in scattering species #" << i_ss
                 << " deviates for f_index " << f << ".";
              throw runtime_error(os.str());
            }
          }
        }
      } else {
        if ((abs(1. - scat_data[i_ss][i_se].f_grid[0] / f_grid[0]) >
             dfrel_threshold) or
            (abs(1. - scat_data[i_ss][i_se].f_grid[0] / f_grid[nf - 1]) >
             dfrel_threshold)) {
          ostringstream os;
          os << "Frequency entry (f=" << scat_data[i_ss][i_se].f_grid[0]
             << "Hz) of scattering element #" << i_se << "\n"
             << "in scattering species #" << i_ss << " is too far (>"
             << dfrel_threshold * 1e2 << "%) from one or more\n"
             << "of the f_grid limits (fmin=" << f_grid[0]
             << "Hz, fmax=" << f_grid[nf - 1] << "Hz).";
          throw runtime_error(os.str());
        }
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
        if (nf_sd != nf_se)  //&& nf_sd != 1)
        {
          ostringstream os;
          os << bs1.str() << "pha_mat_data" << bs2.str() << nf_se << ".";
          throw runtime_error(os.str());
        }
        nf_sd = scat_data[i_ss][i_se].ext_mat_data.nshelves();
        if (nf_sd != nf_se)  //&& nf_sd != 1)
        {
          ostringstream os;
          os << bs1.str() << "ext_mat_data" << bs2.str() << nf_se << ".";
          throw runtime_error(os.str());
        }
        nf_sd = scat_data[i_ss][i_se].abs_vec_data.nshelves();
        if (nf_sd != nf_se)  //&& nf_sd != 1)
        {
          ostringstream os;
          os << bs1.str() << "abs_vec_data" << bs2.str() << nf_se << ".";
          throw runtime_error(os.str());
        }
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
        if (nt_sd != nt_se and nt_sd != 1) {
          ostringstream os;
          os << bs1.str() << "pha_mat_data" << bs2.str() << nt_sd << ".";
          throw runtime_error(os.str());
        }
        nt_sd = scat_data[i_ss][i_se].ext_mat_data.nbooks();
        if (nt_sd != nt_se)  // no need to check for 1 here. since if it is 1,
                             // also T_grid.nelem need to be 1
        {
          ostringstream os;
          os << bs1.str() << "ext_mat_data" << bs2.str() << nt_se << ".";
          throw runtime_error(os.str());
        }
        nt_sd = scat_data[i_ss][i_se].abs_vec_data.nbooks();
        if (nt_sd != nt_se) {
          ostringstream os;
          os << bs1.str() << "abs_vec_data" << bs2.str() << nt_se << ".";
          throw runtime_error(os.str());
        }
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
                     const SpeciesAuxData& isotopologue_ratios,
                     const SpeciesAuxData& partition_functions,
                     const Verbosity&)
{
  checkIsotopologueRatios(abs_species, isotopologue_ratios);
  checkPartitionFunctions(abs_species, partition_functions);
  
  lbl_checked = false;
  
  if (abs_lines_per_species.nelem() not_eq abs_species.nelem()) {
    std::ostringstream os;
    os << "abs_lines_per_species and abs_species must have same length.\n"
    << "Instead len(abs_lines_per_species) = " << abs_lines_per_species.nelem() << ' '
    << "and len(abs_species) = " << abs_species.nelem() << '\n';
    throw std::runtime_error(os.str());
  }
  
  for (Index i=0; i<abs_species.nelem(); i++) {
    auto& specs = abs_species[i];
    auto& lines = abs_lines_per_species[i];
    
    if (not specs.nelem()) {
      if (not lines.nelem()) {
        continue;
      } else {
        throw std::runtime_error("Lines for non-existent species discovered!\n");
      }
    }
    
    const bool any_zeeman = specs[0].Type() == SpeciesTag::TYPE_ZEEMAN;
    for (auto& s: specs) {
      if (s.Type() not_eq SpeciesTag::TYPE_ZEEMAN) {
        if (any_zeeman) {
          std::ostringstream os;
          os << "Zeeman species found but not all sub-species support Zeeman\n";
          os << "Offending tag: " << specs << '\n';
          throw std::runtime_error(os.str());
        }
      }
    }
    
    if (any_zeeman) {
      for (auto& band: lines) {
        for (Index k=0; k<band.NumLines(); k++) {
          auto Ju = band.UpperQuantumNumber(k, QuantumNumberType::J);
          auto Jl = band.LowerQuantumNumber(k, QuantumNumberType::J);
          if (Ju.isDefined() and not is_Wigner3_ready(Ju)) {
            throw std::runtime_error("Bad wigner numbers for upper state J.  Try increasing the Wigner memory allocation.\n");
          } else if (Jl.isDefined() and not is_Wigner3_ready(Jl)) {
            throw std::runtime_error("Bad wigner numbers for lower state J.  Try increasing the Wigner memory allocation.\n");
          }
        }
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
  bool needs_continua = false;
  bool needs_cia = false;
  //bool needs_free_electrons = false;
  bool needs_particles = false;
  bool needs_hxsec = false;

  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    for (Index tgs = 0; tgs < abs_species[sp].nelem(); tgs++) {
      switch (abs_species[sp][tgs].Type()) {
        case SpeciesTag::TYPE_PLAIN:
          needs_lines = true;
          break;
        case SpeciesTag::TYPE_ZEEMAN:
          needs_zeeman = true;
          break;
        case SpeciesTag::TYPE_PREDEF:
          needs_continua = true;
          break;
        case SpeciesTag::TYPE_CIA:
          needs_cia = true;
          break;
        case SpeciesTag::TYPE_FREE_ELECTRONS:
          break;
        case SpeciesTag::TYPE_PARTICLES:
          needs_particles = true;
          break;
        case SpeciesTag::TYPE_HITRAN_XSEC:
          needs_hxsec = true;
          break;
        default:
          ostringstream os;
          os << "Unknown species type: " << abs_species[sp][tgs].Type();
          throw runtime_error(os.str());
          break;
      }
    }
  }

  if ((needs_lines || needs_continua || needs_cia || needs_hxsec) &&
      !(propmat_clearsky_agenda.has_method("propmat_clearskyAddOnTheFly") ||
        propmat_clearsky_agenda.has_method("propmat_clearskyAddFromLookup"))) {
    throw runtime_error(
        "*abs_species* contains line species, CIA species, "
        "hitran xsec species or continua but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddOnTheFly* nor "
        "*propmat_clearskyAddFromLookup*.");
  }

  if (needs_zeeman and
      not propmat_clearsky_agenda.has_method("propmat_clearskyAddZeeman")) {
    throw runtime_error(
        "*abs_species* contains Zeeman species but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddZeeman*.");
  }
  /*
    if (needs_free_electrons
        && !propmat_clearsky_agenda.has_method("propmat_clearskyAddFaraday"))
    {
        throw runtime_error("*abs_species* contains free electrons but *propmat_clearsky_agenda*\n"
                            "does not contain *propmat_clearskyAddFaraday*.");
    }
*/
  if (needs_particles &&
      !(propmat_clearsky_agenda.has_method("propmat_clearskyAddParticles") ||
        propmat_clearsky_agenda.has_method("propmat_clearskyAddParticles2"))) {
    throw runtime_error(
        "*abs_species* contains particles but *propmat_clearsky_agenda*\n"
        "does not contain *propmat_clearskyAddParticles*.");
  }

  propmat_clearsky_agenda_checked = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_checkedCalc(Index& sensor_checked,
                        const Index& atmosphere_dim,
                        const Index& stokes_dim,
                        const Vector& f_grid,
                        const Matrix& sensor_pos,
                        const Matrix& sensor_los,
                        const Matrix& transmitter_pos,
                        const Matrix& mblock_dlos_grid,
                        const Sparse& sensor_response,
                        const Vector& sensor_response_f,
                        const ArrayOfIndex& sensor_response_pol,
                        const Matrix& sensor_response_dlos,
                        const Verbosity&) {
  // Some sizes
  const Index nf = f_grid.nelem();
  const Index nlos = mblock_dlos_grid.nrows();
  const Index n1y = sensor_response.nrows();
  const Index nmblock = sensor_pos.nrows();
  const Index niyb = nf * nlos * stokes_dim;

  // Sensor position and LOS.
  //
  if (!is_increasing(f_grid))
    throw runtime_error("*f_grid* must be a strictly increasing vector.");

  // Sensor position and LOS.
  //
  if (sensor_pos.empty())
    throw runtime_error("*sensor_pos* is empty. This is not allowed.");
  if (sensor_los.empty())
    throw runtime_error("*sensor_los* is empty. This is not allowed.");
  //
  if (sensor_pos.ncols() != atmosphere_dim)
    throw runtime_error(
        "The number of columns of sensor_pos must be "
        "equal to the atmospheric dimensionality.");
  if (atmosphere_dim <= 2 && sensor_los.ncols() != 1)
    throw runtime_error("For 1D and 2D, sensor_los shall have one column.");
  if (atmosphere_dim == 3 && sensor_los.ncols() != 2)
    throw runtime_error("For 3D, sensor_los shall have two columns.");
  if (sensor_los.nrows() != nmblock) {
    ostringstream os;
    os << "The number of rows of sensor_pos and sensor_los must be "
       << "identical, but sensor_pos has " << nmblock << " rows,\n"
       << "while sensor_los has " << sensor_los.nrows() << " rows.";
    throw runtime_error(os.str());
  }
  if (max(sensor_los(joker, 0)) > 180)
    throw runtime_error(
        "First column of *sensor_los* is not allowed to have values above 180.");
  if (atmosphere_dim == 2) {
    if (min(sensor_los(joker, 0)) < -180)
      throw runtime_error(
          "For atmosphere_dim = 2, first column of "
          "*sensor_los* is not allowed to have values below -180.");
  } else {
    if (min(sensor_los(joker, 0)) < 0)
      throw runtime_error(
          "For atmosphere_dim != 2, first column of "
          "*sensor_los* is not allowed to have values below 0.");
  }
  if (atmosphere_dim == 3) {
    if (max(sensor_los(joker, 1)) > 180)
      throw runtime_error(
          "Second column of *sensor_los* is not allowed to have values above 180.");
    else if (min(sensor_los(joker, 1)) < -180)
      throw runtime_error(
          "Second column of *sensor_los* is not allowed to have values below -180.");
  }

  // Transmission position.
  if (transmitter_pos.ncols() > 0 && transmitter_pos.nrows() > 0) {
    if (transmitter_pos.nrows() != sensor_pos.nrows())
      throw runtime_error(
          "*transmitter_pos* must either be empty or have "
          "the same number of rows as *sensor_pos*.");
    if (transmitter_pos.ncols() != max(Index(2), atmosphere_dim))
      throw runtime_error(
          "*transmitter_pos* must either be empty, have "
          "2 for 1D/2D or 3 columns for 3D.");
  }

  // mblock_dlos_grid
  //
  if (mblock_dlos_grid.empty())
    throw runtime_error("*mblock_dlos_grid* is empty.");
  if (mblock_dlos_grid.ncols() > 2)
    throw runtime_error(
        "The maximum number of columns in *mblock_dlos_grid* is two.");
  if (atmosphere_dim < 3) {
    if (mblock_dlos_grid.ncols() != 1)
      throw runtime_error(
          "For 1D and 2D *mblock_dlos_grid* must have exactly one column.");
  }

  // Sensor
  //
  if (sensor_response.ncols() != niyb) {
    ostringstream os;
    os << "The *sensor_response* matrix does not have the right size,\n"
       << "either the method *sensor_responseInit* has not been run or some\n"
       << "of the other sensor response methods has not been correctly\n"
       << "configured.";
    throw runtime_error(os.str());
  }

  // Sensor aux variables
  //
  if (n1y != sensor_response_f.nelem() || n1y != sensor_response_pol.nelem() ||
      n1y != sensor_response_dlos.nrows()) {
    ostringstream os;
    os << "Sensor auxiliary variables do not have the correct size.\n"
       << "The following variables should all have same size:\n"
       << "length of y for one block     : " << n1y << "\n"
       << "sensor_response_f.nelem()     : " << sensor_response_f.nelem()
       << "\nsensor_response_pol.nelem() : " << sensor_response_pol.nelem()
       << "\nsensor_response_dlos.nrows(): " << sensor_response_dlos.nrows()
       << "\n";
    throw runtime_error(os.str());
  }

  // If here, all OK
  sensor_checked = 1;
}
