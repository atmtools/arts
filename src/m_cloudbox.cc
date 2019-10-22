/* Copyright (C) 2002-2017
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
   Claudia Emde     <claudia.emde@dlr.de>
   Cory Davis       <cory.davis@metservice.com>	   
   Jana Mendrok     <jana.mendrok@gmail.com>
   Daniel Kreyling  <daniel.kreyling@nict.go.jp>
   Manfred Brath    <manfred.brath@uni-hamburg.de>
                         
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
  \file   m_cloudbox.cc
  \author Patrick Eriksson, Claudia Emde and Sreerekha T. R.
  \date   2002-05-08 

  \brief  Workspace functions related to the definintion of the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "file.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "microphysics.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric LAT_LON_MIN;
extern const Numeric DENSITY_OF_ICE;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxOff(Index& cloudbox_on,
                 Index& ppath_inside_cloudbox_do,
                 ArrayOfIndex& cloudbox_limits,
                 Agenda& iy_cloudbox_agenda,
                 Tensor4& pnd_field,
                 ArrayOfTensor4& dpnd_field_dx,
                 ArrayOfString& scat_species,
                 ArrayOfArrayOfSingleScatteringData& scat_data,
                 ArrayOfArrayOfSingleScatteringData& scat_data_raw,
                 Index& scat_data_checked,
                 Matrix& particle_masses,
                 const ArrayOfRetrievalQuantity& jacobian_quantities,
                 const Verbosity&) {
  cloudbox_on = 0;
  ppath_inside_cloudbox_do = 0;
  cloudbox_limits.resize(0);
  iy_cloudbox_agenda = Agenda();
  iy_cloudbox_agenda.set_name("iy_cloudbox_agenda");
  pnd_field.resize(0, 0, 0, 0);
  // we need to size dpnd_field to be consistent with jacobian_quantities.
  dpnd_field_dx.resize(jacobian_quantities.nelem());
  scat_data.resize(0);
  scat_species.resize(0);
  // remove scat_data_raw resizing once all scat solvers have been convert to
  // use of (new-type) scat_data
  scat_data_raw.resize(0);
  scat_data_checked = 0;
  particle_masses.resize(0, 0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetAutomatically(  // WS Output:
    //Workspace& /* ws */,
    Index& cloudbox_on,
    ArrayOfIndex& cloudbox_limits,
    //Agenda&  iy_cloudbox_agenda,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor4& particle_field,
    // Control Parameters
    const ArrayOfIndex& cloudbox_limits_old,
    const Numeric& cloudbox_margin,
    const Verbosity& verbosity) {
  // Check existing WSV
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  // includes p_grid chk_if_decresasing
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);
  // Set cloudbox_on
  cloudbox_on = 1;

  if (atmosphere_dim > 1) {
    ostringstream os;
    os << "cloudboxSetAutomatically not yet available for 2D and 3D cases.";
    throw runtime_error(os.str());
  }

  Index np = p_grid.nelem();

  // Allocate cloudbox_limits
  cloudbox_limits.resize(atmosphere_dim * 2);

  bool cb_preset = (min(cloudbox_limits_old) > -1);
  if (cb_preset) {
    if (cloudbox_limits_old.nelem() != atmosphere_dim * 2) {
      ostringstream os;
      os << "The array *cloudbox_limits_old* has incorrect length.\n"
         << "For atmospheric dim. = " << atmosphere_dim
         << " the length shall be " << atmosphere_dim * 2 << " but it is "
         << cloudbox_limits_old.nelem() << ".";
      throw runtime_error(os.str());
    }
  }

  // Initialize boundary counters
  Index p1 = 0, p2 = 0;
  if (cloudbox_margin != -1) {
    if (cb_preset)
      p1 = cloudbox_limits_old[0] + 1;
    else
      p1 = np - 1;
  }
  if (cb_preset) {
    p2 = cloudbox_limits_old[1] - 1;
  }

  // OLE: Commented out until code that uses it at the end of this function is commented back in
  /*
  if ( atmosphere_dim > 1 )
    {
      Index lat1 = particle_field.nrows()-1;
      Index lat2 = 0;
    }
  if ( atmosphere_dim > 2 )
    {
      Index lon1 = particle_field.ncols()-1;
      Index lon2 = 0;
    }
*/

  bool one_not_empty = false;
  bool any_not_empty = false;

  if (!particle_field.empty()) {
    Index nss = particle_field.nbooks();

    //--------- Start loop over fields in particle_field------------------------
    for (Index l = 0; l < nss; l++) {
      //not_empty is set to true, if a single value of particle_field
      //is unequal zero (and not NaN), i.e. if we actually have some amount of
      //these scattering species in the atmosphere.
      chk_scat_species_field(one_not_empty,
                             particle_field(l, joker, joker, joker),
                             "particle_field",
                             atmosphere_dim,
                             p_grid,
                             lat_grid,
                             lon_grid);

      //if particles found, enter detailed search
      if (one_not_empty) {
        any_not_empty = true;
        find_cloudlimits(p1,
                         p2,
                         particle_field(l, joker, joker, joker),
                         atmosphere_dim,
                         cloudbox_margin);
      }
    }
  }

  if (any_not_empty || cb_preset) {
    // decrease lower cb limit by one to ensure that linear interpolation of
    // particle number densities is possible.
    p1 = max(p1 - 1, Index(0));

    Numeric p_margin1;

    // alter lower cloudbox_limit by cloudbox_margin, using barometric
    // height formula
    p_margin1 = barometric_heightformula(p_grid[p1], cloudbox_margin);
    while ((p_grid[p1] < p_margin1) && (p1 > 0)) p1--;
    cloudbox_limits[0] = p1;

    // increase upper cb limit by one to ensure that linear interpolation of
    // particle number densities is possible.
    p2 = min(p2 + 1, np - 1);
    // set upper cloudbox_limit
    // if cloudbox reaches to the upper most pressure level
    if (p2 >= np - 1) {
      CREATE_OUT2;
      out2 << "The cloud reaches to TOA!\n"
           << "Check your *particle_field* data, if realistic!\n";
    }
    cloudbox_limits[1] = p2;

    // assert keyword arguments

    // The pressure in *p1* must be bigger than the pressure in *p2*.
    assert(p_grid[p1] > p_grid[p2]);
    // The pressure in *p1* must be larger than the last value in *p_grid*.
    assert(p_grid[p1] > p_grid[p_grid.nelem() - 1]);
    // The pressure in *p2* must be smaller than the first value in *p_grid*."
    assert(p_grid[p2] < p_grid[0]);

    /*
    if ( atmosphere_dim >= 2 )
    {
      // The latitude in *lat2* must be bigger than the latitude in *lat1*.
      assert ( lat_grid[lat2] > lat_grid[lat1] );
      // The latitude in *lat1* must be >= the second value in *lat_grid*.
      assert ( lat_grid[lat1] >= lat_grid[1] );
      // The latitude in *lat2* must be <= the next to last value in *lat_grid*.
      assert ( lat_grid[lat2] <= lat_grid[lat_grid.nelem()-2] );
    }
    if ( atmosphere_dim == 3 )
    {
      // The longitude in *lon2* must be bigger than the longitude in *lon1*.
      assert ( lon_grid[lon2] > lon_grid[lon1] );
      // The longitude in *lon1* must be >= the second value in *lon_grid*.
      assert ( lon_grid[lon1] >= lon_grid[1] );
      // The longitude in *lon2* must be <= the next to last value in *lon_grid*.
      assert ( lon_grid[lon2] <= lon_grid[lon_grid.nelem()-2] );
    }
    */
  }

  else
  // if all particle fields are zero at each level and cloudbox was not preset,
  // switch cloudbox off.
  {
    CREATE_OUT0;
    cloudbox_on = 0;
    cloudbox_limits[1] =
        -1;  // just for consistency with cloudboxSetAutomatically
    out0 << "Cloudbox is switched off!\n";
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetFullAtm(  //WS Output
    Index& cloudbox_on,
    ArrayOfIndex& cloudbox_limits,
    // WS Input
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Verbosity&) {
  cloudbox_on = 1;
  cloudbox_limits.resize(2 * atmosphere_dim);

  cloudbox_limits[0] = 0;
  cloudbox_limits[1] = p_grid.nelem() - 1;

  if (atmosphere_dim > 1) {
    Index last_lat = lat_grid.nelem() - 1;

    // find minimum lat_grid point i with lat_grid[i]-lat_grid[0]>=LAT_LON_MIN
    Index i = 1;
    while ((i < last_lat - 1) && (lat_grid[i] - lat_grid[0] < LAT_LON_MIN)) i++;
    if (i == last_lat - 1) {
      ostringstream os;
      os << "Can not define lower latitude edge of cloudbox:\n"
         << "Extend of atmosphere too small. Distance to minimum latitude\n"
         << "has to be at least " << LAT_LON_MIN << "deg, but only "
         << lat_grid[i - 1] - lat_grid[0] << " available here.";
      throw runtime_error(os.str());
    }
    cloudbox_limits[2] = i;

    // find maximum lat_grid point j with lat_grid[-1]-lat_grid[j]>=LAT_LON_MIN
    // and j>i
    Index j = last_lat - 1;
    while ((j > i) && (lat_grid[last_lat] - lat_grid[j] < LAT_LON_MIN)) j--;
    if (j == i) {
      ostringstream os;
      os << "Can not define upper latitude edge of cloudbox:\n"
         << "Extend of atmosphere too small. Distance to maximum latitude\n"
         << "has to be at least " << LAT_LON_MIN << "deg, but only "
         << lat_grid[last_lat] - lat_grid[j + 1] << " available here.";
      throw runtime_error(os.str());
    }
    cloudbox_limits[3] = j;
  }

  if (atmosphere_dim > 2) {
    const Numeric latmax = max(abs(lat_grid[cloudbox_limits[2]]),
                               abs(lat_grid[cloudbox_limits[3]]));
    const Numeric lfac = 1 / cos(DEG2RAD * latmax);
    Index last_lon = lon_grid.nelem() - 1;

    // find minimum lon_grid point i with lon_grid[i]-lon_grid[0]>=LAT_LON_MIN/lfac
    Index i = 1;
    while ((i < last_lon - 1) &&
           (lon_grid[i] - lon_grid[0] < LAT_LON_MIN / lfac))
      i++;
    if (i == last_lon - 1) {
      ostringstream os;
      os << "Can not define lower longitude edge of cloudbox:\n"
         << "Extend of atmosphere too small. Distance to minimum longitude\n"
         << "has to be at least " << LAT_LON_MIN / lfac << "deg, but only "
         << lon_grid[i - 1] - lon_grid[0] << " available here.";
      throw runtime_error(os.str());
    }
    cloudbox_limits[4] = i;

    // find maximum lon_grid point j with lon_grid[-1]-lon_grid[j]>=LAT_LON_MIN/lfac
    // and j>i
    Index j = last_lon - 1;
    while ((j > i) && (lon_grid[last_lon] - lon_grid[j] < LAT_LON_MIN / lfac))
      j--;
    if (j == i) {
      ostringstream os;
      os << "Can not define upper longitude edge of cloudbox:\n"
         << "Extend of atmosphere too small. Distance to maximum longitude\n"
         << "has to be at least " << LAT_LON_MIN / lfac << "deg, but only "
         << lon_grid[last_lon] - lon_grid[j + 1] << " available here.";
      throw runtime_error(os.str());
    }
    cloudbox_limits[5] = j;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManually(  // WS Output:
    Index& cloudbox_on,
    ArrayOfIndex& cloudbox_limits,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // Control Parameters:
    const Numeric& p1,
    const Numeric& p2,
    const Numeric& lat1,
    const Numeric& lat2,
    const Numeric& lon1,
    const Numeric& lon2,
    const Verbosity&) {
  // Check existing WSV
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  // Check keyword arguments
  if (p1 <= p2)
    throw runtime_error(
        "The pressure in *p1* must be bigger than the "
        "pressure in *p2*.");
  if (p1 <= p_grid[p_grid.nelem() - 1])
    throw runtime_error(
        "The pressure in *p1* must be larger than the "
        "last value in *p_grid*.");
  if (p2 >= p_grid[0])
    throw runtime_error(
        "The pressure in *p2* must be smaller than the "
        "first value in *p_grid*.");
  if (atmosphere_dim >= 2) {
    if (lat2 <= lat1)
      throw runtime_error(
          "The latitude in *lat2* must be bigger than the "
          "latitude in *lat1*.");
    if (lat1 < lat_grid[1])
      throw runtime_error(
          "The latitude in *lat1* must be >= the "
          "second value in *lat_grid*.");
    if (lat2 > lat_grid[lat_grid.nelem() - 2])
      throw runtime_error(
          "The latitude in *lat2* must be <= the "
          "next to last value in *lat_grid*.");
  }
  if (atmosphere_dim == 3) {
    if (lon2 <= lon1)
      throw runtime_error(
          "The longitude in *lon2* must be bigger than the "
          "longitude in *lon1*.");
    if (lon1 < lon_grid[1])
      throw runtime_error(
          "The longitude in *lon1* must be >= the "
          "second value in *lon_grid*.");
    if (lon2 > lon_grid[lon_grid.nelem() - 2])
      throw runtime_error(
          "The longitude in *lon2* must be <= the "
          "next to last value in *lon_grid*.");
  }

  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize(atmosphere_dim * 2);

  // Pressure limits
  if (p1 > p_grid[1]) {
    cloudbox_limits[0] = 0;
  } else {
    for (cloudbox_limits[0] = 1; p_grid[cloudbox_limits[0] + 1] >= p1;
         cloudbox_limits[0]++) {
    }
  }
  if (p2 < p_grid[p_grid.nelem() - 2]) {
    cloudbox_limits[1] = p_grid.nelem() - 1;
  } else {
    for (cloudbox_limits[1] = p_grid.nelem() - 2;
         p_grid[cloudbox_limits[1] - 1] <= p2;
         cloudbox_limits[1]--) {
    }
  }

  // Latitude limits
  if (atmosphere_dim >= 2) {
    for (cloudbox_limits[2] = 1; lat_grid[cloudbox_limits[2] + 1] <= lat1;
         cloudbox_limits[2]++) {
    }
    for (cloudbox_limits[3] = lat_grid.nelem() - 2;
         lat_grid[cloudbox_limits[3] - 1] >= lat2;
         cloudbox_limits[3]--) {
    }
  }

  // Longitude limits
  if (atmosphere_dim == 3) {
    for (cloudbox_limits[4] = 1; lon_grid[cloudbox_limits[4] + 1] <= lon1;
         cloudbox_limits[4]++) {
    }
    for (cloudbox_limits[5] = lon_grid.nelem() - 2;
         lon_grid[cloudbox_limits[5] - 1] >= lon2;
         cloudbox_limits[5]--) {
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudboxSetManuallyAltitude(  // WS Output:
    Index& cloudbox_on,
    ArrayOfIndex& cloudbox_limits,
    // WS Input:
    const Index& atmosphere_dim,
    const Tensor3& z_field,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // Control Parameters:
    const Numeric& z1,
    const Numeric& z2,
    const Numeric& lat1,
    const Numeric& lat2,
    const Numeric& lon1,
    const Numeric& lon2,
    const Verbosity&) {
  // Check existing WSV
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);

  // Check keyword arguments
  if (z1 >= z2)
    throw runtime_error(
        "The altitude in *z1* must be smaller than the "
        "altitude in *z2*.");
  /* These cases are in fact handled
  if( z1 <= z_field(0, 0, 0) )
    throw runtime_error( "The altitude in *z1* must be larger than the "
                         "first value in *z_field*." );
  if( z2 >= z_field(z_field.npages()-1, 0, 0) )
    throw runtime_error( "The altitude in *z2* must be smaller than the "
                         "last value in *z_field*." );
  */
  if (atmosphere_dim == 3) {
    if (lat2 <= lat1)
      throw runtime_error(
          "The latitude in *lat2* must be bigger than the "
          " latitude in *lat1*.");
    if (lat1 < lat_grid[1])
      throw runtime_error(
          "The latitude in *lat1* must be >= the "
          "second value in *lat_grid*.");
    if (lat2 > lat_grid[lat_grid.nelem() - 2])
      throw runtime_error(
          "The latitude in *lat2* must be <= the "
          "next to last value in *lat_grid*.");
    if (lon2 <= lon1)
      throw runtime_error(
          "The longitude in *lon2* must be bigger than the "
          "longitude in *lon1*.");
    if (lon1 < lon_grid[1])
      throw runtime_error(
          "The longitude in *lon1* must be >= the "
          "second value in *lon_grid*.");
    if (lon2 > lon_grid[lon_grid.nelem() - 2])
      throw runtime_error(
          "The longitude in *lon2* must be <= the "
          "next to last value in *lon_grid*.");
  }

  // Set cloudbox_on
  cloudbox_on = 1;

  // Allocate cloudbox_limits
  cloudbox_limits.resize(atmosphere_dim * 2);

  // Pressure/altitude limits
  if (z1 < z_field(1, 0, 0)) {
    cloudbox_limits[0] = 0;
  } else {
    for (cloudbox_limits[0] = 1; z_field(cloudbox_limits[0] + 1, 0, 0) <= z1;
         cloudbox_limits[0]++) {
    }
  }
  if (z2 > z_field(z_field.npages() - 2, 0, 0)) {
    cloudbox_limits[1] = z_field.npages() - 1;
  } else {
    for (cloudbox_limits[1] = z_field.npages() - 2;
         z_field(cloudbox_limits[1] - 1, 0, 0) >= z2;
         cloudbox_limits[1]--) {
    }
  }

  // Latitude limits
  if (atmosphere_dim >= 2) {
    for (cloudbox_limits[2] = 1; lat_grid[cloudbox_limits[2] + 1] <= lat1;
         cloudbox_limits[2]++) {
    }
    for (cloudbox_limits[3] = lat_grid.nelem() - 2;
         lat_grid[cloudbox_limits[3] - 1] >= lat2;
         cloudbox_limits[3]--) {
    }
  }

  // Longitude limits
  if (atmosphere_dim == 3) {
    for (cloudbox_limits[4] = 1; lon_grid[cloudbox_limits[4] + 1] <= lon1;
         cloudbox_limits[4]++) {
    }
    for (cloudbox_limits[5] = lon_grid.nelem() - 2;
         lon_grid[cloudbox_limits[5] - 1] >= lon2;
         cloudbox_limits[5]--) {
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyInterpCloudboxField(Matrix& iy,
                           const Tensor7& doit_i_field,
                           const Vector& rte_pos,
                           const Vector& rte_los,
                           const Index& jacobian_do,
                           const Index& cloudbox_on,
                           const ArrayOfIndex& cloudbox_limits,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Tensor3& z_field,
                           const Matrix& z_surface,
                           const Index& stokes_dim,
                           const Vector& scat_za_grid,
                           const Vector& scat_aa_grid,
                           const Vector& f_grid,
                           //const Index&    p_interp_order,
                           const Index& za_interp_order,
                           const Index& za_restrict,
                           const Index& cos_za_interp,
                           const Numeric& za_extpolfac,
                           const Index& aa_interp_order,
                           const Verbosity&) {

  //--- Check input -----------------------------------------------------------
  if (!(atmosphere_dim == 1 || atmosphere_dim == 3))
    throw runtime_error("The atmospheric dimensionality must be 1 or 3.");
  if (jacobian_do)
    throw runtime_error(
        "This method does not support jacobians (jacobian_do must be 0)");
  if (!cloudbox_on)
    throw runtime_error(
        "The cloud box is not activated and no outgoing "
        "field can be returned.");
  if (cloudbox_limits.nelem() != 2 * atmosphere_dim)
    throw runtime_error(
        "*cloudbox_limits* is a vector which contains the upper and lower\n"
        "limit of the cloud for all atmospheric dimensions.\n"
        "So its length must be 2 x *atmosphere_dim*");
  if (scat_za_grid.nelem() == 0)
    throw runtime_error(
        "The variable *scat_za_grid* is empty. Are dummy "
        "values from *cloudboxOff used?");
  if (!(za_interp_order < scat_za_grid.nelem()))
    throw runtime_error(
        "Zenith angle interpolation order *za_interp_order*"
        " must be smaller\n"
        "than number of angles in *scat_za_grid*.");
  if (atmosphere_dim > 1 && !(aa_interp_order < scat_aa_grid.nelem()))
    throw runtime_error(
        "Azimuth angle interpolation order *aa_interp_order*"
        " must be smaller\n"
        "than number of angles in *scat_aa_grid*.");
  if (doit_i_field.nlibraries() != f_grid.nelem())
    throw runtime_error(
        "Inconsistency in size between f_grid and doit_i_field! "
        "(This method does not yet handle dispersion type calculations.)");
  //---------------------------------------------------------------------------

  // At each lat/lon, one value of the intensity field is defined at the
  // surface. Simplest way to handle this is to insert z_surface into z_field
  Tensor3 z_with_surface = z_field;
  for (Index ilat=0; ilat < z_surface.nrows(); ilat++) {
    for (Index ilon=0; ilon < z_surface.ncols(); ilon++) {
      Index ip = 0;
      for (; z_surface(ilat,ilon) >= z_field(ip+1,ilat,ilon); ip++) {}
      z_with_surface(ip,ilat,ilon) = z_surface(ilat,ilon);
    }    
  }
  
  // Convert rte_pos to grid positions
  GridPos gp_p, gp_lat, gp_lon;
  rte_pos2gridpos(gp_p,
                  gp_lat,
                  gp_lon,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid,
                  z_with_surface,
                  rte_pos);

  //--- Determine if at border or inside of cloudbox (or outside!)
  //
  // Let us introduce a number coding for cloud box borders.
  // Borders have the same number as position in *cloudbox_limits*.
  // Inside cloud box is coded as 99, and outside as > 100.
  Index border = 999;
  //
  //- Check if at any border
  if (is_gridpos_at_index_i(gp_p, cloudbox_limits[0], false)) {
    border = 0;
  } else if (is_gridpos_at_index_i(gp_p, cloudbox_limits[1], false)) {
    border = 1;
  } else if (atmosphere_dim > 1) {
    if (is_gridpos_at_index_i(gp_lat, cloudbox_limits[2], false)) {
      border = 2;
    } else if (is_gridpos_at_index_i(gp_lat, cloudbox_limits[3], false)) {
      border = 3;
    } else if (atmosphere_dim > 2) {
      if (is_gridpos_at_index_i(gp_lon, cloudbox_limits[4], false)) {
        border = 4;
      } else if (is_gridpos_at_index_i(gp_lon, cloudbox_limits[5], false)) {
        border = 5;
      }
    }
  }

  //
  //- Check if inside (when border<100 here, it means we are on a cloudbox border)
  if (border > 100) {
    // Assume inside as it is easiest to detect if outside (can be detected
    // check in one dimension at the time)
    bool inside = true;
    Numeric fgp;

    // Check in pressure dimension
    fgp = fractional_gp(gp_p);
    if (fgp < Numeric(cloudbox_limits[0]) ||
        fgp > Numeric(cloudbox_limits[1])) {
      inside = false;
    }

    // Check in lat and lon dimensions
    if (atmosphere_dim == 3 && inside) {
      fgp = fractional_gp(gp_lat);
      if (fgp < Numeric(cloudbox_limits[2]) ||
          fgp > Numeric(cloudbox_limits[3])) {
        inside = false;
      }
      fgp = fractional_gp(gp_lon);
      if (fgp < Numeric(cloudbox_limits[4]) ||
          fgp > Numeric(cloudbox_limits[5])) {
        inside = false;
      }
    }

    if (inside) {
      border = 99;
    }
  }

  // If outside, something is wrong
  if (border > 100) {
    throw runtime_error(
        "Given position has been found to be outside the cloudbox.");
  }

  //- Sizes
  const Index nf = f_grid.nelem();
  DEBUG_ONLY(const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1);
  const Index nza = scat_za_grid.nelem();
  const Index naa = doit_i_field.nrows();

  //- Resize *iy*
  iy.resize(nf, stokes_dim);
  //- Make space for spatially interpolated field (we perfrom angle
  //interpolation separately on this then)
  Tensor4 i_field_local(nf, nza, naa, stokes_dim);

  // Index of the p/lat/lon slice (first or last element slice in this
  // dimension), where rte_pos is located. if located on a border.
  Index border_index;
  if (border < 99) {
    if (border % 2)  // odd border id, ie top or north or east
      border_index = cloudbox_limits[border] - cloudbox_limits[border - 1];
    else  // even border id, ie bottom or south or west
      border_index = 0;
  }

  // Sensor points inside the cloudbox
  if (border == 99) {
    if (atmosphere_dim == 3) {
      throw runtime_error(
          "Radiation extraction for a position inside cloudbox\n"
          "is not yet implemented for 3D cases.\n");
    } else {
      assert(atmosphere_dim == 1);

      assert(is_size(doit_i_field, nf, np, 1, 1, nza, 1, stokes_dim));

      // Grid position in *p_grid* 
      gp_p.idx = gp_p.idx - cloudbox_limits[0];
      gridpos_upperend_check(gp_p, cloudbox_limits[1] - cloudbox_limits[0]);
      Vector itw_p(2);
      interpweights(itw_p, gp_p);

      for (Index is = 0; is < stokes_dim; is++)
        for (Index iv = 0; iv < nf; iv++)
          for (Index i_za = 0; i_za < nza; i_za++)
            i_field_local(iv, i_za, 0, is) =
                interp(itw_p, doit_i_field(iv, joker, 0, 0, i_za, 0, is), gp_p);
    }
  }

  // Sensor outside the cloudbox

  // --- 1D ------------------------------------------------------------------
  else if (atmosphere_dim == 1) {
    i_field_local =
        doit_i_field(joker, border_index, 0, 0, joker, joker, joker);
  }

  // --- 3D ------------------------------------------------------------------
  else {
    assert(is_size(doit_i_field,
                   nf,
                   doit_i_field.nvitrines(),
                   doit_i_field.nshelves(),
                   doit_i_field.nbooks(),
                   scat_za_grid.nelem(),
                   scat_aa_grid.nelem(),
                   stokes_dim));

    // Interpolation weights (for 2D "red" interpolation)
    Vector itw(4);

    // Outgoing from cloudbox top or bottom, i.e. from a pressure level
    if (border <= 1) {
      // Lat and lon grid positions with respect to cloud box
      GridPos cb_gp_lat, cb_gp_lon;
      cb_gp_lat = gp_lat;
      cb_gp_lon = gp_lon;
      cb_gp_lat.idx -= cloudbox_limits[2];
      cb_gp_lon.idx -= cloudbox_limits[4];
      //
      gridpos_upperend_check(cb_gp_lat,
                             cloudbox_limits[3] - cloudbox_limits[2]);
      gridpos_upperend_check(cb_gp_lon,
                             cloudbox_limits[5] - cloudbox_limits[4]);

      interpweights(itw, cb_gp_lat, cb_gp_lon);

      for (Index is = 0; is < stokes_dim; is++)
        for (Index iv = 0; iv < nf; iv++)
          for (Index i_za = 0; i_za < nza; i_za++)
            for (Index i_aa = 0; i_aa < naa; i_aa++)
              i_field_local(iv, i_za, i_aa, is) = interp(
                  itw,
                  doit_i_field(iv, border_index, joker, joker, i_za, i_aa, is),
                  cb_gp_lat,
                  cb_gp_lon);
    }

    // Outgoing from cloudbox north or south border, i.e. from a latitude level
    else if (border <= 3) {
      // Pressure and lon grid positions with respect to cloud box
      GridPos cb_gp_p, cb_gp_lon;
      cb_gp_p = gp_p;
      cb_gp_lon = gp_lon;
      cb_gp_p.idx -= cloudbox_limits[0];
      cb_gp_lon.idx -= cloudbox_limits[4];
      //
      gridpos_upperend_check(cb_gp_p, cloudbox_limits[1] - cloudbox_limits[0]);
      gridpos_upperend_check(cb_gp_lon,
                             cloudbox_limits[5] - cloudbox_limits[4]);

      interpweights(itw, cb_gp_p, cb_gp_lon);

      for (Index is = 0; is < stokes_dim; is++)
        for (Index iv = 0; iv < nf; iv++)
          for (Index i_za = 0; i_za < nza; i_za++)
            for (Index i_aa = 0; i_aa < naa; i_aa++)
              i_field_local(iv, i_za, i_aa, is) = interp(
                  itw,
                  doit_i_field(iv, joker, border_index, joker, i_za, i_aa, is),
                  cb_gp_p,
                  cb_gp_lon);
    }

    // Outgoing from cloudbox east or west border, i.e. from a longitude level
    else {
      // Pressure and lat grid positions with respect to cloud box
      GridPos cb_gp_p, cb_gp_lat;
      cb_gp_p = gp_p;
      cb_gp_lat = gp_lat;
      cb_gp_p.idx -= cloudbox_limits[0];
      cb_gp_lat.idx -= cloudbox_limits[2];
      //
      gridpos_upperend_check(cb_gp_p, cloudbox_limits[1] - cloudbox_limits[0]);
      gridpos_upperend_check(cb_gp_lat,
                             cloudbox_limits[3] - cloudbox_limits[2]);

      interpweights(itw, cb_gp_p, cb_gp_lat);

      for (Index is = 0; is < stokes_dim; is++)
        for (Index iv = 0; iv < nf; iv++)
          for (Index i_za = 0; i_za < nza; i_za++)
            for (Index i_aa = 0; i_aa < naa; i_aa++)
              i_field_local(iv, i_za, i_aa, is) = interp(
                  itw,
                  doit_i_field(iv, joker, joker, border_index, i_za, i_aa, is),
                  cb_gp_p,
                  cb_gp_lat);
    }
  }

  // now, do Nth-oder polynomial interpoaltion of angle(s)
  //
  // a bunch of options for testing:
  //   a) interpolation over full zenith angle grid vs over hemispheres
  //      separately.
  //   b) interpolation in plain zenith angles vs. in cosines of zenith angle.


  // find range of scat_za_grid that we will do interpolation over.
  Index za_start = 0;
  Index za_extend = scat_za_grid.nelem();
  if (za_restrict) {
    // which hemisphere do we need?
    if (is_same_within_epsilon(rte_los[0], 90., 1e-6)) {
      //FIXME: we should allow this though, if scat_za_grid has a grid point
      //at 90deg.
      throw runtime_error(
          "Hemisphere-restricted zenith angle interpolation not allowed\n"
          "for 90degree views.");
    }
    else if (rte_los[0] > 90) {
      // upwelling, i.e. second part of scat_za_grid. that is, we need to find
      // the first point in scat_za_grid where za>90. and update za_start
      // accordingly.
      while (za_start < scat_za_grid.nelem() && scat_za_grid[za_start] < 90.) {
        za_start++;
      }
      if (za_start == scat_za_grid.nelem())
        throw runtime_error(
            "No scat_za_grid grid point found in 90-180deg hemisphere.\n"
            "No hemispheric interpolation possible.");
      za_extend -= za_start;
    } else {
      // downwelling, i.e. first part of scat_za_grid. that is, we need to
      // find the last point in scat_za_grid where za<90. and update za_extend
      // accordingly.
      while (za_extend > 0 && scat_za_grid[za_extend - 1] > 90.) {
        za_extend--;
      }
      if (za_extend == 0)
        throw runtime_error(
            "No scat_za_grid grid point found in 0-90deg hemisphere.\n"
            "No hemispheric interpolation possible.");
    }
    if (!(za_interp_order < za_extend)) {
      ostringstream os;
      os << "Zenith angle interpolation order *za_interp_order* ("
         << za_interp_order << ") must be smaller\n"
         << "than number of angles in respective hemisphere (" << za_extend
         << ").";
      throw runtime_error(os.str());
    }
  }
  
  // Grid position in *scat_za_grid*
  GridPosPoly gp_za, gp_aa;

  if (cos_za_interp) {
    Vector cosza_grid(za_extend);
    const Numeric cosza = cos(rte_los[0] * DEG2RAD);

    for (Index i_za = 0; i_za < za_extend; i_za++)
      cosza_grid[i_za] = cos(scat_za_grid[i_za + za_start] * DEG2RAD);

    // OBS: cosza is a decreasing grid!
    const Numeric cosza_min =
        cosza_grid[0] + za_extpolfac * (cosza_grid[0] - cosza_grid[1]);
    if (cosza > cosza_min) {
      ostringstream os;
      os << "Zenith angle " << rte_los[0] << "deg is outside the range"
         << " covered by scat_za_grid.\n"
         << "Lower limit of allowed range is " << acos(cosza_min) * RAD2DEG
         << ".\n"
         << "Increase za_extpolfac (now=" << za_extpolfac << ") or"
         << " use wider scat_za_grid.\n";
      throw runtime_error(os.str());
    }
    const Numeric cosza_max =
        cosza_grid[za_extend - 1] -
        za_extpolfac * (cosza_grid[za_extend - 2] - cosza_grid[za_extend - 1]);
    if (cosza < cosza_max) {
      ostringstream os;
      os << "Zenith angle " << rte_los[0] << "deg is outside the range"
         << " covered by scat_za_grid.\n"
         << "Upper limit of allowed range is " << acos(cosza_max) * RAD2DEG
         << ".\n"
         << "Increase za_extpolfac (now=" << za_extpolfac << ") or"
         << " use wider scat_za_grid.\n";
      throw runtime_error(os.str());
    }

    gridpos_poly(gp_za, cosza_grid, cosza, za_interp_order, za_extpolfac);
  } else {
    Vector za_grid = scat_za_grid[Range(za_start, za_extend)];
    const Numeric za = rte_los[0];

    const Numeric za_min =
        za_grid[0] - za_extpolfac * (za_grid[1] - za_grid[0]);
    if (za < za_min) {
      ostringstream os;
      os << "Zenith angle " << rte_los[0] << "deg is outside the range"
         << " covered by scat_za_grid.\n"
         << "Lower limit of allowed range is " << za_min << ".\n"
         << "Increase za_extpolfac (now=" << za_extpolfac << ") or"
         << " use wider scat_za_grid.\n";
      throw runtime_error(os.str());
    }
    const Numeric za_max =
        za_grid[za_grid.nelem() - 1] +
        za_extpolfac * (za_grid[za_extend - 1] - za_grid[za_extend - 2]);
    if (za > za_max) {
      ostringstream os;
      os << "Zenith angle " << za << "deg is outside the range"
         << " covered by scat_za_grid.\n"
         << "Upper limit of allowed range is " << za_max << ".\n"
         << "Increase za_extpolfac (now=" << za_extpolfac << ") or"
         << " use wider scat_za_grid.\n";
      throw runtime_error(os.str());
    }
    gridpos_poly(gp_za, za_grid, za, za_interp_order, za_extpolfac);
  }

  if (atmosphere_dim > 1) {
    gridpos_poly_cyclic_longitudinal(
        gp_aa, scat_aa_grid, rte_los[1], aa_interp_order);
  } else {
    gp_aa.idx.resize(1);
    gp_aa.idx[0] = 0;
    gp_aa.w.resize(1);
    gp_aa.w[0] = 1;
  }

  // Corresponding interpolation weights
  Vector itw_angs(gp_za.idx.nelem() * gp_aa.idx.nelem());
  interpweights(itw_angs, gp_za, gp_aa);

  for (Index is = 0; is < stokes_dim; is++) {
    for (Index iv = 0; iv < nf; iv++) {
      iy(iv, is) =
          interp(itw_angs,
                 i_field_local(iv, Range(za_start, za_extend), joker, is),
                 gp_za,
                 gp_aa);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_fieldCleanup(  //WS Output:
    Tensor4& particle_field_out,
    //WS Input:
    const Tensor4& particle_field_in,
    const Numeric& threshold,
    const Verbosity&) {
  if (&particle_field_out != &particle_field_in) {
    particle_field_out = particle_field_in;
  }

  // Check that particle_field contains realistic values
  //(values smaller than the threshold will be set to 0)
  for (Index i = 0; i < particle_field_out.nbooks(); i++) {
    for (Index j = 0; j < particle_field_out.npages(); j++) {
      for (Index k = 0; k < particle_field_out.nrows(); k++) {
        for (Index l = 0; l < particle_field_out.ncols(); l++) {
          if (particle_field_out(i, j, k, l) < threshold) {
            particle_field_out(i, j, k, l) = 0.0;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesInit(  //WS Output:
    ArrayOfString& scat_species,
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    Index& scat_data_checked,
    ArrayOfGriddedField3& pnd_field_raw,
    const Verbosity&) {
  scat_species.resize(0);
  scat_data_raw.resize(0);
  scat_meta.resize(0);
  pnd_field_raw.resize(0);
  scat_data_checked = 0;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsPndAndScatAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfGriddedField3& pnd_field_raw,
    // WS Input (needed for checking the datafiles):
    const Index& atmosphere_dim,
    // Keywords:
    const ArrayOfString& scat_data_files,
    const ArrayOfString& pnd_field_files,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  //--- Reading the data ---------------------------------------------------

  if (scat_data_files.nelem() != pnd_field_files.nelem()) {
    ostringstream os;
    os << "Number of elements in scat_data and pnd_field filelists is"
       << "inconsistent.";
    throw runtime_error(os.str());
  }

  Index last_species = scat_data_raw.nelem() - 1;
  if (last_species == -1) {
    scat_data_raw.resize(1);
    last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data_raw* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;

  for (Index i = 0; i < scat_data_files.nelem(); i++) {
    // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors.
    scat_data_raw[last_species].push_back(scat_data_single);
    pnd_field_raw.push_back(pnd_field_data);

    out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
    xml_read_from_file(
        scat_data_files[i],
        scat_data_raw[last_species][scat_data_raw[last_species].nelem() - 1],
        verbosity);

    out2 << "  Read particle number density field\n";
    if (pnd_field_files[i].nelem() < 1) {
      CREATE_OUT1;
      out1 << "Warning: No pnd_field_file specified. Ignored here,\n"
           << "but user HAS TO add that later on!. \n";
    } else {
      xml_read_from_file(pnd_field_files[i],
                         pnd_field_raw[pnd_field_raw.nelem() - 1],
                         verbosity);

      chk_pnd_data(pnd_field_raw[pnd_field_raw.nelem() - 1],
                   pnd_field_files[i],
                   atmosphere_dim,
                   verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesPndAndScatAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfGriddedField3& pnd_field_raw,
    // WS Input(needed for checking the datafiles):
    const Index& atmosphere_dim,
    // Keywords:
    const ArrayOfString& scat_data_files,
    const String& pnd_fieldarray_file,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  //chk_atm_grids ( atmosphere_dim, p_grid, lat_grid, lon_grid );

  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  arr_ssd.resize(scat_data_files.nelem());

  for (Index i = 0; i < scat_data_files.nelem(); i++) {
    out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
    xml_read_from_file(scat_data_files[i], arr_ssd[i], verbosity);
  }

  // append as new scattering species
  if (scat_data_raw.nelem() == 0) {
    scat_data_raw.resize(1);
    scat_data_raw[0] = arr_ssd;
  } else
    scat_data_raw.push_back(arr_ssd);

  out2 << "  Read particle number density data \n";
  ArrayOfGriddedField3 pnd_tmp;
  xml_read_from_file(pnd_fieldarray_file, pnd_tmp, verbosity);

  chk_pnd_raw_data(pnd_tmp, pnd_fieldarray_file, atmosphere_dim, verbosity);

  // append to pnd_field_raw
  for (Index i = 0; i < pnd_tmp.nelem(); ++i)
    pnd_field_raw.push_back(pnd_tmp[i]);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsToabs_speciesAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfGriddedField3& vmr_field_raw,
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    Index& abs_xsec_agenda_checked,
    // WS Input (needed for checking the datafiles):
    const Index& atmosphere_dim,
    const Vector& f_grid,
    // Keywords:
    const ArrayOfString& scat_data_files,
    const ArrayOfString& pnd_field_files,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  //--- Check input ---------------------------------------------------------

  // Atmosphere
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  //chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Frequency grid
  if (f_grid.empty()) throw runtime_error("The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  //--- Reading the data ---------------------------------------------------

  if (scat_data_files.nelem() != pnd_field_files.nelem()) {
    ostringstream os;
    os << "Number of elements in scat_data and pnd_field filelists is"
       << "inconsistent.";
    throw runtime_error(os.str());
  }

  Index last_species = scat_data_raw.nelem() - 1;
  if (last_species == -1) {
    scat_data_raw.resize(1);
    last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data_raw* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;
  ArrayOfString species(1);
  species[0] = "particles";

  for (Index i = 0; i < scat_data_files.nelem(); i++) {
    // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors.
    scat_data_raw[last_species].push_back(scat_data_single);
    vmr_field_raw.push_back(pnd_field_data);

    out2 << "  Read single scattering data file " << scat_data_files[i] << "\n";
    xml_read_from_file(
        scat_data_files[i],
        scat_data_raw[last_species][scat_data_raw[last_species].nelem() - 1],
        verbosity);

    out2 << "  Check single scattering properties\n";
    chk_interpolation_grids(
        "scat_data_single.f_grid to f_grid",
        scat_data_raw[last_species][scat_data_raw[last_species].nelem() - 1]
            .f_grid,
        f_grid);

    out2 << "  Read particle number density field\n";
    if (pnd_field_files[i].nelem() < 1) {
      CREATE_OUT1;
      out1 << "Warning: No pnd_field_file specified. Ignored here,\n"
           << "but user HAS TO add that later on!\n";
    } else {
      try {
        xml_read_from_file(pnd_field_files[i],
                           vmr_field_raw[vmr_field_raw.nelem() - 1],
                           verbosity);
      } catch (...) {
        ArrayOfGriddedField3 tmp;
        try {
          xml_read_from_file(pnd_field_files[i], tmp, verbosity);
          if (tmp.nelem() == 1) {
            vmr_field_raw[vmr_field_raw.nelem() - 1] = tmp[0];
          } else {
            std::ostringstream os;
            os << "The file " << pnd_field_files[i] << "\n"
               << "is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n";
            throw std::runtime_error(os.str());
          }
        } catch (...) {
          std::ostringstream os;
          os << "The file " << pnd_field_files[i] << " does not exist or\n"
             << "its type is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n";
          throw std::runtime_error(os.str());
        }
      }

      chk_pnd_data(vmr_field_raw[vmr_field_raw.nelem() - 1],
                   pnd_field_files[i],
                   atmosphere_dim,
                   verbosity);
    }

    out2 << "  Append 'particle' field to abs_species\n";
    abs_speciesAdd(abs_species,
                   propmat_clearsky_agenda_checked,
                   abs_xsec_agenda_checked,
                   species,
                   verbosity);
  }
  scat_dataCheck(scat_data_raw, "sane", 1e-2, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesScatAndMetaRead(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    // Keywords:
    const ArrayOfString& scat_data_files,
    const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  ArrayOfScatteringMetaData arr_smd;

  arr_ssd.resize(scat_data_files.nelem());
  arr_smd.resize(scat_data_files.nelem());

  Index meta_naming_conv = 0;

  for (Index i = 0; i < 1 && i < scat_data_files.nelem(); i++) {
    out3 << "  Read single scattering data file " << scat_data_files[i] << "\n";
    xml_read_from_file(scat_data_files[i], arr_ssd[i], verbosity);

    // make meta data name from scat data name
    ArrayOfString strarr;
    String scat_meta_file;

    if (i == 0) {
      scat_data_files[i].split(strarr, ".xml");
      scat_meta_file = strarr[0] + ".meta.xml";

      try {
        find_xml_file(scat_meta_file, verbosity);
      } catch (const runtime_error&) {
      }

      if (file_exists(scat_meta_file)) {
        out3 << "  Read scattering meta data\n";

        xml_read_from_file(scat_meta_file, arr_smd[i], verbosity);

        meta_naming_conv = 1;
      } else {
        try {
          scat_data_files[i].split(strarr, "scat_data");
          if (strarr.nelem() < 2) {
            throw std::runtime_error(
                "Splitting scattering data filename up at 'scat_data' also failed.");
          }
          scat_meta_file = strarr[0] + "scat_meta" + strarr[1];

          out3 << "  Read scattering meta data\n";
          xml_read_from_file(scat_meta_file, arr_smd[i], verbosity);

          meta_naming_conv = 2;
        } catch (const std::runtime_error& e) {
          ostringstream os;
          os << "No meta data file following one of the allowed naming "
             << "conventions was found.\n"
             << "Allowed are "
             << "*.meta.xml from *.xml and "
             << "*scat_meta* from *scat_data*\n"
             << "Scattering meta data file not found: " << scat_meta_file
             << "\n"
             << e.what();
          throw runtime_error(os.str());
        }
      }
    }
  }

  ArrayOfString fail_msg;

#pragma omp parallel for if (!arts_omp_in_parallel() &&                       \
                             scat_data_files.nelem() > 1)                     \
    num_threads(arts_omp_get_max_threads() > 16 ? 16                          \
                                                : arts_omp_get_max_threads()) \
        shared(out3, scat_data_files, arr_ssd, arr_smd)
  for (Index i = 1; i < scat_data_files.nelem(); i++) {
    // make meta data name from scat data name
    ArrayOfString strarr;
    String scat_meta_file;
    SingleScatteringData ssd;
    ScatteringMetaData smd;

    try {
      out3 << "  Read single scattering data file " << scat_data_files[i]
           << "\n";
      xml_read_from_file(scat_data_files[i], ssd, verbosity);

      scat_data_files[i].split(strarr, ".xml");
      scat_meta_file = strarr[0] + ".meta.xml";

      if (meta_naming_conv == 1) {
        scat_data_files[i].split(strarr, ".xml");
        scat_meta_file = strarr[0] + ".meta.xml";

        out3 << "  Read scattering meta data\n";
        xml_read_from_file(scat_meta_file, smd, verbosity);
      } else {
        scat_data_files[i].split(strarr, "scat_data");
        scat_meta_file = strarr[0] + "scat_meta" + strarr[1];

        out3 << "  Read scattering meta data\n";
        xml_read_from_file(scat_meta_file, smd, verbosity);
      }
    } catch (const std::exception& e) {
      ostringstream os;
      os << "Run-time error reading scattering data : \n" << e.what();
#pragma omp critical(ybatchCalc_push_fail_msg)
      fail_msg.push_back(os.str());
    }

    //FIXME: currently nothing is done in chk_scattering_meta_data!
    chk_scattering_meta_data(smd, scat_meta_file, verbosity);

#pragma omp critical(ScatSpeciesScatAndMetaRead_assign_ssd)
    arr_ssd[i] = std::move(ssd);
#pragma omp critical(ScatSpeciesScatAndMetaRead_assign_smd)
    arr_smd[i] = std::move(smd);
  }

  if (fail_msg.nelem()) {
    ostringstream os;
    for (auto& msg : fail_msg) os << msg << '\n';

    throw runtime_error(os.str());
  }

  // check if arrays have same size
  chk_scattering_data(arr_ssd, arr_smd, verbosity);

  // append as new scattering species
  scat_data_raw.push_back(std::move(arr_ssd));
  scat_meta.push_back(std::move(arr_smd));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsSelect(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    // WS Input:
    const ArrayOfString& scat_species,
    const String& species,
    const String& sizeparam,
    const Numeric& sizemin,
    const Numeric& sizemax,
    const Numeric& tolerance,
    const String& delim,
    const Verbosity&) {
  // first check that sizes of scat_species and scat_data_raw/scat_meta agree
  Index nspecies = scat_species.nelem();
  if (nspecies != scat_data_raw.nelem() || nspecies != scat_meta.nelem()) {
    ostringstream os;
    os << "Number of scattering species specified by scat_species does\n"
       << "not agree with number of scattering species in\n"
       << "scat_data_raw or scat_meta:\n"
       << "scat_species has " << nspecies
       << " entries, while scat_data_raw has " << scat_data_raw.nelem()
       << " and scat_meta has " << scat_meta.nelem() << ".";
    throw runtime_error(os.str());
  }

  // create temporary containers for selected elements
  ArrayOfSingleScatteringData scat_data_raw_tmp;
  ArrayOfScatteringMetaData scat_meta_tmp;

  String partfield_name;
  //find the species to handle: compare 'species' to 'partfield' part of
  //scat_species tags
  Index i_ss = -1;
  for (Index i = 0; i < scat_species.nelem(); i++) {
    parse_partfield_name(partfield_name, scat_species[i], delim);
    if (partfield_name == species) i_ss = i;
  }
  if (i_ss < 0) {
    ostringstream os;
    os << "Scattering species " << species << " not found among scat_species.";
    throw runtime_error(os.str());
  }

  // choosing the specified SingleScatteringData and ScatteringMetaData
  if (sizeparam == "diameter_max")
    for (Index i_se = 0; i_se < scat_meta[i_ss].nelem(); i_se++) {
      // scattering element diameter is extracted from the
      // scattering element's meta data and checked whether it's within size
      // selected range (sizemax < 0 check follows from wildcard usage and
      // means consider all sizes on the upper end)
      if (scat_meta[i_ss][i_se].diameter_max > sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance > scat_meta[i_ss][i_se].diameter_max ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else if (sizeparam == "diameter_volume_equ")
    for (Index i_se = 0; i_se < scat_meta[i_ss].nelem(); i_se++) {
      if (scat_meta[i_ss][i_se].diameter_volume_equ >
              sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance >
               scat_meta[i_ss][i_se].diameter_volume_equ ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else if (sizeparam == "diameter_area_equ_aerodynamical")
    for (Index i_se = 0; i_se < scat_meta[i_ss].nelem(); i_se++) {
      if (scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical >
              sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance >
               scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else {
    ostringstream os;
    os << "Size parameter " << sizeparam << "is unknown.";
    throw runtime_error(os.str());
  }

  // To use a particle species field without associated scattering element
  // data poses a high risk of accidentially neglecting these species. That's
  // unlikely what the user intends. Hence throw error.
  if (scat_meta_tmp.nelem() < 1) {
    ostringstream os;
    os << "For scattering species " << species << " no scattering "
       << "element matching the requested size range found.\n"
       << "Check *scat_data_raw* and *scat_meta* input as well as your size limit "
       << "selection!";
    throw runtime_error(os.str());
  }

  scat_meta[i_ss] = std::move(scat_meta_tmp);
  scat_data_raw[i_ss] = std::move(scat_data_raw_tmp);

  // check if array is empty. should never apply (since we checked the re-worked
  // data before and that error should also catch cases that are empty from the
  // beginning).
  assert(TotalNumberOfElements(scat_meta));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesExtendTemperature(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    // Keywords:
    const ArrayOfString& scat_species,
    const String& species,
    const String& scat_species_delim,
    const Numeric& T_low,
    const Numeric& T_high,
    const Verbosity&) {
  const bool do_tl = (T_low >= 0.);
  const bool do_th = (T_high >= 0.);

  if (do_tl || do_th) {
    Index i_ss = -1;
    if (species == "") {
      i_ss = scat_data_raw.nelem() - 1;
      if (i_ss == -1) {
        ostringstream os;
        os << "No *scat_data* available. Can not extend temperature range on "
           << "inexistent data.";
        throw runtime_error(os.str());
      }
    } else {
      // first check that sizes of scat_species and scat_data_raw agree
      Index nspecies = scat_species.nelem();
      if (nspecies != scat_data_raw.nelem()) {
        ostringstream os;
        os << "Number of scattering species specified by scat_species does\n"
           << "not agree with number of scattering species in *scat_data*:\n"
           << "scat_species has " << nspecies
           << " entries, while *scat_data* has " << scat_data_raw.nelem()
           << ".";
        throw runtime_error(os.str());
      }
      String partfield_name;
      //find the species to handle: compare 'species' to 'partfield' part of
      //scat_species tags
      for (Index i = 0; i < scat_species.nelem(); i++) {
        parse_partfield_name(
            partfield_name, scat_species[i], scat_species_delim);
        if (partfield_name == species) i_ss = i;
      }
      if (i_ss < 0) {
        ostringstream os;
        os << "Scattering species " << species
           << " not found among scat_species.";
        throw runtime_error(os.str());
      }
    }

    for (Index i_se = 0; i_se < scat_data_raw[i_ss].nelem(); i_se++) {
      const SingleScatteringData& ssdo = scat_data_raw[i_ss][i_se];
      const Index nTo = ssdo.T_grid.nelem();
      Index nTn = nTo;
      bool do_htl, do_hth;
      if (nTo > 1) {
        do_htl = (do_tl && (T_low < ssdo.T_grid[0]));
        do_hth = (do_th && (T_high > last(ssdo.T_grid)));
      } else {
        do_htl = false;
        do_hth = false;
      }

      if (do_htl || do_hth) {
        // create new instance of SingleScatteringData
        SingleScatteringData ssdn;
        Index iToff;

        // determine new temperature grid
        if (do_htl) nTn += 1;
        if (do_hth) nTn += 1;
        Vector T_grid_new(nTn);
        if (do_htl) {
          T_grid_new[0] = T_low;
          iToff = 1;
        } else {
          iToff = 0;
        }
        for (Index iT = 0; iT < nTo; iT++)
          T_grid_new[iT + iToff] = scat_data_raw[i_ss][i_se].T_grid[iT];
        if (do_hth) T_grid_new[nTo + iToff] = T_high;
        ssdn.T_grid = std::move(T_grid_new);

        // copy grids and other descriptive data that is to remain identical
        ssdn.ptype = ssdo.ptype;
        ostringstream description;
        description << ssdo.description;  // here just copy. we append further
                                          // info below if applicable.
        ssdn.f_grid = ssdo.f_grid;
        ssdn.za_grid = ssdo.za_grid;
        ssdn.aa_grid = ssdo.aa_grid;

        // determine size of current optical property data
        const Index nf = ssdo.f_grid.nelem();
        const Index nzas = ssdo.pha_mat_data.nshelves();
        const Index naas = ssdo.pha_mat_data.nbooks();
        const Index nzai = ssdo.pha_mat_data.npages();
        const Index naai = ssdo.pha_mat_data.nrows();
        const Index nmep = ssdo.pha_mat_data.ncols();
        const Index nmee = ssdo.ext_mat_data.ncols();
        const Index nvea = ssdo.abs_vec_data.ncols();

        // create containers for extended optical property data
        ssdn.pha_mat_data.resize(nf, nTn, nzas, naas, nzai, naai, nmep);
        ssdn.ext_mat_data.resize(nf, nTn, nzai, naai, nmee);
        ssdn.abs_vec_data.resize(nf, nTn, nzai, naai, nvea);

        // copy optical property data
        for (Index iT = 0; iT < nTo; iT++) {
          ssdn.pha_mat_data(
              joker, iT + iToff, joker, joker, joker, joker, joker) =
              ssdo.pha_mat_data(joker, iT, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, iT + iToff, joker, joker, joker) =
              ssdo.ext_mat_data(joker, iT, joker, joker, joker);
          ssdn.abs_vec_data(joker, iT + iToff, joker, joker, joker) =
              ssdo.abs_vec_data(joker, iT, joker, joker, joker);
        }

        // duplicate optical property data on T-edges if applicable
        if (do_htl) {
          ssdn.pha_mat_data(joker, 0, joker, joker, joker, joker, joker) =
              ssdn.pha_mat_data(joker, 1, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, 0, joker, joker, joker) =
              ssdn.ext_mat_data(joker, 1, joker, joker, joker);
          ssdn.abs_vec_data(joker, 0, joker, joker, joker) =
              ssdn.abs_vec_data(joker, 1, joker, joker, joker);
          description << "\n"
                      << "Low temperature limit extended by"
                      << " duplicating previous low temperature limit"
                      << " single scattering properties.";
        }
        if (do_hth) {
          ssdn.pha_mat_data(joker, nTn - 1, joker, joker, joker, joker, joker) =
              ssdn.pha_mat_data(
                  joker, nTn - 2, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, nTn - 1, joker, joker, joker) =
              ssdn.ext_mat_data(joker, nTn - 2, joker, joker, joker);
          ssdn.abs_vec_data(joker, nTn - 1, joker, joker, joker) =
              ssdn.abs_vec_data(joker, nTn - 2, joker, joker, joker);
          description << "\n"
                      << "High temperature limit extended by"
                      << " duplicating previous high temperature limit"
                      << " single scattering properties.";
        }
        ssdn.description = description.str();
        scat_data_raw[i_ss][i_se] = std::move(ssdn);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFrompnd_field_raw(  //WS Output:
    Tensor4& pnd_field,
    ArrayOfTensor4& dpnd_field_dx,
    //WS Input
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const ArrayOfGriddedField3& pnd_field_raw,
    const Index& atmosphere_dim,
    const ArrayOfIndex& cloudbox_limits,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Index& zeropadding,
    const Verbosity& verbosity) {
  // Basic checks of input variables
  //
  // Particle number density data
  //
  if (pnd_field_raw.empty()) {
    ostringstream os;
    os << "No particle number density data given. Please use WSMs\n"
       << "*ParticleTypeInit* and *ParticleTypeAdd(All)* for reading\n"
       << "scattering element data.\n";
    throw runtime_error(os.str());
  }

  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);
  ArrayOfIndex cloudbox_limits_tmp;
  /*if ( cloudbox_limits.empty() )
    {
      //If no limits set, the cloud(box) is supposed to cover the
      //complete atmosphere. This particularly to facilitate use of
      //scat_data_single&pnd_fields for non-scatt, greybody cloud calculations.
      //Hence, set the limits to first & last elements of the different grids.
      //Note: no margin left at lat/lon_grid edges!.
      cloudbox_limits_tmp.resize(2*atmosphere_dim);

      // Pressure limits
      cloudbox_limits_tmp[0] = 0;
      cloudbox_limits_tmp[1] = p_grid.nelem() - 1;
      // Latitude limits
      if( atmosphere_dim >= 2 )
        {
          cloudbox_limits_tmp[2] = 0;
          cloudbox_limits_tmp[3] = lat_grid.nelem() - 1;
        }
      // Longitude limits
      if( atmosphere_dim == 3 )
        {
          cloudbox_limits_tmp[4] = 0;
          cloudbox_limits_tmp[5] = lon_grid.nelem() - 1;
        }
    }
  else */
  if (cloudbox_limits.nelem() != 2 * atmosphere_dim)
    throw runtime_error(
        "*cloudbox_limits* is a vector which contains the"
        "upper and lower limit of the cloud for all "
        "atmospheric dimensions. So its dimension must"
        "be 2 x *atmosphere_dim*");
  else
    cloudbox_limits_tmp = cloudbox_limits;

  // Check that pnd_field_raw has at least 2 grid-points in each dimension.
  // Otherwise, interpolation further down will fail with assertion.

  for (Index d = 0; d < atmosphere_dim; d++) {
    for (Index i = 0; i < pnd_field_raw.nelem(); i++) {
      if (pnd_field_raw[i].get_grid_size(d) < 2) {
        ostringstream os;
        os << "Error in pnd_field_raw data. ";
        os << "Dimension " << d << " (name: \"";
        os << pnd_field_raw[i].get_grid_name(d);
        os << "\") has only ";
        os << pnd_field_raw[i].get_grid_size(d);
        os << " element";
        os << ((pnd_field_raw[i].get_grid_size(d) == 1) ? "" : "s");
        os << ". Must be at least 2.";
        throw runtime_error(os.str());
      }
    }
  }
  const Index Np_cloud = cloudbox_limits_tmp[1] - cloudbox_limits_tmp[0] + 1;

  ConstVectorView p_grid_cloud =
      p_grid[Range(cloudbox_limits_tmp[0], Np_cloud)];

  // Check that no scatterers exist outside the cloudbox
  chk_pnd_field_raw_only_in_cloudbox(atmosphere_dim,
                                     pnd_field_raw,
                                     p_grid,
                                     lat_grid,
                                     lon_grid,
                                     cloudbox_limits_tmp);

  //==========================================================================
  if (atmosphere_dim == 1) {
    ArrayOfGriddedField3 pnd_field_tmp;

    GriddedFieldPRegrid(
        pnd_field_tmp, p_grid_cloud, pnd_field_raw, 1, zeropadding, verbosity);

    FieldFromGriddedField(pnd_field,
                          p_grid_cloud,
                          pnd_field_tmp[0].get_numeric_grid(1),
                          pnd_field_tmp[0].get_numeric_grid(2),
                          pnd_field_tmp,
                          verbosity);
  } else if (atmosphere_dim == 2) {
    const Index Nlat_cloud =
        cloudbox_limits_tmp[3] - cloudbox_limits_tmp[2] + 1;

    ConstVectorView lat_grid_cloud =
        lat_grid[Range(cloudbox_limits_tmp[2], Nlat_cloud)];

    if (zeropadding) {
      // FIXME: OLE: Implement this
      CREATE_OUT0;
      out0 << "WARNING: zeropadding currently only supported for 1D.";
    }

    //Resize variables
    pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, 1);

    // Gridpositions:
    ArrayOfGridPos gp_p(Np_cloud);
    ArrayOfGridPos gp_lat(Nlat_cloud);

    // Interpolate pnd_field.
    // Loop over the scattering element:
    for (Index i = 0; i < pnd_field_raw.nelem(); ++i) {
      // Calculate grid positions:
      p2gridpos(gp_p,
                pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                p_grid_cloud);
      gridpos(gp_lat,
              pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
              lat_grid_cloud);

      // Interpolation weights:
      Tensor3 itw(Np_cloud, Nlat_cloud, 4);
      interpweights(itw, gp_p, gp_lat);

      // Interpolate:
      interp(pnd_field(i, joker, joker, 0),
             itw,
             pnd_field_raw[i].data(joker, joker, 0),
             gp_p,
             gp_lat);
    }
  } else {
    const Index Nlat_cloud =
        cloudbox_limits_tmp[3] - cloudbox_limits_tmp[2] + 1;
    const Index Nlon_cloud =
        cloudbox_limits_tmp[5] - cloudbox_limits_tmp[4] + 1;

    if (zeropadding) {
      // FIXME: OLE: Implement this
      CREATE_OUT0;
      out0 << "WARNING: zeropadding currently only supported for 1D.";
    }

    ConstVectorView lat_grid_cloud =
        lat_grid[Range(cloudbox_limits_tmp[2], Nlat_cloud)];
    ConstVectorView lon_grid_cloud =
        lon_grid[Range(cloudbox_limits_tmp[4], Nlon_cloud)];

    //Resize variables
    pnd_field.resize(pnd_field_raw.nelem(), Np_cloud, Nlat_cloud, Nlon_cloud);

    // Gridpositions:
    ArrayOfGridPos gp_p(Np_cloud);
    ArrayOfGridPos gp_lat(Nlat_cloud);
    ArrayOfGridPos gp_lon(Nlon_cloud);

    // Interpolate pnd_field.
    // Loop over the scattering element types:
    for (Index i = 0; i < pnd_field_raw.nelem(); ++i) {
      // Calculate grid positions:
      p2gridpos(gp_p,
                pnd_field_raw[i].get_numeric_grid(GFIELD3_P_GRID),
                p_grid_cloud);
      gridpos(gp_lat,
              pnd_field_raw[i].get_numeric_grid(GFIELD3_LAT_GRID),
              lat_grid_cloud);
      gridpos(gp_lon,
              pnd_field_raw[i].get_numeric_grid(GFIELD3_LON_GRID),
              lon_grid_cloud);

      // Interpolation weights:
      Tensor4 itw(Np_cloud, Nlat_cloud, Nlon_cloud, 8);
      interpweights(itw, gp_p, gp_lat, gp_lon);

      // Interpolate:
      interp(pnd_field(i, joker, joker, joker),
             itw,
             pnd_field_raw[i].data,
             gp_p,
             gp_lat,
             gp_lon);
    }
  }

  // no (cloudy) Jacobians with this WSM, hence no calc.
  // but we need to size dpnd_field to be consistent with jacobian_quantities.
  dpnd_field_dx.resize(jacobian_quantities.nelem());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldExpand1D(Tensor4& pnd_field,
                       const Index& atmosphere_dim,
                       const Index& cloudbox_on,
                       const ArrayOfIndex& cloudbox_limits,
                       const Index& nzero,
                       const Verbosity&) {
  /*  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
*/
  if (atmosphere_dim == 1) {
    throw runtime_error("No use in calling this method for 1D.");
  }
  if (!cloudbox_on) {
    throw runtime_error("No use in calling this method with cloudbox off.");
  }

  if (nzero < 1) {
    throw runtime_error("The argument *nzero* must be > 0.");
  }

  // Sizes
  const Index npart = pnd_field.nbooks();
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index nlon = 1;
  if (atmosphere_dim == 3) {
    nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  }

  if (pnd_field.npages() != np || pnd_field.nrows() != 1 ||
      pnd_field.ncols() != 1) {
    throw runtime_error(
        "The input *pnd_field* is either not 1D or does not "
        "match pressure size of cloudbox.");
  }

  // Temporary container
  Tensor4 pnd_temp = pnd_field;

  // Resize and fill
  pnd_field.resize(npart, np, nlat, nlon);
  pnd_field = 0;
  //
  for (Index ilon = nzero; ilon < nlon - nzero; ilon++) {
    for (Index ilat = nzero; ilat < nlat - nzero; ilat++) {
      for (Index ip = 0; ip < np; ip++) {
        for (Index is = 0; is < npart; is++) {
          pnd_field(is, ip, ilat, ilon) = pnd_temp(is, ip, 0, 0);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldZero(  //WS Output:
    Tensor4& pnd_field,
    ArrayOfTensor4& dpnd_field_dx,
    ArrayOfArrayOfSingleScatteringData& scat_data,
    //WS Input:
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const ArrayOfIndex& cloudbox_limits,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Verbosity&) {
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);

  if (cloudbox_limits.nelem() != 2 * atmosphere_dim)
    throw runtime_error(
        "*cloudbox_limits* is a vector which contains the"
        "upper and lower limit of the cloud for all "
        "atmospheric dimensions. So its dimension must"
        "be 2 x *atmosphere_dim*");

  // Resize pnd_field and set it to 0:
  Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index nlat = 1, nlon = 1;
  if (atmosphere_dim > 1) {
    nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    if (atmosphere_dim > 2) {
      nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    }
  }

  // no (cloudy) Jacobians with this WSM, hence no setting.
  // but we need to size dpnd_field to be consistent with jacobian_quantities.
  dpnd_field_dx.resize(jacobian_quantities.nelem());

  // Do only reset scat_data if it has not been set yet.
  // There's no need otherwise, and it's rather unpractical for testing when
  // doing so: we might want to do some actual calcs with the scat_data later
  // on. So why throw it away?
  const Index N_se = TotalNumberOfElements(scat_data);
  if (N_se > 0) {
    pnd_field.resize(N_se, np, nlat, nlon);
  } else {
    pnd_field.resize(1, np, nlat, nlon);

    //Resize scat_data and set it to 0:
    // Number of scattering elements
    scat_data.resize(1);
    scat_data[0].resize(1);
    scat_data[0][0].ptype = PTYPE_TOTAL_RND;
    scat_data[0][0].description = " ";
    // Grids which contain full ranges which one wants to calculate
    Index nf = f_grid.nelem();
    scat_data[0][0].f_grid.resize(nf);
    scat_data[0][0].f_grid = f_grid;
    Index nT = 1;
    scat_data[0][0].T_grid.resize(nT);
    scat_data[0][0].T_grid = 270.;
    Index nza = 5;
    nlinspace(scat_data[0][0].za_grid, 0, 180, nza);
    // Resize the data arrays
    scat_data[0][0].pha_mat_data.resize(nf, nT, nza, 1, 1, 1, 6);
    scat_data[0][0].pha_mat_data = 0.;
    scat_data[0][0].ext_mat_data.resize(nf, nT, 1, 1, 1);
    scat_data[0][0].ext_mat_data = 0.;
    scat_data[0][0].abs_vec_data.resize(nf, nT, 1, 1, 1);
    scat_data[0][0].abs_vec_data = 0.;
  }

  pnd_field = 0.;
}


