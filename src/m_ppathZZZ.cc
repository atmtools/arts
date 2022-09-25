/* Copyright (C) 2021 Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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
 * @file   m_ppath.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-08-11 
 *
 * @brief  Workspace functions releated to calculation of propagation paths.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "geodeticZZZ.h"
#include "math_funcs.h"
#include "ppathZZZ.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCheckStartPoint(const Ppath& ppath,
                          const Index& background,
                          const Index& np,
                          const Numeric& altitude,
                          const Numeric& daltitude,
                          const Numeric& latitude,
                          const Numeric& dlatitude,
                          const Numeric& longitude,
                          const Numeric& dlongitude,
                          const Numeric& zenith_angle,
                          const Numeric& dzenith_angle,
                          const Numeric& azimuth_angle,
                          const Numeric& dazimuth_angle,
                          const Verbosity&) {
  // pos and los to check
  ConstVectorView pos = ppath.start_pos, los = ppath.start_los;
  
  if (background >= 0 && ppath.backgroundZZZ != background) {
    ARTS_USER_ERROR("Radiative background not as expected!\n"
                    "  background in ppath: ", ppath.backgroundZZZ, "\n"
                    "  background expected: ", background );
  }
  if (np >= 0 && ppath.np != np) {
    ARTS_USER_ERROR("Number of ppath points not as expected!\n"
                    "  number in ppath: ", ppath.np, "\n"
                    "  number expected: ", np );
  }
  if (daltitude >= 0 && abs(pos[0]-altitude) > daltitude) {
    ARTS_USER_ERROR("Start altitude not as expected!\n"
                    "  altitude in ppath: ", pos[0], "\n"
                    "  altitude expected: ", altitude, "\n"
                    "         difference: ", abs(pos[0]-altitude), "\n"
                    "      set tolarance: ", daltitude );
  }
  if (dlatitude >= 0 && abs(pos[1]-latitude) > dlatitude) {
    ARTS_USER_ERROR("Start latitude not as expected!\n"
                    "  latitude in ppath: ", pos[1], "\n"
                    "  latitude expected: ", latitude, "\n"
                    "         difference: ", abs(pos[1]-latitude), "\n"
                    "      set tolarance: ", dlatitude );
  }
  if (dlongitude >= 0 && abs(pos[2]-longitude) > dlongitude) {
    ARTS_USER_ERROR("Start longitude not as expected!\n"
                    "  longitude in ppath: ", pos[2], "\n"
                    "  longitude expected: ", longitude, "\n"
                    "          difference: ", abs(pos[2]-longitude), "\n"
                    "       set tolarance: ", dlongitude );
  }
  if (dzenith_angle >= 0 && abs(los[0]-zenith_angle) > dzenith_angle) {
    ARTS_USER_ERROR("Start zenith angle not as expected!\n"
                    "  zenith angle in ppath: ", los[0], "\n"
                    "  zenith angle expected: ", zenith_angle, "\n"
                    "             difference: ", abs(los[0]-zenith_angle), "\n"
                    "          set tolarance: ", dzenith_angle );
  }
  if (dazimuth_angle >= 0 && abs(los[1]-azimuth_angle) > dazimuth_angle) {
    ARTS_USER_ERROR("Start azimuth angle not as expected!\n"
                    "  azimuth angle in ppath: ", los[1], "\n"
                    "  azimuth angle expected: ", azimuth_angle, "\n"
                    "              difference: ", abs(los[1]-azimuth_angle), "\n"
                    "           set tolarance: ", dazimuth_angle );
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathGeometric(Ppath& ppath,
                    const Vector& rte_pos,
                    const Vector& rte_los,
                    const Index& atmosphere_dim,
                    const Vector& refellipsoid,
                    const Vector& z_grid,
                    const Vector& lat_grid,
                    const Vector& lon_grid,
                    const Index& cloudbox_on,
                    const ArrayOfIndex& cloudbox_limits,
                    const GriddedField2& surface_elevation,
                    const Numeric& ppath_lmax,
                    const Numeric& ppath_stop_distance,
                    const Numeric& l_accuracy,
                    const Index& safe_surface_search,
                    const Verbosity&) {  
  // Convert rte_pos/los to ECEF
  Vector ecef(3), decef(3);
  geodetic_los2ecef(ecef, decef, rte_pos, rte_los, refellipsoid);

  // Number of points of ppath (-1 signifies still not known)
  Index np = -1;
  
  // Some initial checks and calculations
  Numeric l2toa;
  Vector pos_toa(3), los_toa(2);
  enum PpathBackground background = ppath_init_calc(l2toa,
                                                    pos_toa,
                                                    los_toa,
                                                    rte_pos,
                                                    rte_los,
                                                    ecef,
                                                    decef,
                                                    atmosphere_dim,
                                                    refellipsoid,
                                                    z_grid,
                                                    lat_grid,
                                                    lon_grid,
                                                    cloudbox_on,
                                                    cloudbox_limits);
  if (background == PPATH_BACKGROUND_SPACE)
    np = 0;
  if (background == PPATH_BACKGROUND_CLOUDBOX)
    np = 1;
  // Length from sensor to start of ppath inside the atmosphere
  const Numeric l_end = l2toa>0 ? l2toa : 0; 
  
  // If np < 0 we have a path to determine
  Numeric l_start = -1;
  if (np < 0) {
    // Length from sensor to end of ppath inside the atmosphere or
    // to tangent point if it's closer
    Numeric l_test;

    // Distance along ppath to the surface (negative if no intersection)
    // The function also checks that rte_pos is above the surface
    l_test = find_crossing_with_surface_z(rte_pos,
                                          rte_los,
                                          ecef,
                                          decef,
                                          atmosphere_dim,
                                          refellipsoid,
                                          surface_elevation,
                                          l_accuracy,
                                          safe_surface_search);
    if (l_test >= 0) {
      l_start = l_test;
      background = PPATH_BACKGROUND_SURFACE;
    }
    
    // Determine length to cloudbox (negative if no intersection)
    l_test = find_crossing_cloudbox(rte_pos,
                                    rte_los,
                                    ecef,
                                    decef,
                                    atmosphere_dim,
                                    refellipsoid,
                                    z_grid,
                                    lat_grid,
                                    lon_grid,
                                    cloudbox_on,
                                    cloudbox_limits,
                                    true);
    if (l_test >= 0 && (l_start < 0 || l_test < l_start)) {
      l_start = l_test;
      background = PPATH_BACKGROUND_CLOUDBOX;
    }

    // If still no intersection, l_start must be at TOA
    if (l_start < 0) {
      if (l2toa < 0) {
        // Observation from within the atmosphere
        l_start = intersection_altitude(ecef, decef, refellipsoid, last(z_grid));
      } else {
        // Observation from outside the atmosphere
        Vector ecef_toa(3);
        ecef_at_distance(ecef_toa, ecef, decef, l2toa);
        // Ignore lengths < 1m to find exit point, and not the entrance point
        // from which we start
        l_start = l2toa +
          intersection_altitude(ecef_toa, decef, refellipsoid, last(z_grid), 1);
      }
      background = PPATH_BACKGROUND_SPACE;
    }
    // Apply total length criterion?
    if (ppath_stop_distance > 0) {
      l_start = l2toa > 0 ? l2toa+ppath_stop_distance : ppath_stop_distance;
      background = PPATH_BACKGROUND_STOP_DISTANCE;
    }
    ARTS_ASSERT (background);

    // Determine np
    ARTS_ASSERT (l_start >= l_end);
    np = 1 + Index(ceil((l_start-l_end)/ppath_lmax));
  }

  // Fill ppath
  ppath.dim = atmosphere_dim;
  ppath.np  = np;
  ppath.backgroundZZZ = background;
  ppath.start_lstep = 0; 
  ppath.end_pos = rte_pos; 
  ppath.end_los = rte_los;
  ppath.end_lstep = l_end; 
  ppath.nreal = Vector(np, 1.0);
  ppath.ngroup = Vector(np, 1.0);
  ppath.pos.resize(np,3);
  ppath.los.resize(np,2);
  //
  if (np == 0) {
    ppath.lstep.resize(0);
  } else  if (np == 1) {
    ppath.lstep.resize(0);
    if (l2toa >= 0) {
      ppath.pos(0,joker) = pos_toa; 
      ppath.los(0,joker) = los_toa;
    } else{
      ppath.pos(0,joker) = rte_pos; 
      ppath.los(0,joker) = rte_los;
    }
  } else {
    // Create equidistant length vector and fill pos and los
    Vector l;
    nlinspace(l, l_end, l_start, np);
    for (Index i=0; i<np; i++) {
      poslos_at_distance(ppath.pos(i,joker),
                         ppath.los(i,joker),
                         ecef,
                         decef,
                         refellipsoid,
                         l[i]);
    }
    ppath.lstep.resize(np-1);
    ppath.lstep = l[1] - l[0];
  }
  if (np == 0) {
    ppath.start_pos = ppath.end_pos;
    ppath.start_los = ppath.end_los;
  } else {
    ppath.start_pos = ppath.pos(ppath.np-1,joker);
    ppath.start_los = ppath.los(ppath.np-1,joker);
  }

  // Check that ppath is fully inside the atmosphere before doing grid positions
  if (atmosphere_dim > 1 && ppath.np > 1 &&
      is_pos_outside_atmosphere(ppath.pos(ppath.np-1,joker),
                                atmosphere_dim,
                                last(z_grid),
                                lat_grid,
                                lon_grid)) {
    ConstVectorView pos = ppath.pos(ppath.np-1,joker);
    if (atmosphere_dim == 2) {
      ARTS_USER_ERROR("The propagation path is not fully inside the atmosphere.\n"
                      "The model atmosphere has this extend\n"
                      "  altitude: ", z_grid[0], "-", last(z_grid), "\n"
                      "  latitude: ", lat_grid[0], "-", last(lat_grid), "\n"
                      "while the position of *ppath* furthest away from the sensor\n"
                      "is at (h,lat): (", pos[0], ",", pos[1], ")\n"
                      "You need to extend the model atmosphere to cover this point.");
    } else {
      ARTS_USER_ERROR("The propagation path is not fully inside the atmosphere.\n"
                      "The model atmosphere has this extend\n"
                      "  altitude: ", z_grid[0], "-", last(z_grid), "\n"
                      "  latitude: ", lat_grid[0], "-", last(lat_grid), "\n"
                      " longitude: ", lon_grid[0], "-", last(lon_grid), "\n"
                      "while the position of *ppath* furthest away from the sensor\n"
                      "is at (h,lat,lon): (", pos[0], ",", pos[1], ",", pos[2], ")\n"
                      "You need to extend the model atmosphere to cover this point.");
    }
  }

  // Adjust longitudes (if 3D) and calculate grid positions
  ppath_fix_lon_and_gp(ppath, z_grid, lat_grid, lon_grid);
}
