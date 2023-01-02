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

#include "check_input.h"
#include "geodeticZZZ.h"
#include "ppathZZZ.h"
#include "variousZZZ.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathAddGridCrossings(Ppath& ppath,
                           const Vector& refellipsoid,
                           const Numeric& l_step_max,
                           const Vector& z_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Verbosity&)
{
  chk_if_positive("l_step_max", l_step_max);
  
  ppath_add_grid_crossings(ppath,
                           refellipsoid,
                           z_grid,
                           lat_grid,
                           lon_grid,
                           l_step_max);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCheckEndPoint(const Ppath& ppath,
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
                        const Verbosity&)
{
  // pos and los to check
  ConstVectorView pos = ppath.end_pos, los = ppath.end_los;

  ARTS_USER_ERROR_IF(background >= 0 && ppath.backgroundZZZ != background,
                     "Radiative background not as expected!\n"
                     "  background in ppath: ", ppath.backgroundZZZ,
                     "\n  background expected: ", background);
 
  ARTS_USER_ERROR_IF(np >= 0 && ppath.np != np,
                     "Number of ppath points not as expected!\n"
                     "  number in ppath: ", ppath.np,
                     "\n  number expected: ", np);

  ARTS_USER_ERROR_IF(daltitude >= 0 && abs(pos[0] - altitude) > daltitude,
                     "End altitude not as expected!\n"
                     "  altitude in ppath: ", pos[0],
                     "\n  altitude expected: ", altitude,
                     "\n         difference: ", abs(pos[0] - altitude),
                     "\n      set tolarance: ", daltitude);
  ARTS_USER_ERROR_IF(dlatitude >= 0 && abs(pos[1] - latitude) > dlatitude,
                     "End latitude not as expected!\n"
                     "  latitude in ppath: ", pos[1],
                     "\n  latitude expected: ", latitude,
                     "\n         difference: ", abs(pos[1] - latitude),
                     "\n      set tolarance: ", dlatitude);
  
  ARTS_USER_ERROR_IF(dlongitude >= 0 && abs(pos[2] - longitude) > dlongitude,
                     "End longitude not as expected!\n"
                     "  longitude in ppath: ", pos[2],
                     "\n  longitude expected: ", longitude,
                     "\n          difference: ", abs(pos[2] - longitude),
                     "\n       set tolarance: ", dlongitude);

  ARTS_USER_ERROR_IF(dzenith_angle >= 0 && abs(los[0] - zenith_angle) > dzenith_angle,
                     "End zenith angle not as expected!\n"
                     "  zenith angle in ppath: ", los[0],
                     "\n  zenith angle expected: ", zenith_angle,
                     "\n             difference: ", abs(los[0] - zenith_angle),
                     "\n          set tolarance: ", dzenith_angle);
  ARTS_USER_ERROR_IF(dazimuth_angle >= 0 && abs(los[1] - azimuth_angle) > dazimuth_angle,
                     "End azimuth angle not as expected!\n"
                     "  azimuth angle in ppath: ", los[1],
                     "\n  azimuth angle expected: ", azimuth_angle,
                     "\n              difference: ", abs(los[1] - azimuth_angle),
                     "\n           set tolarance: ", dazimuth_angle);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathGeometric(Ppath& ppath,
                    const Vector& rte_pos,
                    const Vector& rte_los,
                    const Vector& refellipsoid,
                    const GriddedField2& surface_elevation,
                    const Numeric& z_toa,
                    const Numeric& l_step_max,
                    const Numeric& l_total_max,
                    const Numeric& surface_search_accuracy,
                    const Index& surface_search_safe,
                    const Verbosity&)
{
  chk_rte_pos("rte_pos", rte_pos);
  chk_rte_los("rte_los", rte_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_surface_elevation(surface_elevation);
  chk_if_positive("z_toa", z_toa);
  chk_if_positive("l_step_max", l_step_max);
  chk_if_positive("surface_search_accuracy", surface_search_accuracy);
  chk_if_bool("surface_search_safe", surface_search_safe);

  // Convert rte_pos/los to ECEF
  Vector ecef(3), decef(3);
  geodetic_los2ecef(ecef, decef, rte_pos, rte_los, refellipsoid);

  // Relate rte_pos to TOA
  Numeric l_outside, l_inside = -1;
  const bool start_in_space = ppath_l2toa_from_above(l_outside,
                                                     rte_pos,
                                                     rte_los,
                                                     ecef,
                                                     decef,
                                                     refellipsoid,
                                                     z_toa);

  // Number of points in ppath and radiative background
  Index np = -1;  // -1 flags not yet known
  enum PpathBackground background = PPATH_BACKGROUND_UNDEFINED;
  
  // No ppath if above and looking outside of atmosphere
  if (start_in_space && l_outside < 0) {
    np = 0;
    background = PPATH_BACKGROUND_SPACE;

  // We have a path! 
  } else {
    // Distance to the surface (negative if no intersection)
    // This is from TOA if sensor outside
    // The function also checks that rte_pos is above the surface
    l_inside = find_crossing_with_surface_z(rte_pos,
                                            rte_los,
                                            ecef,
                                            decef,
                                            refellipsoid,
                                            surface_elevation,
                                            surface_search_accuracy,
                                            surface_search_safe);
    l_inside -= l_outside;

    // If intersection with surface, we have found end
    if (l_inside > 0) {
      background = PPATH_BACKGROUND_SURFACE;
    // If not, end must be TOA, but
    } else {
      if (start_in_space) {
        // We have a limb sounding from space.
        // We need to calculate from where ppath enters the atmosphere
        Vector ecef_toa(3);
        ecef_at_distance(ecef_toa, ecef, decef, l_outside);
        // Ignore lengths < 1m to find exit point, and not the entrance point
        // from which we start
        l_inside = intersection_altitude(ecef_toa, decef, refellipsoid, z_toa, 1.0);
      } else {
        // We have upward or limb, both from within the atmosphere
        l_inside = intersection_altitude(ecef, decef, refellipsoid, z_toa);
      }
      background = PPATH_BACKGROUND_SPACE;
    }

    // Consider l_total_max
    if (l_total_max > 0 && l_inside > l_total_max) {
      l_inside = l_total_max;
      background = PPATH_BACKGROUND_STOP_DISTANCE;
    }

    // Determine np and l_step
    ARTS_ASSERT(l_inside > 0);
    np = 1 + Index(ceil(l_inside / l_step_max));
  }

  // Fill ppath
  ppath.np = np;
  ARTS_ASSERT(background != PPATH_BACKGROUND_UNDEFINED);
  ppath.backgroundZZZ = background;
  ppath.start_lstep = 0;
  ppath.start_pos = rte_pos;
  ppath.start_los = rte_los;
  ppath.start_lstep = l_outside > 0 ? l_outside : 0;
  ppath.end_lstep = 0.0;
  ppath.nreal = Vector(np, 1.0);
  ppath.ngroup = Vector(np, 1.0);
  ppath.pos.resize(np, 3);
  ppath.los.resize(np, 2);
  //
  if (np == 0) {
    ppath.lstep.resize(0);
  } else {
    // Create an equidistant length vector and fill pos and los
    Vector l;
    nlinspace(l, l_outside, l_outside + l_inside, np);
    for (Index i = 0; i < np; i++) {
      poslos_at_distance(ppath.pos(i, joker),
                         ppath.los(i, joker),
                         ecef,
                         decef,
                         refellipsoid,
                         l[i]);
    }
    ppath.lstep.resize(np - 1);
    ppath.lstep = l[1] - l[0];
  }
  if (np == 0) {
    ppath.end_pos = ppath.start_pos;
    ppath.end_los = ppath.start_los;
  } else {
    ppath.end_pos = ppath.pos(ppath.np - 1, joker);
    ppath.end_los = ppath.los(ppath.np - 1, joker);
  }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void ppathRefracted(Ppath& ppath,
                    const Vector& rte_pos,
                    const Vector& rte_los,
                    const Vector& refellipsoid,
                    const GriddedField2& surface_elevation,
                    const Numeric& z_toa,
                    const Numeric& l_step_max,
                    const Numeric& l_total_max,
                    const Numeric& l_raytrace,
                    const Verbosity&)
{
  chk_rte_pos("rte_pos", rte_pos);
  chk_rte_los("rte_los", rte_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_surface_elevation(surface_elevation);
  chk_if_positive("z_toa", z_toa);
  chk_if_positive("l_step_max", l_step_max);

  // Convert rte_pos/los to ECEF
  Vector ecef(3), decef(3);
  geodetic_los2ecef(ecef, decef, rte_pos, rte_los, refellipsoid);

  // Relate rte_pos to TOA
  Numeric l_outside;
  const bool start_in_space = ppath_l2toa_from_above(l_outside,
                                                     rte_pos,
                                                     rte_los,
                                                     ecef,
                                                     decef,
                                                     refellipsoid,
                                                     z_toa);

  // Number of points in ppath and radiative background
  Index np = -1;  // -1 flags not yet known
  enum PpathBackground background = PPATH_BACKGROUND_UNDEFINED;

  // Containers for found pos, los and lstep
  Array<Vector> pos_a;
  Array<Vector> los_a;
  Array<Numeric> lstep_a;

  // No ppath if above and looking outside of atmosphere
  if (start_in_space && l_outside < 0) {
    np = 0;
    background = PPATH_BACKGROUND_SPACE;

  // We have a path! 
  } else {
    // Variables representing latest known position of ppath 
    Vector pos0(3), los0(2), ecef0(3);

    // Init these variables
    if (start_in_space) {
        // We need to calculate from where ppath enters the atmosphere
        ecef_at_distance(ecef0, ecef, decef, l_outside);
        ecef2geodetic_los(pos0, los0, ecef0, decef, refellipsoid);
    } else {
      pos0 = rte_pos;
      los0 = rte_los;
      ecef0 = ecef;
    }

    // Append pos0 and los0 to array variables (no lstep to add yet)
    np = 1;
    pos_a.push_back(pos0);
    los_a.push_back(los0);

    // Actual ray tracing length
    const Index n_rt_per_step =
      l_raytrace > 0 ? Index(ceil(l_step_max / l_raytrace)) : 1;
    const Numeric l_rt = l_step_max / Numeric(n_rt_per_step);
    
    // Variables for next ray tracing step
    Vector pos_try(3), los_try(2), ecef_try(3);
    Numeric l2pos0;  // Actual length of step

    // Help length variables
    Numeric l_from_start = 0.0;
    Numeric l_this_step = 0.0;
    Index n_this_step = 0;
    
    // Loop as long we are inside the atmosphere
    //
    bool inside = true;
    //
    while (inside) {
      
      // Move forward with l_rt
      ecef_at_distance(ecef_try, ecef0, decef, l_rt);
      l2pos0 = l_rt;  // Can be changed below
      ecef2geodetic_los(pos_try, los_try, ecef_try, decef, refellipsoid);

      // Check if we still are inside. If not, determine end point
      // Above TOA?
      if (pos_try[0] >= z_toa) {
        inside = false;
        background = PPATH_BACKGROUND_SPACE;
        l2pos0 = intersection_altitude(ecef0, decef, refellipsoid, z_toa);
        ecef_at_distance(ecef0, ecef0, decef, l2pos0);
        ecef2geodetic_los(pos0, los0, ecef0, decef, refellipsoid);
      }

      // Passed active l_total_max?
      else if (l_total_max > 0 && l_total_max <= l_from_start + l2pos0) {
        // Fill and extend if condition
        inside = false;
        background = PPATH_BACKGROUND_STOP_DISTANCE;
        l2pos0 = l_total_max - l_from_start;
        ecef_at_distance(ecef0, ecef0, decef, l2pos0);
        ecef2geodetic_los(pos0, los0, ecef0, decef, refellipsoid);

        // Below surface?
      } else {
        const Numeric z_surface = surface_z_at_pos(pos_try, surface_elevation);
        if (pos_try[0] <= z_surface) {
          inside = false;
          background = PPATH_BACKGROUND_SURFACE;
          l2pos0 = find_crossing_with_surface_z(pos0,
                                                los0,
                                                ecef0,
                                                decef,
                                                refellipsoid,
                                                surface_elevation,
                                                l_raytrace/100,
                                                0);
          ecef_at_distance(ecef0, ecef0, decef, l2pos0);
          ecef2geodetic_los(pos0, los0, ecef0, decef, refellipsoid);
        }
      }  

      // If inside, then we take try values
      if (inside) {
        pos0 = pos_try;
        los0 = los_try;
        ecef0 = ecef_try;
        l2pos0 = l_rt;
      }

      // Update step variables
      l_this_step += l2pos0;
      n_this_step++;

      // Add to arrays if ready with either full path or step
      if (!inside || n_this_step == n_rt_per_step) {
        pos_a.push_back(pos0);
        los_a.push_back(los0);
        lstep_a.push_back(l_this_step);
        ++np;
        l_this_step = 0.0;
        n_this_step = 0;
      }
    }
  }

  // Fill ppath
  ppath.np = np;
  ARTS_ASSERT(background != PPATH_BACKGROUND_UNDEFINED);
  ppath.backgroundZZZ = background;
  ppath.start_lstep = 0;
  ppath.start_pos = rte_pos;
  ppath.start_los = rte_los;
  ppath.start_lstep = l_outside > 0 ? l_outside : 0;
  ppath.end_lstep = 0.0;
  ppath.nreal = Vector(np, 1.0);  // !!!
  ppath.ngroup = Vector(np, 1.0); // !!!
  ppath.pos.resize(np, 3);
  ppath.los.resize(np, 2);
  //
  if (np == 0) {
    ppath.lstep.resize(0);
  } else {
    ppath.lstep.resize(np - 1);
    for (Index i=0; i<np; ++i) {
      ppath.pos(i, joker) = pos_a[i];
      ppath.los(i, joker) = los_a[i];
      if (i < np - 1)
        ppath.lstep[i] = lstep_a[i];
    }
  }
  if (np == 0) {
    ppath.end_pos = ppath.start_pos;
    ppath.end_los = ppath.start_los;
  } else {
    ppath.end_pos = ppath.pos(ppath.np - 1, joker);
    ppath.end_los = ppath.los(ppath.np - 1, joker);
  }
}

