/**
    @file    m_ppath.cc
    @author  Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date    2023-01-14
 
    @brief   Workspace methods related to complete propagation paths
             (ppath), including to extract *geo_pos*.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "ppath.h"
#include "ppath_struct.h"
#include "variousZZZ.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posEndOfPpath(Vector& geo_pos,
                       const Ppath& ppath,
                       const Verbosity& verbosity) {
  geo_pos.resize(5);

  if (ppath.np) {
    geo_pos[Range(0, 3)] = ppath.pos(ppath.np - 1, joker);
    geo_pos[Range(3, 2)] = ppath.los(ppath.np - 1, joker);

  } else {
    geo_pos = NAN;
  }
  
  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posLowestAltitudeOfPpath(Vector& geo_pos,
                                  const Ppath& ppath,
                                  const Verbosity& verbosity) {
  geo_pos.resize(5);

  if (ppath.np) {
    // Take first point of ppath as first guess
    geo_pos[Range(0, 3)] = ppath.pos(0, joker);
    geo_pos[Range(3, 2)] = ppath.los(0, joker);

    for (Index i = 1; i < ppath.np; i++) {
      if (ppath.pos(i, 0) < geo_pos[0]) {
        geo_pos[Range(0, 3)] = ppath.pos(i, joker);
        geo_pos[Range(3, 2)] = ppath.los(i, joker);
      }
    }

  } else {
    geo_pos = NAN;
  }

  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posWhereAltitudeIsPassed(Vector& geo_pos,
                                  const Ppath& ppath,
                                  const Numeric& altitude,
                                  const Verbosity& verbosity) {
  geo_pos.resize(5);
  geo_pos = NAN;

  if (ppath.np) {
    bool found = false;
    Index ihit = 0;
    bool above = false;

    if (ppath.pos(0, 0) >= altitude) {
      above = true;
    }

    while (!found && ihit < ppath.np - 1) {
      ihit += 1;
      if (above && ppath.pos(ihit, 0) < altitude) {
        found = true;
      } else if (!above && ppath.pos(ihit, 0) >= altitude) {
        found = true;
      }
    }
    
    if (found) {
      geo_pos[0] = altitude;
      
      // Make a simple linear interpolation to determine lat, lon, za and aa
      const Numeric w = (altitude - ppath.pos(ihit - 1, 0)) /
        (ppath.pos(ihit, 0) - ppath.pos(ihit - 1, 0));
      
      geo_pos[1] = w * ppath.pos(ihit, 1) + (1 - w) * ppath.pos(ihit - 1, 1);
      geo_pos[2] = w * ppath.pos(ihit, 2) + (1 - w) * ppath.pos(ihit - 1, 2);
      geo_pos[3] = w * ppath.los(ihit, 0) + (1 - w) * ppath.los(ihit - 1, 0);
      geo_pos[4] = w * ppath.los(ihit, 1) + (1 - w) * ppath.los(ihit - 1, 1);
    }
  }
  
  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathAddGridCrossings(Ppath& ppath,
                           const Numeric& ppath_lstep,
                           const Vector& refellipsoid,
                           const Vector& z_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Verbosity&)
{
  chk_if_positive("ppath_lstep", ppath_lstep);

  ppath_add_grid_crossings(ppath,
                           refellipsoid,
                           z_grid,
                           lat_grid,
                           lon_grid,
                           ppath_lstep);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCheckEndPoint(const Ppath& ppath,
                        const String& background,
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
  const PpathBackground ppath_background = Options::toPpathBackgroundOrThrow(background);

  ARTS_USER_ERROR_IF(ppath_background != PpathBackground::Undefined &&
                     ppath.backgroundZZZ != ppath_background,
      "Radiative background not as expected!\n"
      "  background in ppath: ",
      ppath.backgroundZZZ, "\n  background expected: ", ppath_background);

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
void ppathCheckInsideDomain(const Ppath& ppath,
                            const Numeric& lat_min,
                            const Numeric& lat_max,
                            const Numeric& lon_min,
                            const Numeric& lon_max,
                            const Verbosity&)
{
  // Start from end, as end point likely most outside
  for (Index i=ppath.np - 1; i>=0; --i)
    ARTS_USER_ERROR_IF(ppath.pos(i, 1) < lat_min ||
                       ppath.pos(i, 1) > lat_max ||
                       ppath.pos(i, 2) < lon_min ||
                       ppath.pos(i, 2) > lon_max,
        "At least one point of propagation path outside of specified domain!\n",
        "     latitude limits of domain: [", lat_min, ", ", lat_max, "]\n"
        "    longitude limits of domain: [", lon_min, ", ", lon_max, "]\n"
        "  point found at (z, lat, lon): (", ppath.pos(i, 0), ", ",
                       ppath.pos(i, 1), ", ", ppath.pos(i, 2), ")\n");
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCheckInsideGrids(const Ppath& ppath,
                           const Vector& latitude_grid,
                           const Vector& longitude_grid,
                           const Verbosity& verbosity)
{
  ppathCheckInsideDomain(ppath,
                         latitude_grid[0],
                         last(latitude_grid),
                         longitude_grid[0],
                         last(longitude_grid),
                         verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathGeometric(Ppath& ppath,
                    const Vector& rte_pos,
                    const Vector& rte_los,
                    const Numeric& ppath_lstep,
                    const Numeric& ppath_ltotal,
                    const Vector& refellipsoid,
                    const GriddedField2& surface_elevation,
                    const Numeric& surface_search_accuracy,
                    const Index& surface_search_safe,
                    const Numeric& z_toa,
                    const Index& include_specular_ppath,
                    const Verbosity& verbosity)
{
  chk_rte_pos("rte_pos", rte_pos);
  chk_rte_los("rte_los", rte_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_surface_elevation(surface_elevation);
  chk_if_positive("z_toa", z_toa);
  chk_if_positive("ppath_lstep", ppath_lstep);
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
  PpathBackground background = PpathBackground::Undefined;

  // No ppath if above and looking outside of atmosphere
  if (start_in_space && l_outside < 0) {
    np = 0;
    background = PpathBackground::Space;

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
      background = PpathBackground::Surface;
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
      background = PpathBackground::Space;
    }

    // Consider ppath_ltotal
    if (ppath_ltotal > 0 && l_inside > ppath_ltotal) {
      l_inside = ppath_ltotal;
      background = PpathBackground::StopDistance;
    }

    // Determine np and l_step
    ARTS_ASSERT(l_inside > 0);
    np = 1 + Index(ceil(l_inside / ppath_lstep));
  }

  // Fill ppath
  ppath.np = np;
  ARTS_ASSERT(background != PpathBackground::Undefined);
  ppath.backgroundZZZ = background;
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

  // If surface intersection, include part beyond?
  if (include_specular_ppath && ppath.backgroundZZZ == PpathBackground::Surface) {

    Vector pos{ppath.pos(ppath.np-1, joker)};
    Vector los(2);
    specular_los_calc(los,
                      refellipsoid,
                      surface_elevation,
                      pos[Range(1, 2)],
                      ppath.los(ppath.np-1, joker));
    
    Ppath ppath2;
    ppathGeometric(ppath2,
                   pos,
                   los,
                   ppath_lstep,
                   ppath_ltotal - sum(ppath.lstep),
                   refellipsoid,
                   surface_elevation,
                   surface_search_accuracy,
                   surface_search_safe,
                   z_toa,
                   include_specular_ppath,
                   verbosity);

    ppath_extend(ppath, ppath2);
  }  
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathRefracted(Workspace& ws,
                    Ppath& ppath,
                    const Agenda& refr_index_air_ZZZ_agenda,
                    const Vector& rte_pos,
                    const Vector& rte_los,
                    const Numeric& ppath_lstep,
                    const Numeric& ppath_ltotal,
                    const Numeric& ppath_lraytrace,
                    const Vector& refellipsoid,
                    const GriddedField2& surface_elevation,
                    const Numeric& surface_search_accuracy,
                    const Numeric& z_toa,
                    const Index& do_horizontal_gradients,
                    const Index& do_twosided_perturb,
                    const Index& include_specular_ppath,
                    const Verbosity& verbosity)
{
  chk_rte_pos("rte_pos", rte_pos);
  chk_rte_los("rte_los", rte_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_surface_elevation(surface_elevation);
  chk_if_positive("z_toa", z_toa);
  chk_if_positive("ppath_lstep", ppath_lstep);

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
  PpathBackground background = PpathBackground::Undefined;

  // Containers for found pos, los, lstep and refractive indices
  Array<Vector> pos_a;
  Array<Vector> los_a;
  Array<Numeric> lstep_a, nreal_a, ngroup_a;

  // No ppath if above and looking outside of atmosphere
  if (start_in_space && l_outside < 0) {
    np = 0;
    background = PpathBackground::Space;

  // We have a path!
  } else {
    // Variables representing latest known position of ppath
    Vector pos0(3), los0(2), ecef0(3), decef0(3);

    // Init these variables
    if (start_in_space) {
        // We need to calculate from where ppath enters the atmosphere
        ecef_at_distance(ecef0, ecef, decef, l_outside);
        decef0 = decef;
        ecef2geodetic_los(pos0, los0, ecef0, decef0, refellipsoid);
    } else {
      pos0 = rte_pos;
      los0 = rte_los;
      ecef0 = ecef;
      decef0 = decef;
    }

    // Refractive index and its gradients at pos0
    Numeric n_real, n_group, dndz, dndlat, dndlon;
    refr_index_and_its_gradients(n_real,
                                 n_group,
                                 dndz,
                                 dndlat,
                                 dndlon,
                                 ws,
                                 refr_index_air_ZZZ_agenda,
                                 pos0,
                                 do_horizontal_gradients,
                                 do_twosided_perturb);
    
    // Append pos0 and los0 to array variables (no lstep to add yet)
    np = 1;
    pos_a.push_back(pos0);
    los_a.push_back(los0);
    nreal_a.push_back(n_real);
    ngroup_a.push_back(n_group);
    
    // Actual ray tracing length
    const Index n_rt_per_step =
      ppath_lraytrace > 0 ? Index(ceil(ppath_lstep / ppath_lraytrace)) : 1;
    const Numeric l_rt = ppath_lstep / Numeric(n_rt_per_step);

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
      ecef_at_distance(ecef_try, ecef0, decef0, l_rt);
      l2pos0 = l_rt;  // Can be changed below
      ecef2geodetic_los(pos_try, los_try, ecef_try, decef0, refellipsoid);
      
      // Check if we still are inside. If not, determine end point
      // Above TOA?
      if (pos_try[0] >= z_toa) {
        inside = false;
        background = PpathBackground::Space;
        l2pos0 = intersection_altitude(ecef0, decef0, refellipsoid, z_toa);
        ecef_at_distance(ecef0, ecef0, decef0, l2pos0);
        ecef2geodetic_los(pos0, los0, ecef0, decef0, refellipsoid);
      }

      // Passed active ppath_ltotal?
      else if (ppath_ltotal > 0 && ppath_ltotal <= l_from_start + l2pos0) {
        // Fill and extend if condition
        inside = false;
        background = PpathBackground::StopDistance;
        l2pos0 = ppath_ltotal - l_from_start;
        ecef_at_distance(ecef0, ecef0, decef0, l2pos0);
        ecef2geodetic_los(pos0, los0, ecef0, decef0, refellipsoid);

        // Below surface?
      } else {
        const Numeric z_surface = interp_gfield2(surface_elevation,
                                                 Vector{pos_try[Range(1, 2)]});
        if (pos_try[0] <= z_surface) {
          inside = false;
          background = PpathBackground::Surface;
          l2pos0 = find_crossing_with_surface_z(pos0,
                                                los0,
                                                ecef0,
                                                decef0,
                                                refellipsoid,
                                                surface_elevation,
                                                surface_search_accuracy,
                                                0);
          ecef_at_distance(ecef0, ecef0, decef0, l2pos0);
          ecef2geodetic_los(pos0, los0, ecef0, decef0, refellipsoid);
        }
      }

      // If inside, then we take try values and apply refraction
      if (inside) {
        pos0 = pos_try;
        los0 = los_try;
        ecef0 = ecef_try;
        l2pos0 = l_rt;
        //
        // Calculate n and its gradients for new pos0
        refr_index_and_its_gradients(n_real,
                                     n_group,
                                     dndz,
                                     dndlat,
                                     dndlon,
                                     ws,
                                     refr_index_air_ZZZ_agenda,
                                     pos0,
                                     do_horizontal_gradients,
                                     do_twosided_perturb);
        //
        // Update LOS angles
        //
        // Tried to use mean n and dndz over the step, but got basically
        // identical result (as dndz is constant in each layer)
        //
        // The expressions below are found in ARTS theory document, but
        // are also commented below this function
        //
        if (do_horizontal_gradients) {
          const Numeric sinza = sin(DEG2RAD * los0[0]);
          const Numeric cosza = cos(DEG2RAD * los0[0]);
          const Numeric sinaa = sin(DEG2RAD * los0[1]);
          const Numeric cosaa = cos(DEG2RAD * los0[1]);
          const Numeric r = norm2( ecef0 );
          const Numeric dndlatp =  dndlat / r;
          // Make sure that we don't divide with zero (if lat = +-90)
          const Numeric dndlonp =  dndlon / (r * max(cos(DEG2RAD * pos0[1]), 1e-6));
          const Numeric fac = (RAD2DEG * l_rt / n_real);
          //
          los0[0] += fac * (-sinza * dndz +
                            cosza * (cosaa * dndlatp + sinaa * dndlonp));
          los0[1] += fac * sinza * (-sinaa * dndlatp + cosaa * dndlonp);
          // Make sure we are inside [0,180] and [-180,180]
          if (los0[0] < 0) {
            los0[0] = -los0[0];
            los0[1] += 180;
          } else if (los0[0] > 180) {
            los0[0] = 360 - los0[0];
            los0[1] += 180;
          }
          if (los0[1] < -180)
            los0[1] += 360;
          else if (los0[1] > 180)
            los0[1] -= 360;

        // Just vertical gradient to consider:
        } else {
          const Numeric sinza = sin(DEG2RAD * los0[0]);
          //
          los0[0] -= (RAD2DEG * l_rt / n_real) * (sinza * dndz);
        }
        // Don't forget to update decef0! (efecf0 recalculated)
        geodetic_los2ecef(ecef0, decef0, pos0, los0, refellipsoid);
        
      // If not inside, we just need to determine refractive index
      } else {
        refr_index_air_ZZZ_agendaExecute(ws,
                                         n_real,
                                         n_group,
                                         pos0,
                                         refr_index_air_ZZZ_agenda);

      }

      // Update step variables
      l_from_start += l2pos0;
      l_this_step += l2pos0;
      n_this_step++;

      // Add to arrays if ready with either full path or step
      if (!inside || n_this_step == n_rt_per_step) {
        //
        pos_a.push_back(pos0);
        los_a.push_back(los0);
        lstep_a.push_back(l_this_step);
        nreal_a.push_back(n_real);
        ngroup_a.push_back(n_group);
        //
        ++np;
        l_this_step = 0.0;
        n_this_step = 0;
      }
    }
  }

  // Fill ppath
  ppath.np = np;
  ARTS_ASSERT(background != PpathBackground::Undefined);
  ppath.backgroundZZZ = background;
  ppath.start_pos = rte_pos;
  ppath.start_los = rte_los;
  ppath.start_lstep = l_outside > 0 ? l_outside : 0;
  ppath.end_lstep = 0.0;
  ppath.nreal = Vector(np);  
  ppath.ngroup = Vector(np);
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
      ppath.nreal[i] = nreal_a[i];
      ppath.ngroup[i] = ngroup_a[i];
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

  // If surface intersection, include part beyond?
  if (include_specular_ppath && ppath.backgroundZZZ == PpathBackground::Surface) {

    Vector pos{ppath.pos(ppath.np-1, joker)};
    Vector los(2);
    specular_los_calc(los,
                      refellipsoid,
                      surface_elevation,
                      pos[Range(1, 2)],
                      ppath.los(ppath.np-1, joker));
    
    Ppath ppath2;
    ppathRefracted(ws,
                   ppath2,
                   refr_index_air_ZZZ_agenda,
                   pos,
                   los,
                   ppath_lstep,
                   ppath_ltotal - sum(ppath.lstep),
                   ppath_lraytrace,
                   refellipsoid,
                   surface_elevation,
                   surface_search_accuracy,
                   z_toa,
                   do_horizontal_gradients,
                   do_twosided_perturb,
                   include_specular_ppath,
                   verbosity);

    ppath_extend(ppath, ppath2);
  }
}
// Comments on expressions for effect of refraction:
//
// The expressions used are an extension of Eq 9.33 in Rodgers book
// "Inverse methods for atmospheric sounding". That equation deals
// with dza/dl for a 2D geometry. Here we assume that all angles are
// in radians and we don't care about DEG2RAD and RAD2DEG terms.
// l, r, za and aa are used for length radius, zenith angle and
// azimuth angle, respectively.
//
// Without horizontal gradients, the expressiuon applied comes
// directly from Rodgers equation.
//
// A first consideration for horizontal gradients is to convert them
// from change per angle to change per meter. If we denote the
// converted gradients with prime, these conversions are:
//   dn/dlat' = dn/dlat / r 
//   dn/dlon' = dn/dlon / (r * cos(lat))
//
// The rest is just about applying angles correctly (note that n is
// placed on left side to keep expressions more clean):
//
// n * dza/dl = -sin(za) * dn/dz + 
//               cos(za) * cos(aa) * dn/dlat' +
//               cos(za) * sin(aa) * dn/dlon'
// 
// n * daa/dl = -sin(za) * sin(aa) * dn/dlat' +
//               sin(za) * cos(aa) * dn/dlon'
//
// Hopefully I have got it right! The expression for dza/dl agrees
// with Rodgers equation for aa=0. And some tests made and deviations
// from geomtrical path had at least the correct sign.
//
// Patrick 230106
