/**
    @file    ppath.cc
    @author  Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date    2023-01-01

    @brief   Functions releated to calculation of propagation paths (ppath).
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include <algorithm>

#include "auto_md.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "ppath.h"
#include "ppath_struct.h"
#include "surface.h"
#include "variousZZZ.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);


Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& surface_search_accuracy,
                                     const Index& surface_search_safe)
{
  // Find min and max surface altitude
  const Numeric z_min = min(surface_elevation.data);
  const Numeric z_max = max(surface_elevation.data);

  // Catch upward looking cases that can not have a surface intersection
  if (rte_pos[0] >= z_max && rte_los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (rte_pos[0] < z_max) {
    Numeric z_surf = interp_gfield2(surface_elevation, rte_pos[Range(1, 2)]);
    if (rte_pos[0] < z_surf - surface_search_accuracy)
      ARTS_USER_ERROR(
          "The sensor is below the surface. Not allowed!\n"
          "The sensor altitude is at ", rte_pos[0], " m\n"
          "The surface altitude is ", z_surf, " m\n"
          "The position is (lat,lon): (", rte_pos[1], ",", rte_pos[2], ")");
  }

  // Constant surface altitude (in comparison to *surface_search_accuracy*)
  if (z_max - z_min < surface_search_accuracy / 100) {
    // Catch cases with position on the ground, as they can fail if
    // intersection_altitude is used
    if (rte_pos[0] <= z_max) {
      return 0.0;
    } else {
      return intersection_altitude(ecef, decef, refellipsoid, z_min);
    }

    // The general case
  } else {
    // Find a distance that is guaranteed above or at surface
    // If below z_max, this distance is 0. Otherwise given by z_max
    Numeric l_min;
    if (rte_pos[0] <= z_max)
      l_min = 0;
    else {
      l_min = intersection_altitude(ecef, decef, refellipsoid, z_max);
      // No intersection if not even z_max is reached
      if (l_min < 0) return -1;
    }
    // Find max distance for search.
    // If below z_max and upward, given by z_max
    // Otherwise in general given by z_min. If z_min not reached, the distance
    // is instead given by tangent point
    Numeric l_max;
    bool l_max_could_be_above_surface = false;
    if (rte_pos[0] <= z_max && rte_los[0] <= 90) {
      l_max = intersection_altitude(ecef, decef, refellipsoid, z_max);
      l_max_could_be_above_surface = true;
    } else {
      l_max = intersection_altitude(ecef, decef, refellipsoid, z_min);
    }
    if (l_max < 0) {
      Vector ecef_tan(3);
      approx_geometrical_tangent_point(ecef_tan, ecef, decef, refellipsoid);
      l_max = ecef_distance(ecef, ecef_tan);
      // To not miss intersections just after the tangent point, we add a
      // a distance that depends om planet radius (for Earth 111 km).
      l_max += refellipsoid[0] * sin(DEG2RAD);
      l_max_could_be_above_surface = true;
    }

    // Safe but slow approach
    // ----------------------
    if (surface_search_safe) {
      Numeric l_test =
          l_min - surface_search_accuracy / 2;  // Remove l/2 to get exact result
      bool above_surface = true;   // if true l_test is 0
      while (above_surface && l_test < l_max) {
        l_test += surface_search_accuracy;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf = interp_gfield2(surface_elevation, pos[Range(1, 2)]);
        if (pos[0] < z_surf) above_surface = false;
      }
      if (above_surface) {
        return -1;
      } else {
        return l_test - surface_search_accuracy / 2;
      }

      // Bisection search
      // ----------------------
    } else {
      // If l_max matches a point above the surface, we have no intersection
      // according to this search algorithm. And the search fails. So we need
      // to check that point if status unclear
      if (l_max_could_be_above_surface) {
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_max);
        Numeric z_surf = interp_gfield2(surface_elevation, pos[Range(1, 2)]);
        if (pos[0] > z_surf) return -1;
      }
      // Start bisection
      while (l_max - l_min > 2 * surface_search_accuracy) {
        const Numeric l_test = (l_min + l_max) / 2;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf = interp_gfield2(surface_elevation, pos[Range(1, 2)]);
        if (pos[0] >= z_surf)
          l_min = l_test;
        else
          l_max = l_test;
      }
      return (l_min + l_max) / 2;
    }
  }
}


void ppath_add_grid_crossings(Ppath& ppath,
                              const Vector& refellipsoid,
                              const Vector& z_grid,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const Numeric& ppath_lstep)
{
  const Index nz = z_grid.nelem();
  const Index nlat = lat_grid.nelem();
  const Index nlon = lon_grid.nelem();

  // Nothing to do if there is no ppath step, or all grids empty
  if (ppath.np < 2  || !(nz || nlat || nlon)) {
    return;
  }

  // Extend grids to make sure they cover all points of ppath
  Vector zgrid2(nz ? nz + 2 : 0);
  if (nz) {
    zgrid2[0] = -9e6;
    zgrid2[Range(1, nz)] = z_grid;
    zgrid2[nz + 1] = 9e6;
  }
  Vector latgrid2(nlat ? nlat + 2 : 0);
  if (nlat) {
    latgrid2[0] = -91;
    latgrid2[Range(1, nlat)] = lat_grid;
    latgrid2[nlat + 1] = 91;
  }
  Vector longrid2(nlon ? nlon + 2 : 0);
  if (nlon) {
    longrid2[0] = -190;
    longrid2[Range(1, nlon)] = lon_grid;
    longrid2[nlon + 1] = 370;
  }

  // l means distance from ppath pos[0]
  // dl means distance from some other ppath point
  // Excpetion: ppath_lstep is still a local length

  // Process ppath to set up some help variables
  Vector l_acc_ppath(ppath.np);     // Accumulated length along ppath
  Vector ngp_z(nz ? ppath.np : 0);  // Grid positions as Numeric, i.e. idx+fd[0]
  Vector ngp_lat(nlat ? ppath.np : 0);
  Vector ngp_lon(nlat ? ppath.np : 0);
  //
  ArrayOfGridPos gp_z(ngp_z.nelem());
  if (nz)
    gridpos(gp_z, zgrid2, ppath.pos(joker, 0));
  ArrayOfGridPos gp_lat(ngp_lat.nelem());
  if (nlat)
    gridpos(gp_lat, latgrid2, ppath.pos(joker, 1));
  ArrayOfGridPos gp_lon(ngp_lon.nelem());
  if (nlon)
    gridpos(gp_lon, longrid2, ppath.pos(joker, 2));
  //
  for (Index ip = 0; ip < ppath.np; ++ip) {
    if (ip == 0) {
      l_acc_ppath[ip] = 0;
    } else {
      l_acc_ppath[ip] = l_acc_ppath[ip - 1] + ppath.lstep[ip - 1];
    }
    if (nz)
      ngp_z[ip] = Numeric(gp_z[ip].idx) + gp_z[ip].fd[0];
    if (nlat)
      ngp_lat[ip] = Numeric(gp_lat[ip].idx) + gp_lat[ip].fd[0];
    if (nlon)
      ngp_lon[ip] = Numeric(gp_lon[ip].idx) + gp_lon[ip].fd[0];
  }

  // Total length of ppath, minus a small distance to avoid that end point gets repeated
  const Numeric l2end = l_acc_ppath[ppath.np - 1] - 1.0e-3;

  // Containers for new ppath points (excluding start and end points, that
  // always are taken from original ppath)
  ArrayOfIndex istart_array(0);
  ArrayOfNumeric l_array(0);

  // Loop ppath steps to find grid crossings
  Numeric l_last_inserted = 0;
  for (Index ip = 0; ip < ppath.np - 1; ++ip) {

    // Length to grid crossings inside ppath step
    ArrayOfNumeric dl_from_ip(0);

    // Change in integer grid position for each dimension
    const Index dgp_z =
      nz ? n_int_between(ngp_z[ip], ngp_z[ip + 1]) : 0;
    const Index dgp_lat =
      nlat ? n_int_between(ngp_lat[ip], ngp_lat[ip + 1]) : 0;
    const Index dgp_lon =
      nlon ? n_int_between(ngp_lon[ip], ngp_lon[ip + 1]) : 0;

    if (dgp_z || dgp_lat || dgp_lon) {
      // ECEF at start end of ppath step
      Vector ecef(3), decef(3);
      geodetic_los2ecef(ecef,
                        decef,
                        ppath.pos(ip, joker),
                        ppath.los(ip, joker),
                        refellipsoid);

      // Crossing(s) of z_grid
      for (Index i = 1; i <= abs(dgp_z); ++i) {

        const Numeric dl_test = intersection_altitude(
            ecef,
            decef,
            refellipsoid,
            zgrid2[int_at_step(ngp_z[ip], sign(dgp_z) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Crossing(s) of lat_grid
      for (Index i = 1; i <= abs(dgp_lat); ++i) {
        const Numeric dl_test = intersection_latitude(
            ecef,
            decef,
            ppath.pos(ip, joker),
            ppath.los(ip, joker),
            refellipsoid,
            latgrid2[int_at_step(ngp_lat[ip], sign(dgp_lat) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Crossing(s) of lon_grid
      for (Index i = 1; i <= abs(dgp_lon); ++i) {
        const Numeric dl_test = intersection_longitude(
            ecef,
            decef,
            ppath.pos(ip, joker),
            ppath.los(ip, joker),
            longrid2[int_at_step(ngp_lon[ip], sign(dgp_lon) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Sort dl_from_ip
      std::sort(dl_from_ip.begin(), dl_from_ip.end());

      // Move to overall arrays and add points if ppath_lstep that requires
      for (Index i = 0; i < dl_from_ip.nelem(); ++i) {
        // Some useful lengths
        const Numeric l_next = l_acc_ppath[ip] + dl_from_ip[i];
        const Numeric dl = l_next - l_last_inserted;
        // Number of points needed to fulfill ppath_lstep
        if (dl > ppath_lstep) {
          const Index n_extra = Index(std::floor(dl / ppath_lstep));
          const Numeric dl_step = dl / Numeric(n_extra + 1);
          for (Index extra = 0; extra < n_extra; ++extra) {
            const Numeric l_extra = l_last_inserted + dl_step;
            istart_array.push_back(ip);
            l_array.push_back(l_extra);
            l_last_inserted = l_extra;
          }
        }
        // Add grid crossing point
        istart_array.push_back(ip);
        l_array.push_back(l_next);
        l_last_inserted = l_next;
      }
    }  // if dgp_p
  }    // ip loop

  // The distance between last grid crossing and end point can exceed ppath_lstep
  // Fix (largely same code as above)!
  const Numeric dl = l_acc_ppath[ppath.np - 1] - l_last_inserted;
  if (dl > ppath_lstep) {
    const Index n_extra = Index(std::floor(dl / ppath_lstep));
    const Numeric dl_step = dl / Numeric(n_extra + 1);
    for (Index extra = 0; extra < n_extra; ++extra) {
      const Numeric l_extra = l_last_inserted + dl_step;
      istart_array.push_back(ppath.np - 2);
      l_array.push_back(l_extra);
      l_last_inserted = l_extra;
    }
  }
  //-------------------------------------------------------------

  // Make copies of data in ppath that will change, but we need
  Index np = ppath.np;
  Vector nreal = ppath.nreal;
  Vector ngroup = ppath.ngroup;
  Matrix pos = ppath.pos;
  Matrix los = ppath.los;

  // New size of ppath
  const Index nl = l_array.nelem();

  ppath.np = nl + 2;
  ppath.nreal = Vector(ppath.np, 1.0);  // We guess on no refraction
  ppath.ngroup = Vector(ppath.np, 1.0);
  ppath.lstep.resize(ppath.np - 1);
  ppath.pos.resize(ppath.np, 3);
  ppath.los.resize(ppath.np, 2);

  // Pos and los at end points
  ppath.pos(0, joker) = pos(0, joker);
  ppath.los(0, joker) = los(0, joker);
  ppath.pos(ppath.np - 1, joker) = pos(np - 1, joker);
  ppath.los(ppath.np - 1, joker) = los(np - 1, joker);

  // Calculate and insert new pos and los, and do lstep in parallel
  Vector l_array_as_vector(nl);
  if (nl) {
    for (Index i = 0; i < nl; ++i) {
      l_array_as_vector[i] = l_array[i];
      Vector ecef(3), decef(3);
      geodetic_los2ecef(ecef,
                        decef,
                        pos(istart_array[i], joker),
                        los(istart_array[i], joker),
                        refellipsoid);
      poslos_at_distance(ppath.pos(i + 1, joker),
                         ppath.los(i + 1, joker),
                         ecef,
                         decef,
                         refellipsoid,
                         l_array[i] - l_acc_ppath[istart_array[i]]);
      if (i > 0) {
        ppath.lstep[i] = l_array[i] - l_array[i - 1];
        //ARTS_ASSERT(ppath.lstep[i] > 0)
      }
    }
    ppath.lstep[0] = l_array[0];
    ppath.lstep[nl] = l_acc_ppath[np - 1] - l_array[nl - 1];
    ARTS_ASSERT(ppath.lstep[0] > 0)
    ARTS_ASSERT(ppath.lstep[nl] > 0)
  } else {
    ppath.lstep[0] = l_acc_ppath[np - 1];
  }

  // New refractive indices, mainly set by interpolation
  if (max(nreal) > 1.0) {
    ppath.nreal[0] = nreal[0];
    ppath.ngroup[0] = ngroup[0];
    ppath.nreal[ppath.np - 1] = nreal[np - 1];
    ppath.ngroup[ppath.np - 1] = ngroup[np - 1];
    //
    if (nl) {
      ArrayOfGridPos gp(nl);
      gridpos(gp, l_acc_ppath, l_array_as_vector);
      Matrix itw(nl, 2);
      interpweights(itw, gp);
      interp(ppath.nreal[Range(1, nl)], itw, nreal, gp);
      interp(ppath.ngroup[Range(1, nl)], itw, ngroup, gp);
    }
  }
}


void ppath_extend(Ppath& ppath,
                  const Ppath& ppath2)
{
  // A crude check that last point in ppath is first in ppath2
  ARTS_ASSERT(fabs(ppath.pos(ppath.np-1, 0) - ppath2.pos(0, 0)) < 0.1);
  ARTS_ASSERT(fabs(ppath.pos(ppath.np-1, 1) - ppath2.pos(0, 1)) < 0.01);
  // We don't assert longitude, as we can have shifts of 360 deg and
  // longitude is undefined at the poles

  // Make partial copy of ppath
  Ppath ppath1;
  ppath1.np = ppath.np;
  ppath1.nreal = ppath.nreal;
  ppath1.ngroup = ppath.ngroup;
  ppath1.pos = ppath.pos;
  ppath1.los = ppath.los;
  ppath1.lstep = ppath.lstep;

  // Ranges for extended ppath
  Range part1 = Range(0, ppath1.np);
  Range part2 = Range(ppath1.np, ppath2.np);

  // Create extended ppath
  ppath.np = ppath1.np + ppath2.np;
  ppath.backgroundZZZ = ppath2.backgroundZZZ;
  // Start pos/los_lstep kept as given
  // but end pos/los_lstep should be taken from ppath2
  ppath.end_pos = ppath2.end_pos;
  ppath.end_los = ppath2.end_los;
  ppath.end_lstep = ppath2.end_lstep;
  //
  ppath.pos.resize(ppath.np, 3);
  ppath.pos(part1, joker) = ppath1.pos;
  ppath.pos(part2, joker) = ppath2.pos;
  //
  ppath.los.resize(ppath.np, 2);
  ppath.los(part1, joker) = ppath1.los;
  ppath.los(part2, joker) = ppath2.los;
  //
  ppath.nreal.resize(ppath.np);
  ppath.nreal[part1] = ppath1.nreal;
  ppath.nreal[part2] = ppath2.nreal;
  //
  ppath.ngroup.resize(ppath.np);
  ppath.ngroup[part1] = ppath1.ngroup;
  ppath.ngroup[part2] = ppath2.ngroup;
  //
  // ppath.lstep is a special case
  ppath.lstep.resize(ppath.np-1);
  ppath.lstep[Range(0,ppath1.np-1)] = ppath1.lstep;
  ppath.lstep[ppath1.np-1] = 0;
  ppath.lstep[Range(ppath1.np,ppath2.np-1)] = ppath2.lstep;
}


bool ppath_l2toa_from_above(Numeric& l2toa,
                            ConstVectorView rte_pos,
                            ConstVectorView rte_los,
                            ConstVectorView ecef,
                            ConstVectorView decef,
                            const Vector& refellipsoid,
                            const Numeric& z_toa)
{
  // Cases that are inside atmosphere
  if (rte_pos[0] < z_toa || (rte_pos[0] == z_toa && rte_los[0] > 90)) {
    l2toa = 0;
    return false;

  // Outside of the atmosphere
  } else {
    // No need to check if upward looking
    if (rte_los[0] <= 90) {
      l2toa = -1;
    } else {
      // Returns negative if no intersection
      l2toa = intersection_altitude(ecef, decef, refellipsoid, z_toa);
    }
    return true;
  }
}


void refracted_link_basic(Workspace& ws,
                          Ppath& ppath,
                          const Agenda& refr_index_air_ZZZ_agenda,
                          const Numeric& ppath_lstep,
                          const Numeric& ppath_lraytrace,
                          const Vector& refellipsoid,
                          const GriddedField2& surface_elevation,
                          const Numeric& surface_search_accuracy,
                          const Numeric& z_toa,
                          const Index& do_horizontal_gradients,
                          const Index& do_twosided_perturb,
                          const Vector& start_pos,
                          const Vector& target_pos,
                          const Numeric& target_dl,
                          const Index& max_iterations,
                          const Index& robust)
{
  // Start and target positions as ECEF
  Vector ecef_start(3), ecef_target(3);
  geodetic2ecef(ecef_start, start_pos, refellipsoid);
  geodetic2ecef(ecef_target, target_pos, refellipsoid);

  // The estimated geometric "false" target position
  // We take target_pos as first guess
  Vector target_false = target_pos;
  Vector ecef_false = ecef_target;

  // Surface elevation to use if downward looking
  GriddedField2 surface_elevation_m10km;
  LatLonFieldSetConstant(surface_elevation_m10km, -10e3, "", Verbosity());

  // Various variables used below
  Numeric l2false;
  Vector rte_los(2);
  Vector ecef(3), decef(3);  // Used for several purposes
  bool downward = false;
  Index iteration = 0;
  bool ready = false;
  bool any_failure = false;

  // Loop until ready or any failure
  while (!ready) {

    // Geometric LOS to false target
    Vector dummy(3);
    ecef_vector_distance(decef, ecef_start, ecef_false);
    decef /= norm2(decef);
    ecef2geodetic_los(dummy, rte_los, ecef_start, decef, refellipsoid);

    // Geometric distance to false target
    l2false = ecef_distance(ecef_start, ecef_false);

    // If iteration 0, check if risk for downward view. If yes we
    // iterate with surface at -10 km, and if convergence check
    // afterwards if there is an intersection
    if (!iteration && rte_los[0] > 90)
      downward = true;
    
    // Refracted with same los and ltotal
    ppathRefracted(ws,
                   ppath,
                   refr_index_air_ZZZ_agenda,
                   start_pos,
                   rte_los,
                   ppath_lstep,
                   l2false,
                   ppath_lraytrace,
                   refellipsoid,
                   downward ? surface_elevation_m10km : surface_elevation,
                   surface_search_accuracy,
                   z_toa,
                   do_horizontal_gradients,
                   do_twosided_perturb,
                   0,
                   Verbosity());

    // Intersection with surface?
    if (ppath.backgroundZZZ == PpathBackground::Surface) {
      if (robust) {
        any_failure = true;
        break;
      } else {
        ARTS_USER_ERROR("Ended up with surface intersection during iteration!\n"
                        "Target position seems to be behind the horizon\n"
                        "Giving up ....");
      }
    }
    
    // Extend ppath into space?
    if (ppath.backgroundZZZ == PpathBackground::Space) {
      const Numeric l2end = l2false - ppath.start_lstep - ppath.lstep.sum();
      if (l2end > 0) {
        geodetic_los2ecef(ecef, decef, ppath.end_pos, ppath.end_los, refellipsoid);
        ecef_at_distance(ecef, ecef, decef, l2end);
        ecef2geodetic_los(ppath.end_pos, ppath.end_los, ecef, decef, refellipsoid);
        ppath.end_lstep = l2end;
      }
    }

    // Distance to target
    geodetic2ecef(ecef, ppath.end_pos, refellipsoid);
    ecef_vector_distance(decef, ecef_target, ecef);
    Numeric dl = norm2(decef);
    if (iteration && dl < target_dl)
      ready = true;

    // Update target_geom by moving it with negative miss vector
    if (!ready) {
      decef *= -1.0;
      ecef_false += decef;

      iteration++;

      if (iteration > max_iterations) {
        if (robust) {
          any_failure = true;
          break;
        } else {
          ARTS_USER_ERROR("Maximum number of iterations reached!\n"
                          "Giving up ....");
        }
      }
    }
  }

  // Post-processing part:
  
  // If we have worked with negative surface altitude, we now make a
  // test if the path works with actual surface altitude
  if (!any_failure && downward) {
    ppathRefracted(ws,
                   ppath,
                   refr_index_air_ZZZ_agenda,
                   start_pos,
                   rte_los,
                   ppath_lstep,
                   l2false,
                   ppath_lraytrace,
                   refellipsoid,
                   surface_elevation,
                   surface_search_accuracy,
                   z_toa,
                   do_horizontal_gradients,
                   do_twosided_perturb,
                   0,
                   Verbosity());
    
    if (ppath.backgroundZZZ == PpathBackground::Surface) {
      if (robust) {
        any_failure = true;
      } else {
        ARTS_USER_ERROR("Target position can not be reached due to "
                        "surface intersection!\n");
      }
    }
  }
  
  if (!any_failure) {
    // Just remains to set background
    ppath.backgroundZZZ = PpathBackground::Transmitter;

  } else {
    // Set empty ppath
    ppath.backgroundZZZ = PpathBackground::Undefined;
    ppath.np = 0;
    ppath.pos.resize(0, 0);
    ppath.los.resize(0, 0);
    ppath.lstep.resize(0);
    ppath.nreal.resize(0);
    ppath.ngroup.resize(0);
    ppath.end_pos = ppath.start_pos;
    ppath.end_los = ppath.start_los;
    ppath.end_lstep = 0;
  }
}


void specular_los_calc(VectorView los_new,
                       const Vector& refellipsoid,
                       const GriddedField2& surface_elevation,
                       ConstVectorView pos2D,
                       ConstVectorView los,
                       const bool& ignore_topography)
{
  ARTS_ASSERT(los_new.nelem() == 2);
  ARTS_ASSERT(pos2D.nelem() == 2);
  ARTS_ASSERT(los.nelem() == 2);

  // Ignore surface tilt if told so or surface_elevation.data has size (1,1)
  if (ignore_topography || (surface_elevation.data.nrows() == 1 &&
                            surface_elevation.data.ncols() == 1)) {
    los_new[0] = 180 - los[0];
    los_new[1] = los[1];

  } else {
    // Determine surface normal
    Vector pos(3), ecef(3), decef(3);
    surface_normal(pos, ecef, decef, refellipsoid, surface_elevation, pos2D);

    // Convert los to ECEF direction (ECEF recalculated!)
    Vector decef_los(3);
    geodetic_los2ecef(ecef, decef_los, pos, los, refellipsoid);

    // Dot product between normal and los should be negative.
    ARTS_USER_ERROR_IF(decef * decef_los > 0,
        "The incoming direction is from below the surface!");

    // If v is incoming los, n is the normal and the reflected
    // direction, w; is: w = v - 2 * (v * n) * n
    // See e.g.
    // http://www.sunshine2k.de/articles/coding/vectorreflection/
    // vectorreflection.html
    Vector change = decef;
    change *= -2.0 * (decef_los * decef);
    Vector decef_new_los = decef_los;
    decef_new_los += change;

    // Convert to LOS (pos recalculated and we use this for a rough assert)
    ecef2geodetic_los(pos, los_new, ecef, decef_new_los, refellipsoid);
    ARTS_ASSERT(fabs(pos[1] - pos2D[0]) < 0.01);
    // We don't assert longitude to avoid problems at lat=90 and
    // possible shifts of 360 deg
  }
}
