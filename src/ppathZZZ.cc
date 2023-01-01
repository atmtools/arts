/* Copyright (C) 2021 Patrick Eriksson <patrick.eriksson@chalmers.se>

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
 * @file   ppath.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2023-01-01
 *
 * @brief  Functions releated to calculation of propagation paths.
 *
 * The term propagation path is here shortened to ppath.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include <algorithm>

#include "geodeticZZZ.h"
#include "ppathZZZ.h"
#include "variousZZZ.h"


void ppath_add_grid_crossings(Ppath& ppath,
                              const Vector& refellipsoid,
                              const Vector& z_grid,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const Numeric& l_step_max)
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
  // Excpetion: l_step_max is still a local length
  
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

      // Move to overall arrays and add points if l_step_max that requires
      for (Index i = 0; i < dl_from_ip.nelem(); ++i) {
        // Some useful lengths
        const Numeric l_next = l_acc_ppath[ip] + dl_from_ip[i];
        const Numeric dl = l_next - l_last_inserted;
        // Number of points needed to fulfill l_step_max
        if (dl > l_step_max) {
          const Index n_extra = Index(std::floor(dl / l_step_max));
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

  // The distance between last grid crossing and end point can exceed l_step_max
  // Fix (largely same code as above)!
  const Numeric dl = l_acc_ppath[ppath.np - 1] - l_last_inserted;
  if (dl > l_step_max) {
    const Index n_extra = Index(std::floor(dl / l_step_max));
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
