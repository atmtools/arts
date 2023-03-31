/* Copyright (C) 2002-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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
 * All in this file should be revised, to either be updated or deleted.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic_OLD.h"
#include "lin_alg.h"
#include "m_general.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath_OLD.h"
#include "refraction.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCalc(Workspace& ws,
               Ppath& ppath,
               const Agenda& ppath_agenda,
               const Numeric& ppath_lmax,
               const Numeric& ppath_lraytrace,
               const Index& atmgeom_checked,
               const Vector& f_grid,
               const Index& cloudbox_on,
               const Index& cloudbox_checked,
               const Index& ppath_inside_cloudbox_do,
               const Vector& rte_pos,
               const Vector& rte_los,
               const Vector& rte_pos2,
               const Verbosity&) {
  // Basics
  //
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");

  ppath_agendaExecute(ws,
                      ppath,
                      ppath_lmax,
                      ppath_lraytrace,
                      rte_pos,
                      rte_los,
                      rte_pos2,
                      cloudbox_on,
                      ppath_inside_cloudbox_do,
                      f_grid,
                      ppath_agenda);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppathFromRtePos2(Workspace& ws,
                      Ppath& ppath,
                      Vector& rte_los,
                      Numeric& ppath_lraytrace,
                      const Agenda& ppath_step_agenda,
                      const Index& atmosphere_dim,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Tensor3& z_field,
                      const Vector& f_grid,
                      const Vector& refellipsoid,
                      const Matrix& z_surface,
                      const Vector& rte_pos,
                      const Vector& rte_pos2,
                      const Numeric& ppath_lmax,
                      const Numeric& za_accuracy,
                      const Numeric& pplrt_factor,
                      const Numeric& pplrt_lowest,
                      const Verbosity& verbosity) {
  //--- Check input -----------------------------------------------------------
  ARTS_USER_ERROR_IF (atmosphere_dim == 2,
        "2D atmospheres not yet handled. Support for negative"
        " zenith angles needed. Remind me (Patrick) to fix this.");
  //---------------------------------------------------------------------------

  // Geometric LOS from rte_pos to rte_pos2
  Vector rte_los_geom;
  rte_losGeometricFromRtePosToRtePos2(rte_los_geom,
                                      atmosphere_dim,
                                      lat_grid,
                                      lon_grid,
                                      refellipsoid,
                                      rte_pos,
                                      rte_pos2,
                                      verbosity);

  // Radius of rte_pos and rte_pos2
  const Numeric r1 =
      pos2refell_r(atmosphere_dim, refellipsoid, lat_grid, lon_grid, rte_pos) +
      rte_pos[0];
  const Numeric r2 =
      pos2refell_r(atmosphere_dim, refellipsoid, lat_grid, lon_grid, rte_pos2) +
      rte_pos2[0];

  // Geometric distance between rte_pos and rte_pos2, effective 2D-lat for
  // rte_pos and and Cartesian coordinates of rte_pos:
  Numeric l12, lat1 = 0, x1, y1 = 0, z1;
  if (atmosphere_dim <= 2) {
    if (atmosphere_dim == 2) {
      lat1 = rte_pos[1];
    }
    distance2D(l12, r1, lat1, r2, rte_pos2[1]);
    pol2cart(x1, z1, r1, lat1);
  } else {
    distance3D(l12, r1, rte_pos[1], rte_pos[2], r2, rte_pos2[1], rte_pos2[2]);
    sph2cart(x1, y1, z1, r1, rte_pos[1], rte_pos[2]);
  }

  // Define remaining variables used in the while-loop below
  //
  // Basic bookkeeping variables
  Numeric za_upp_limit = 180;
  Numeric za_low_limit = 0;
  //
  // Various variables associated with the ppath, and the point of the path
  // closest to the transmitter
  Ppath ppt;                     // "Test ppath"
  Index ip = -999;               // Index of closest ppath point
  Numeric xip, yip = 0, zip;     // Cartesian coords. of the closest ppath point
  Numeric dxip, dyip = 0, dzip;  // Cartesian LOS of the closest ppath point
  //
  // Data for the intersection of the l12-sphere
  Vector posc(max(Index(2), atmosphere_dim));
  Numeric rc, xc, yc = 0, zc;

  CREATE_OUT2;
  CREATE_OUT3;

  const Index maxiter = 99;
  Vector t_za(maxiter, -999), t_dza(maxiter, -999);
  Index it = -1;

  // Keep trying until ready or ground intersetion determined
  //
  bool ground = false;
  bool failed = false;
  Index ntries = 0;
  //
  while (true) {
    // Path for present rte_los (no cloudbox!)
    ppath_calc(ws,
               ppt,
               ppath_step_agenda,
               atmosphere_dim,
               p_grid,
               lat_grid,
               lon_grid,
               z_field,
               f_grid,
               refellipsoid,
               z_surface,
               0,
               ArrayOfIndex(0),
               rte_pos,
               rte_los,
               ppath_lmax,
               ppath_lraytrace,
               0,
               verbosity);

    // Find the point closest to rte_pos2, on the side towards rte_pos.
    // We do this by looking at the distance to rte_pos, that should be
    // as close to l12 as possible, but not exceed it.
    Numeric lip = 99e99;
    ip = ppt.np;
    //
    while (lip >= l12 && ip > 0) {
      ip--;
      if (atmosphere_dim <= 2) {
        distance2D(lip, r1, lat1, ppt.r[ip], ppt.pos(ip, 1));
      } else {
        distance3D(lip,
                   r1,
                   rte_pos[1],
                   rte_pos[2],
                   ppt.r[ip],
                   ppt.pos(ip, 1),
                   ppt.pos(ip, 2));
      }
    }

    Numeric za_new, daa = 0;

    // Surface intersection:
    // Not OK if the ground position is too far from rte_pos2.
    // (30 km selected to allow misses of smaller size when rte_pos2 is at
    // surface level, but surface interference never OK if rte_pos above TOA)
    if (ppath_what_background(ppt) == 2 && ip == ppt.np - 1 &&
        l12 - lip > 30e3) {
      za_new = rte_los[0] - 1;
      za_upp_limit = rte_los[0];
    }

    // Ppath OK
    else {
      // Estimate ppath at the distance of l12, and calculate size
      // of "miss" (measured in diffference in geometric angles)
      Vector los;
      Numeric dza;
      if (atmosphere_dim <= 2) {
        // Convert pos and los for point ip to cartesian coordinates
        poslos2cart(
            xip, zip, dxip, dzip, ppt.r[ip], ppt.pos(ip, 1), ppt.los(ip, 0));
        // Find where the extension from point ip crosses the l12
        // sphere: point c
        Numeric latc;
        line_circle_intersect(xc, zc, xip, zip, dxip, dzip, x1, z1, l12);
        cart2pol(rc, latc, xc, zc, ppt.pos(ip, 1), ppt.los(ip, 0));
        posc[1] = latc;
        posc[0] =
            rc - pos2refell_r(
                     atmosphere_dim, refellipsoid, lat_grid, lon_grid, posc);
      } else {
        // Convert pos and los for point ip to cartesian coordinates
        poslos2cart(xip,
                    yip,
                    zip,
                    dxip,
                    dyip,
                    dzip,
                    ppt.r[ip],
                    ppt.pos(ip, 1),
                    ppt.pos(ip, 2),
                    ppt.los(ip, 0),
                    ppt.los(ip, 1));
        // Find where the extension from point ip crosses the l12
        // sphere: point c
        Numeric latc, lonc;
        line_sphere_intersect(
            xc, yc, zc, xip, yip, zip, dxip, dyip, dzip, x1, y1, z1, l12);
        cart2sph(rc,
                 latc,
                 lonc,
                 xc,
                 yc,
                 zc,
                 ppt.pos(ip, 1),
                 ppt.pos(ip, 2),
                 ppt.los(ip, 0),
                 ppt.los(ip, 1));
        posc[1] = latc;
        posc[2] = lonc;
        posc[0] =
            rc - pos2refell_r(
                     atmosphere_dim, refellipsoid, lat_grid, lon_grid, posc);
      }
      //
      rte_losGeometricFromRtePosToRtePos2(los,
                                          atmosphere_dim,
                                          lat_grid,
                                          lon_grid,
                                          refellipsoid,
                                          rte_pos,
                                          posc,
                                          verbosity);
      //
      dza = los[0] - rte_los_geom[0];

      // Update bookkeeping variables
      it++;
      t_za[it] = rte_los[0];
      t_dza[it] = dza;
      //
      if (dza > 0 && rte_los[0] < za_upp_limit) {
        za_upp_limit = rte_los[0];
      } else if (dza < 0 && rte_los[0] > za_low_limit) {
        za_low_limit = rte_los[0];
      }

      // Ready ?
      if (abs(dza) <= za_accuracy) {
        break;
      } else if (za_upp_limit - za_low_limit <= za_accuracy / 10) {
        if (max(t_dza) < -10 * za_accuracy) {
          ground = true;
          out3 << "    Ground intersection determined !!!\n";
          break;
        } else {
          failed = true;
          out3 << "    Zenith angle search range closed !!!\n";
          break;
        }
      }
      // Catch non-convergence (just for extra safety, za-range should be
      // closed quicker than this)
      ntries += 1;
      if (ntries >= maxiter) {
        failed = true;
        out3 << "    Too many iterations !!!\n";
        break;
      }

      // Estimate new angle
      if (it < 1) {
        za_new = rte_los[0] - dza;
      } else {
        // Estimate new angle by linear regression over some of the
        // last calculations
        const Index nfit = min(it + 1, (Index)3);
        const Index i0 = it - nfit + 1;
        Vector p;
        linreg(p, t_za[Range(i0, nfit)], t_dza[Range(i0, nfit)]);
        za_new = -p[0] / p[1];
      }
      //
      if (atmosphere_dim == 3) {
        daa = los[1] - rte_los_geom[1];
      }
    }

    // Update rte_los. Use bisection of za_new is basically
    // identical to old angle, or is outside lower or upper
    // limit. Otherwise use reult of linear reg.
    if (std::isinf(za_new) || std::isnan(za_new) ||
        abs(za_new - rte_los[0]) < 0.99 * za_accuracy ||
        za_new <= za_low_limit || za_new >= za_upp_limit) {

      //Additional exit condition to avoid endless loop.
      if (abs(za_upp_limit-za_low_limit)<za_upp_limit*1e-15){
        ppath_init_structure(ppath, atmosphere_dim, 1);
        ppath_set_background(ppath, 0);
        return;
      }

      rte_los[0] = (za_low_limit + za_upp_limit) / 2;

    } else {
      rte_los[0] = za_new;
      if (atmosphere_dim == 3) {
        rte_los[1] -= daa;
        if (rte_los[1] < -180) {
          rte_los[1] += 360;
        } else if (rte_los[1] > 180) {
          rte_los[1] -= 360;
        }
      }
    }
  }  // while
  //--------------------------------------------------------------------------

  // If failed re-try with a shorter ppath_lraytrace, if not ending up with
  // a too small value.
  if (failed) {
    ppath_lraytrace /= pplrt_factor;

    if (ppath_lraytrace >= pplrt_lowest) {
      out2 << "  Re-start with ppath_lraytrace = " << ppath_lraytrace;
      ppathFromRtePos2(ws,
                       ppath,
                       rte_los,
                       ppath_lraytrace,
                       ppath_step_agenda,
                       atmosphere_dim,
                       p_grid,
                       lat_grid,
                       lon_grid,
                       z_field,
                       f_grid,
                       refellipsoid,
                       z_surface,
                       rte_pos,
                       rte_pos2,
                       ppath_lmax,
                       za_accuracy,
                       pplrt_factor,
                       pplrt_lowest,
                       verbosity);
    } else {
      ppath_init_structure(ppath, atmosphere_dim, 1);
      ppath_set_background(ppath, 0);
    }
    return;  // --->
  }

  // Create final ppath.
  // If ground intersection: Set to length 1 and ground background,
  // to flag non-OK path
  // Otherwise: Fill path and set background to transmitter

  if (ground) {
    ppath_init_structure(ppath, atmosphere_dim, 1);
    ppath_set_background(ppath, 2);
  }

  else {
    // Distance between point ip of ppt and posc
    Numeric ll;
    if (atmosphere_dim <= 2) {
      distance2D(ll, rc, posc[1], ppt.r[ip], ppt.pos(ip, 1));
    } else {
      distance3D(
          ll, rc, posc[1], posc[2], ppt.r[ip], ppt.pos(ip, 1), ppt.pos(ip, 2));
    }

    // Last point of ppt closest to rte_pos2. No point to add, maybe
    // calculate start_lstep and start_los:
    if (ip == ppt.np - 1) {
      ppath_init_structure(ppath, atmosphere_dim, ppt.np);
      ppath_copy(ppath, ppt, -1);
      if (ppath_what_background(ppath) == 1) {
        ppath.start_lstep = ll;
        Numeric d1, d2 = 0, d3;
        if (atmosphere_dim <= 2) {
          cart2poslos(d1,
                      d3,
                      ppath.start_los[0],
                      xc,
                      zc,
                      dxip,
                      dzip,
                      ppt.r[ip] * sin(DEG2RAD * ppt.los(ip, 0)),
                      ppt.pos(ip, 1),
                      ppt.los(ip, 0));
        } else {
          cart2poslos(d1,
                      d2,
                      d3,
                      ppath.start_los[0],
                      ppath.start_los[1],
                      xc,
                      yc,
                      zc,
                      dxip,
                      dyip,
                      dzip,
                      ppt.r[ip] * sin(DEG2RAD * ppt.los(ip, 0)),
                      xip,
                      yip,
                      zip,  // Added 161027,PE
                      ppt.pos(ip, 1),
                      ppt.pos(ip, 2),
                      ppt.los(ip, 0),
                      ppt.los(ip, 1));
        }
      }
    }
    // rte_pos2 inside the atmosphere (posc entered as end point)
    else {
      ppath_init_structure(ppath, atmosphere_dim, ip + 2);
      ppath_copy(ppath, ppt, ip + 1);
      //
      const Index i = ip + 1;
      if (atmosphere_dim <= 2) {
        cart2poslos(ppath.r[i],
                    ppath.pos(i, 1),
                    ppath.los(i, 0),
                    xc,
                    zc,
                    dxip,
                    dzip,
                    ppt.r[ip] * sin(DEG2RAD * ppt.los(ip, 0)),
                    ppt.pos(ip, 1),
                    ppt.los(ip, 0));
      } else {
        cart2poslos(ppath.r[i],
                    ppath.pos(i, 1),
                    ppath.pos(i, 2),
                    ppath.los(i, 0),
                    ppath.los(i, 1),
                    xc,
                    yc,
                    zc,
                    dxip,
                    dyip,
                    dzip,
                    ppt.r[ip] * sin(DEG2RAD * ppt.los(ip, 0)),
                    xip,
                    yip,
                    zip,  // Added 161027,PE
                    ppt.pos(ip, 1),
                    ppt.pos(ip, 2),
                    ppt.los(ip, 0),
                    ppt.los(ip, 1));
      }
      //
      ppath.pos(i, joker) = posc;
      ppath.lstep[i - 1] = ll;
      ppath.start_los = ppath.los(i, joker);

      // n by linear interpolation
      // Gets tripped when ll is very close to (slightly greater than) lstep (ISA)
      ARTS_ASSERT(ll < ppt.lstep[i - 1]);
      const Numeric w = ll / ppt.lstep[i - 1];
      ppath.nreal[i] = (1 - w) * ppt.nreal[i - 1] + w * ppt.nreal[i];
      ppath.ngroup[i] = (1 - w) * ppt.ngroup[i - 1] + w * ppt.ngroup[i];

      // Grid positions
      GridPos gp_lat, gp_lon;
      rte_pos2gridpos(ppath.gp_p[i],
                      gp_lat,
                      gp_lon,
                      atmosphere_dim,
                      p_grid,
                      lat_grid,
                      lon_grid,
                      z_field,
                      ppath.pos(i, Range(0, atmosphere_dim)));
      if (atmosphere_dim >= 2) {
        gridpos_copy(ppath.gp_lat[i], gp_lat);
        if (atmosphere_dim == 3) {
          gridpos_copy(ppath.gp_lon[i], gp_lon);
        }
      }
    }

    // Common stuff
    ppath_set_background(ppath, 9);
    ppath.start_pos = rte_pos2;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathPlaneParallel(Ppath& ppath,
                        const Index& atmosphere_dim,
                        const Tensor3& z_field,
                        const Matrix& z_surface,
                        const Index& cloudbox_on,
                        const ArrayOfIndex& cloudbox_limits,
                        const Index& ppath_inside_cloudbox_do,
                        const Vector& rte_pos,
                        const Vector& rte_los,
                        const Numeric& ppath_lmax,
                        const Verbosity&) {
  // This function is a WSM but it is normally only called from yCalc.
  // For that reason, this function does not repeat input checks that are
  // performed in yCalc, it only performs checks regarding the sensor
  // position and LOS.

  const Numeric z_sensor = rte_pos[0];
  const Numeric za_sensor = rte_los[0];
  const Index nz = z_field.npages();
  const Numeric z_toa = z_field(nz - 1, 0, 0);
  const bool above_toa = z_sensor > z_toa ? true : false;
  const Numeric z_end = above_toa ? z_toa : z_sensor;
  const Numeric dz2dl = abs(1 / cos(DEG2RAD * za_sensor));
  Index background = -99;

  // Basics checks of input
  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
                      "The function can only be used for 1D atmospheres.");
  chk_rte_pos(atmosphere_dim, rte_pos);
  chk_rte_los(atmosphere_dim, rte_los);
  ARTS_USER_ERROR_IF (ppath_inside_cloudbox_do && !cloudbox_on,
        "The WSV *ppath_inside_cloudbox_do* can only be set "
        "to 1 if also *cloudbox_on* is 1.");
  ARTS_USER_ERROR_IF (z_sensor < z_surface(0, 0),
       "The sensor is below the surface."
       "   altitude of sensor  : ", z_sensor, "\n"
       "   altitude of surface : ", z_surface(0, 0))
  ARTS_USER_ERROR_IF (abs(za_sensor - 90) < 0.1,
      "The zenith angle is ", za_sensor, "\n"
      "The method does not allow this. The zenith angle must deviate\n"
      "from 90 deg with at least 0.1 deg. That is, to be outside [89.9,90.1].")

  // Find end grid position
  GridPos gp_end;
  // To avoid compiler warnings, start to assuming above_toa
  gp_end.idx = nz - 2;
  gp_end.fd[0] = 1;
  gp_end.fd[1] = 0;
  if (!above_toa) {
    for (Index i = 0; i < nz - 1; i++) {
      if (z_sensor < z_field(i + 1, 0, 0)) {
        gp_end.idx = i;
        gp_end.fd[0] = (z_sensor - z_field(i, 0, 0)) /
                       (z_field(i + 1, 0, 0) - z_field(i, 0, 0));
        gp_end.fd[1] = 1 - gp_end.fd[0];
        break;
      }
    }
  }

  // Catch cases resulting in a ppath with 1 point
  bool path_to_follow = true;
  if (above_toa && za_sensor < 90) {
    // Path fully in space
    ppath_init_structure(ppath, atmosphere_dim, 1);
    background = 1;
    path_to_follow = false;
  } else if (z_sensor == z_surface(0, 0) && za_sensor > 90) {
    // On ground, looking down
    ppath_init_structure(ppath, atmosphere_dim, 1);
    background = 2;
    path_to_follow = false;
  } else if (cloudbox_on) {
    if (!ppath_inside_cloudbox_do &&
        z_sensor > z_field(cloudbox_limits[0], 0, 0) &&
        z_sensor < z_field(cloudbox_limits[1], 0, 0)) {
      // Inside cloud box
      ppath_init_structure(ppath, atmosphere_dim, 1);
      background = 4;
      path_to_follow = false;
    } else if ((z_sensor == z_field(cloudbox_limits[0], 0, 0) &&
                za_sensor > 90) ||
               (z_sensor == z_field(cloudbox_limits[1], 0, 0) &&
                za_sensor < 90)) {
      // Cloud box boundary
      ppath_init_structure(ppath, atmosphere_dim, 1);
      background = 3;
      path_to_follow = false;
    } else if (above_toa && cloudbox_limits[1] == nz - 1) {
      // Cloud box boundary is at TOA
      ppath_init_structure(ppath, atmosphere_dim, 1);
      background = 3;
      path_to_follow = false;
    }
  }

  // Determine ppath
  if (path_to_follow) {
    const Numeric max_dz = ppath_lmax > 0 ? ppath_lmax / dz2dl : 9e99;

    // Variables to describe each "break-point" of ppath. Point 0 is the end
    // point. Not all nz points are necessarily passed.
    ArrayOfIndex l_idx(nz);
    ArrayOfVector l_fd0(nz);
    ArrayOfVector l_z(nz);
    Index nptot = 0;

    // Determine number of ppath points in each layer
    {
      Numeric z = z_end;
      Index iout = -1;

      // Code similar, but for simplicity, we handle down- and
      // up-ward separately:
      if (za_sensor > 90)  // Downward-looking
      {
        // Here we go down to next pressure level (or the surface) in each
        // step. That is, if above surface, last point of step has fd[0]=0.

        // Put in end point
        iout++;
        nptot++;
        l_fd0[0].resize(1);
        l_z[0].resize(1);
        l_idx[0] = gp_end.idx;
        l_fd0[0][0] = gp_end.fd[0];
        l_z[0][0] = z_end;

        for (Index i = gp_end.idx; i >= 0 && background < 0; i--) {
          // Surface inside layer?
          Numeric dz_step;
          if (z_field(i, 0, 0) > z_surface(0, 0)) {
            dz_step = z - z_field(i, 0, 0);
          } else {
            dz_step = z - z_surface(0, 0);
            background = 2;
          }

          const Index np =
              dz_step <= max_dz ? 1 : Index(ceil(dz_step / max_dz));
          const Numeric dz = dz_step / Numeric(np);
          const Numeric dz_layer = z_field(i + 1, 0, 0) - z_field(i, 0, 0);

          // Update counters and resize
          iout++;
          nptot += np;
          l_fd0[iout].resize(np);
          l_z[iout].resize(np);

          // Intermediate points
          for (Index j = 0; j < np - 1; j++) {
            l_z[iout][j] = z - (Numeric(j) + 1) * dz;
            l_fd0[iout][j] = (l_z[iout][j] - z_field(i, 0, 0)) / dz_layer;
          }

          // End points handled seperately to avoid numerical problems
          l_idx[iout] = i;
          if (background == 2)  // Surface is reached
          {
            l_z[iout][np - 1] = z_surface(0, 0);
            l_fd0[iout][np - 1] =
                (l_z[iout][np - 1] - z_field(i, 0, 0)) / dz_layer;
          } else {
            l_z[iout][np - 1] = z_field(i, 0, 0);
            l_fd0[iout][np - 1] = 0;
            //
            if (cloudbox_on &&
                (i == cloudbox_limits[1] || i == cloudbox_limits[0])) {
              background = 3;
            }
          }

          // Update z
          z = z_field(i, 0, 0);
        }
      } else  // Upward-looking
      {
        // Here we have that first point of step has fd[0]=0, if not at
        // sensor
        for (Index i = gp_end.idx; i < nz && background < 0; i++) {
          Numeric dz_layer;
          Numeric dz_step;
          if (cloudbox_on && i != gp_end.idx &&
              (i == cloudbox_limits[0] ||
               i == cloudbox_limits[1])) {  // At an active cloudbox boundary
            dz_step = 0;
            dz_layer = 1;
            background = 3;
          } else if (i == nz - 1) {  // At TOA
            dz_step = 0;
            dz_layer = 1;
            background = 1;
          } else {
            dz_step = z_field(i + 1, 0, 0) - z;
            dz_layer = z_field(i + 1, 0, 0) - z_field(i, 0, 0);
          }

          const Index np =
              dz_step <= max_dz ? 1 : Index(ceil(dz_step / max_dz));
          const Numeric dz = dz_step / Numeric(np);

          // Update counters and resize
          iout++;
          nptot += np;
          l_fd0[iout].resize(np);
          l_z[iout].resize(np);

          // Start points handled seperately to avoid numerical problems
          if (i == gp_end.idx) {  // At sensor
            l_idx[iout] = i;
            l_z[iout][0] = z_sensor;
            l_fd0[iout][0] = gp_end.fd[0];
          } else if (i == nz - 1) {  // At TOA
            l_idx[iout] = i - 1;
            l_z[iout][0] = z_field(i, 0, 0);
            l_fd0[iout][0] = 1;
          } else {
            l_idx[iout] = i;
            l_z[iout][0] = z_field(i, 0, 0);
            l_fd0[iout][0] = 0;
          }

          // Intermediate points
          for (Index j = 1; j < np; j++) {
            l_z[iout][j] = z + Numeric(j) * dz;
            l_fd0[iout][j] = (l_z[iout][j] - z_field(i, 0, 0)) / dz_layer;
          }

          // Update z
          if (background < 0) {
            z = z_field(i + 1, 0, 0);
          }
        }
      }
    }

    ppath_init_structure(ppath, atmosphere_dim, nptot);

    // Fill ppath.pos(joker,0), ppath.gp_p and ppath.lstep
    Index iout = -1;
    Numeric z_last = -999;
    for (Index i = 0; i < nz; i++) {
      for (Index j = 0; j < l_z[i].nelem(); j++) {
        iout++;
        ppath.pos(iout, 0) = l_z[i][j];
        ppath.gp_p[iout].idx = l_idx[i];
        ppath.gp_p[iout].fd[0] = l_fd0[i][j];
        ppath.gp_p[iout].fd[1] = 1 - l_fd0[i][j];
        if (iout == 0) {
          z_last = ppath.pos(iout, 0);
        } else {
          ppath.lstep[iout - 1] = dz2dl * abs(z_last - l_z[i][j]);
          z_last = l_z[i][j];
        }
      }
    }
  }

  // Remaining data
  ppath_set_background(ppath, background);
  if (ppath.np == 1) {
    ppath.pos(0, 0) = z_end;
    ppath.gp_p[0] = gp_end;
  }
  ppath.pos(joker, 1) = 0;
  ppath.los(joker, 0) = za_sensor;
  ppath.constant = INFINITY;  // Not defined here as r = Inf
  ppath.r = INFINITY;
  ppath.start_pos[0] = ppath.pos(ppath.np - 1, 0);
  ppath.start_pos[1] = 0;
  ppath.start_los[0] = za_sensor;
  ppath.end_pos[0] = z_sensor;
  ppath.end_pos[1] = 0;
  ppath.end_los[0] = za_sensor;
  if (above_toa) {
    ppath.end_lstep = dz2dl * (z_sensor - z_toa);
  }
  ppath.nreal = 1;
  ppath.ngroup = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathStepByStep(Workspace& ws,
                     Ppath& ppath,
                     const Agenda& ppath_step_agenda,
                     const Index& ppath_inside_cloudbox_do,
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& z_field,
                     const Vector& f_grid,
                     const Vector& refellipsoid,
                     const Matrix& z_surface,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Vector& rte_pos,
                     const Vector& rte_los,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Verbosity& verbosity) {
  ppath_calc(ws,
             ppath,
             ppath_step_agenda,
             atmosphere_dim,
             p_grid,
             lat_grid,
             lon_grid,
             z_field,
             f_grid,
             refellipsoid,
             z_surface,
             cloudbox_on,
             cloudbox_limits,
             rte_pos,
             rte_los,
             ppath_lmax,
             ppath_lraytrace,
             ppath_inside_cloudbox_do,
             verbosity);
}

// FIXMEDOC@Richard  TRy to describe the meaning of ppath_field 

/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_fieldFromDownUpLimbGeoms(Workspace& ws,
                                    ArrayOfPpath& ppath_field,
                                    const Agenda& ppath_agenda,
                                    const Numeric& ppath_lmax,
                                    const Numeric& ppath_lraytrace,
                                    const Index& atmgeom_checked,
                                    const Tensor3& z_field,
                                    const Vector& f_grid,
                                    const Index& cloudbox_on,
                                    const Index& cloudbox_checked,
                                    const Index& ppath_inside_cloudbox_do,
                                    const Vector& rte_pos,
                                    const Vector& rte_los,
                                    const Vector& rte_pos2,
                                    const Vector& refellipsoid,
                                    const Index& atmosphere_dim,
                                    const Index& zenith_angles_per_position,
                                    const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (atmosphere_dim not_eq 1,
                      "Only for 1D atmospheres");
  ARTS_USER_ERROR_IF (refellipsoid[1] not_eq 0.0,
                      "Not allowed for non-spherical planets");
  ARTS_USER_ERROR_IF (ppath_lmax >= 0,
                      "Only allowed for long paths (ppath_lmax < 0)");

  // Positions and angles of interest
  const Numeric zmin = z_field(0, 0, 0);
  const Numeric zmax = z_field(z_field.npages() - 1, 0, 0);
  const Numeric r = refellipsoid[0];
  const Numeric above_surface_tangent =
      90 - RAD2DEG * std::acos((r) / (r + zmax)) + 1e-4;
  const Numeric below_surface_tangent =
      90 - RAD2DEG * std::acos((r) / (r + zmax)) - 1e-4;
  const Numeric top_tangent = 90 - 1e-4;

  ppath_field.resize(3 * zenith_angles_per_position);
  Index ppath_field_pos = 0;

  Vector zenith_angles(zenith_angles_per_position);

  // Upwards:
  nlinspace(zenith_angles, 0, 90, zenith_angles_per_position);
  Vector rte_pos_true = rte_pos;
  rte_pos_true[0] = zmin;
  Vector rte_los_true = rte_los;
  for (Index iz = 0; iz < zenith_angles_per_position; iz++) {
    rte_los_true[0] = zenith_angles[iz];

    ppathCalc(ws,
              ppath_field[ppath_field_pos],
              ppath_agenda,
              ppath_lmax,
              ppath_lraytrace,
              atmgeom_checked,
              f_grid,
              cloudbox_on,
              cloudbox_checked,
              ppath_inside_cloudbox_do,
              rte_pos_true,
              rte_los_true,
              rte_pos2,
              verbosity);

    ppath_field_pos++;
  }

  // Limb:
  nlinspace(zenith_angles,
            above_surface_tangent,
            top_tangent,
            zenith_angles_per_position);
  rte_pos_true[0] = zmax;
  for (Index iz = 0; iz < zenith_angles_per_position; iz++) {
    rte_los_true[0] = 180 - zenith_angles[iz];

    ppathCalc(ws,
              ppath_field[ppath_field_pos],
              ppath_agenda,
              ppath_lmax,
              ppath_lraytrace,
              atmgeom_checked,
              f_grid,
              cloudbox_on,
              cloudbox_checked,
              ppath_inside_cloudbox_do,
              rte_pos_true,
              rte_los_true,
              rte_pos2,
              verbosity);

    ppath_field_pos++;
  }

  // Downwards:
  nlinspace(
      zenith_angles, 0, below_surface_tangent, zenith_angles_per_position);
  for (Index iz = 0; iz < zenith_angles_per_position; iz++) {
    rte_los_true[0] = 180 - zenith_angles[iz];

    ppathCalc(ws,
              ppath_field[ppath_field_pos],
              ppath_agenda,
              ppath_lmax,
              ppath_lraytrace,
              atmgeom_checked,
              f_grid,
              cloudbox_on,
              cloudbox_checked,
              ppath_inside_cloudbox_do,
              rte_pos_true,
              rte_los_true,
              rte_pos2,
              verbosity);

    ppath_field_pos++;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepGeometric(  // WS Output:
    Ppath& ppath_step,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    const Vector& refellipsoid,
    const Matrix& z_surface,
    const Numeric& ppath_lmax,
    const Verbosity&) {
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks
  // here.

  // A call with background set, just wants to obtain the refractive index for
  // complete ppaths consistent of a single point.
  if (!ppath_what_background(ppath_step)) {
    if (atmosphere_dim == 1) {
      ppath_step_geom_1d(ppath_step,
                         z_field(joker, 0, 0),
                         refellipsoid,
                         z_surface(0, 0),
                         ppath_lmax);
    }

    else if (atmosphere_dim == 2) {
      ppath_step_geom_2d(ppath_step,
                         lat_grid,
                         z_field(joker, joker, 0),
                         refellipsoid,
                         z_surface(joker, 0),
                         ppath_lmax);
    }

    else if (atmosphere_dim == 3) {
      ppath_step_geom_3d(ppath_step,
                         lat_grid,
                         lon_grid,
                         z_field,
                         refellipsoid,
                         z_surface,
                         ppath_lmax);
    }

    else {
      ARTS_USER_ERROR ( "The atmospheric dimensionality must be 1-3.");
    }
  }

  else {
    ARTS_ASSERT(ppath_step.np == 1);
    ppath_step.nreal[0] = 1;
    ppath_step.ngroup[0] = 1;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepRefractionBasic(Workspace& ws,
                               Ppath& ppath_step,
                               const Agenda& refr_index_air_agenda,
                               const Index& atmosphere_dim,
                               const Vector& p_grid,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const Tensor3& z_field,
                               const Tensor3& t_field,
                               const Tensor4& vmr_field,
                               const Vector& refellipsoid,
                               const Matrix& z_surface,
                               const Vector& f_grid,
                               const Numeric& ppath_lmax,
                               const Numeric& ppath_lraytrace,
                               const Verbosity&) {
  // Input checks here would be rather costly as this function is called
  // many times.
  ARTS_ASSERT(ppath_lraytrace > 0);

  // A call with background set, just wants to obtain the refractive index for
  // complete ppaths consistent of a single point.
  if (!ppath_what_background(ppath_step)) {
    if (atmosphere_dim == 1) {
      ppath_step_refr_1d(ws,
                         ppath_step,
                         p_grid,
                         z_field,
                         t_field,
                         vmr_field,
                         f_grid,
                         refellipsoid,
                         z_surface(0, 0),
                         ppath_lmax,
                         refr_index_air_agenda,
                         "linear_basic",
                         ppath_lraytrace);
    } else if (atmosphere_dim == 2) {
      ppath_step_refr_2d(ws,
                         ppath_step,
                         p_grid,
                         lat_grid,
                         z_field,
                         t_field,
                         vmr_field,
                         f_grid,
                         refellipsoid,
                         z_surface(joker, 0),
                         ppath_lmax,
                         refr_index_air_agenda,
                         "linear_basic",
                         ppath_lraytrace);
    } else if (atmosphere_dim == 3) {
      ppath_step_refr_3d(ws,
                         ppath_step,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         z_field,
                         t_field,
                         vmr_field,
                         f_grid,
                         refellipsoid,
                         z_surface,
                         ppath_lmax,
                         refr_index_air_agenda,
                         "linear_basic",
                         ppath_lraytrace);
    } else {
      ARTS_USER_ERROR ( "The atmospheric dimensionality must be 1-3.");
    }
  }

  else {
    ARTS_ASSERT(ppath_step.np == 1);
    if (atmosphere_dim == 1) {
      get_refr_index_1d(ws,
                        ppath_step.nreal[0],
                        ppath_step.ngroup[0],
                        refr_index_air_agenda,
                        p_grid,
                        refellipsoid,
                        z_field,
                        t_field,
                        vmr_field,
                        f_grid,
                        ppath_step.r[0]);
    } else if (atmosphere_dim == 2) {
      get_refr_index_2d(ws,
                        ppath_step.nreal[0],
                        ppath_step.ngroup[0],
                        refr_index_air_agenda,
                        p_grid,
                        lat_grid,
                        refellipsoid,
                        z_field,
                        t_field,
                        vmr_field,
                        f_grid,
                        ppath_step.r[0],
                        ppath_step.pos(0, 1));
    } else {
      get_refr_index_3d(ws,
                        ppath_step.nreal[0],
                        ppath_step.ngroup[0],
                        refr_index_air_agenda,
                        p_grid,
                        lat_grid,
                        lon_grid,
                        refellipsoid,
                        z_field,
                        t_field,
                        vmr_field,
                        f_grid,
                        ppath_step.r[0],
                        ppath_step.pos(0, 1),
                        ppath_step.pos(0, 2));
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losGeometricFromRtePosToRtePos2(Vector& rte_los,
                                         const Index& atmosphere_dim,
                                         const Vector& lat_grid,
                                         const Vector& lon_grid,
                                         const Vector& refellipsoid,
                                         const Vector& rte_pos,
                                         const Vector& rte_pos2,
                                         const Verbosity&) {
  // Check input
  chk_rte_pos(atmosphere_dim, rte_pos);
  chk_rte_pos(atmosphere_dim, rte_pos2, true);

  // Radius of rte_pos and rte_pos2
  const Numeric r1 =
      pos2refell_r(atmosphere_dim, refellipsoid, lat_grid, lon_grid, rte_pos) +
      rte_pos[0];
  const Numeric r2 =
      pos2refell_r(atmosphere_dim, refellipsoid, lat_grid, lon_grid, rte_pos2) +
      rte_pos2[0];

  // Remaining polar and cartesian coordinates of rte_pos
  Numeric lat1, lon1 = 0, x1, y1 = 0, z1;
  // Cartesian coordinates of rte_pos2
  Numeric x2, y2 = 0, z2;
  //
  if (atmosphere_dim == 1) {
    // Latitude distance implicitly checked by chk_rte_pos
    lat1 = 0;
    pol2cart(x1, z1, r1, lat1);
    pol2cart(x2, z2, r2, rte_pos2[1]);
  } else if (atmosphere_dim == 2) {
    lat1 = rte_pos[1];
    pol2cart(x1, z1, r1, lat1);
    pol2cart(x2, z2, r2, rte_pos2[1]);
  } else {
    lat1 = rte_pos[1];
    lon1 = rte_pos[2];
    sph2cart(x1, y1, z1, r1, lat1, lon1);
    sph2cart(x2, y2, z2, r2, rte_pos2[1], rte_pos2[2]);
  }

  // Geometrical LOS to transmitter
  Numeric za, aa;
  //
  los2xyz(za, aa, r1, lat1, lon1, x1, y1, z1, x2, y2, z2);
  //
  if (atmosphere_dim == 3) {
    rte_los.resize(2);
    rte_los[0] = za;
    rte_los[1] = aa;
  } else {
    rte_los.resize(1);
    rte_los[0] = za;
    if (atmosphere_dim == 2 && aa < 0)  // Should 2D-za be negative?
    {
      rte_los[0] = -za;
    }
  }
}
