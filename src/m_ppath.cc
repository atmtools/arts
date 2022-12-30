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
 * @file   m_ppath.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2002-05-08 
 *
 * @brief  Workspace functions releated to propagation paths variables.
 *
 * The file includes special functions to set the sensor position and LOS,
 * and functions for calculation of propagation paths.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "m_general.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "refraction.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric NAT_LOG_2=Constant::ln_2;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void dlosDiffOfLos(Matrix& dlos,
                   const Vector& ref_los,
                   const Matrix& other_los,
                   const Verbosity&) {
  ARTS_USER_ERROR_IF (ref_los.nelem() != 2,
                      "*ref_los* must have two columns.");
  ARTS_USER_ERROR_IF (other_los.ncols() != 2,
                      "*other_los* must have two columns.");

  const Index nlos = other_los.nrows();

  dlos.resize(nlos, 2);

  for (Index i = 0; i < nlos; i++) {
    diff_za_aa(dlos(i, 0),
               dlos(i, 1),
               ref_los[0],
               ref_los[1],
               other_los(i, 0),
               other_los(i, 1));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dlosGauss(Matrix& dlos,
               Vector& dlos_weight_vector,
               const Numeric& fwhm_deg,
               const Index& ntarget,
               const Index& include_response_in_weight,
               const Verbosity&) {
  const Index n_per_layer = 3;
  
  // Use FWHM and sigma in radians to get solid angles right
  const Numeric fwhm = DEG2RAD *fwhm_deg;
  const Numeric si = fwhm / (2 * sqrt(2 * NAT_LOG_2));

  // Cumulative distribution of Gauss weighted area, as a function of radius x
  Vector xp, cx;
  {
    VectorLinSpace(xp, 0, 1.5*fwhm, 0.02*fwhm, Verbosity());
    Vector gx;
    VectorGaussian(gx, xp, 0, -1.0, fwhm, Verbosity());
    gx *= xp;  // Weight with radius, no need to include pi
    const Index np = gx.nelem();
    cx.resize(np);
    cumsum(cx, gx);
    cx /= cx[np-1];
  }
    
  // Number of layers (not including (0,0)), and total number of points
  const Index nlayers = (Index) round((ntarget - 1) / n_per_layer);
  const Index npoints = 1 + nlayers * n_per_layer;

  // Distribution of the layers w.r.t. cumulative distribution
  Vector cp(nlayers);
  {
    const Numeric nterm = 1 / (Numeric) npoints;
    for (Index i=0; i < nlayers; ++i)
      cp[i] = nterm + (1 - nterm) * ((Numeric)i+0.5)/(Numeric)nlayers;
  }
  
  // Radii of layers, obtained by interpolating xp(cx) to cp
  Vector r(nlayers);
  {
    ArrayOfGridPos gp(nlayers);
    gridpos(gp, cx, cp);
    Matrix itw(nlayers, 2);
    interpweights(itw, gp);
    interp(r, itw, xp, gp);
  }

  // Factor to rescale Gauss(x) from 1D to 2D (along y=0)
  const Numeric scfac = 1 / (sqrt(2 * PI) * si);
  
  // Calculate dlos and weights
  dlos.resize(npoints, 2);
  dlos(0, joker) = 0.0;
  //
  dlos_weight_vector.resize(npoints);
  dlos_weight_vector = 1.0 / (Numeric) npoints;
  // If include_response_in_weight, all weights equal.
  // Otherwise, we need to divide with value of 2D Gauss for radius:
  Vector gv(1);
  if (!include_response_in_weight) {
    VectorGaussian(gv, Vector(1, 0.0), 0, -1.0, fwhm, Verbosity());
    dlos_weight_vector[0] = dlos_weight_vector[0] / (scfac * gv[0]);
    gv.resize(nlayers);
    VectorGaussian(gv, r, 0, -1.0, fwhm, Verbosity());
  }
  //
  const Numeric dalpha = 2 * PI / (Numeric) n_per_layer;
  const ArrayOfIndex shift = {0, 2, 4, 1, 3};
  //
  Index n = 0;
  for (Index i=0; i<nlayers; ++i) {
    Numeric alpha0;
    if (nlayers <= 3)
      alpha0 = dalpha * ((Numeric)(i%2) / 2.0);
    else
      alpha0 = dalpha * ((Numeric) shift[i%5] / 5.0);
    for (Index angle=0; angle<n_per_layer; ++angle) {
      const Numeric alpha = alpha0 + (Numeric) angle * dalpha;
      ++n;
      dlos(n, 0) = r[i] * cos(alpha);
      dlos(n, 1) = r[i] * sin(alpha);
      if (!include_response_in_weight) {
        dlos_weight_vector[n] = dlos_weight_vector[n] / (scfac * gv[i]);
      }
    }
  }
  dlos *= RAD2DEG;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dlosUniform(Matrix& dlos,
                 Vector& dlos_weight_vector,
                 const Numeric& width,
                 const Index& npoints,
                 const Index& crop_circular,
                 const Verbosity&) {
  ARTS_USER_ERROR_IF(npoints < 2, "GIN npoints must be > 1.");

  // Edges of angular grid
  Vector grid_edges;
  Numeric hwidth = width / 2.0;
  nlinspace(grid_edges, -hwidth, hwidth, npoints + 1);

  // Angular grid
  const Numeric spacing = grid_edges[1] - grid_edges[0];
  Vector grid;
  hwidth -= spacing / 2.0;
  nlinspace(grid, -hwidth, hwidth, npoints);

  // Square set
  dlos.resize(npoints * npoints, 2);
  dlos_weight_vector.resize(npoints * npoints);
  //
  grid_edges *= DEG2RAD; 
  const Numeric fac = DEG2RAD * spacing;
  for (Index z = 0; z < npoints; ++z) {
    const Numeric solid_angle = fac * (sin(grid_edges[z+1]) - sin(grid_edges[z])); 
    for (Index a = 0; a < npoints; ++a) {
      const Index i = a * npoints + z;
      dlos(i, 0) = grid[z];
      dlos(i, 1) = grid[a];
      dlos_weight_vector[i] = solid_angle;  
    }
  }

  // Crop to circular?
  if (crop_circular) {
    // Pick out points inside radius (with special treatment of npoints=3)
    Matrix dlos_tmp(dlos.nrows(), 2);
    Vector sa_tmp(dlos_weight_vector.nelem());
    const Numeric r = width / 2.0 * (npoints != 3 ? 1 : 0.8);
    //
    Index n = 0;
    for (Index i = 0; i < npoints * npoints; ++i) {
      if (sqrt(pow(dlos(i,0), 2.0) + pow(dlos(i,1), 2.0)) <= r) {
        dlos_tmp(n, joker) = dlos(i, joker);
        sa_tmp[n] = dlos_weight_vector[i];
        ++n;
      }
    }

    // Reset output variables 
    dlos = dlos_tmp(Range(0, n), joker);
    dlos_weight_vector = sa_tmp[Range(0, n)];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posEndOfPpath(Vector& geo_pos,
                       const Ppath& ppath,
                       const Verbosity& verbosity) {
  geo_pos.resize(5);
  geo_pos = NAN;

  geo_pos[Range(0, ppath.pos.ncols())] =
      ppath.pos(ppath.np - 1, Range(0, ppath.pos.ncols()));
  geo_pos[Range(3, ppath.los.ncols())] =
      ppath.los(ppath.np - 1, Range(0, ppath.los.ncols()));

  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posLowestAltitudeOfPpath(Vector& geo_pos,
                                  const Ppath& ppath,
                                  const Verbosity& verbosity) {
  geo_pos.resize(5);
  geo_pos = NAN;

  // Take first point of ppath as first guess
  geo_pos[Range(0, ppath.pos.ncols())] =
      ppath.pos(0, Range(0, ppath.pos.ncols()));
  geo_pos[Range(3, ppath.los.ncols())] =
      ppath.los(0, Range(0, ppath.los.ncols()));

  for (Index i = 1; i < ppath.np; i++) {
    if (ppath.pos(i, 0) < geo_pos[0]) {
      geo_pos[Range(0, ppath.pos.ncols())] =
          ppath.pos(i, Range(0, ppath.pos.ncols()));
      geo_pos[Range(3, ppath.los.ncols())] =
          ppath.los(i, Range(0, ppath.los.ncols()));
    }
  }

  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void geo_posWherePpathPassesZref(Vector& geo_pos,
                                 const Ppath& ppath,
                                 const Numeric& z_ref,
                                 const Verbosity& verbosity) {
  geo_pos.resize(5);
  geo_pos = NAN;

  bool found = false;
  Index ihit = 0;
  bool above = false;

  if (ppath.pos(0, 0) >= z_ref) {
    above = true;
  }

  while (!found && ihit < ppath.np - 1) {
    ihit += 1;
    if (above && ppath.pos(ihit, 0) < z_ref) {
      found = true;
    } else if (!above && ppath.pos(ihit, 0) >= z_ref) {
      found = true;
    }
  }

  if (found) {
    geo_pos[0] = z_ref;

    // Make a simple linear interpolation to determine lat and lon
    const Numeric w = (z_ref - ppath.pos(ihit - 1, 0)) /
                      (ppath.pos(ihit, 0) - ppath.pos(ihit - 1, 0));

    geo_pos[3] = w * ppath.los(ihit, 0) + (1 - w) * ppath.los(ihit - 1, 0);

    if (ppath.pos.ncols() > 1) {
      geo_pos[1] = w * ppath.pos(ihit, 1) + (1 - w) * ppath.pos(ihit - 1, 1);

      if (ppath.pos.ncols() > 2) {
        geo_pos[2] = w * ppath.pos(ihit, 2) + (1 - w) * ppath.pos(ihit - 1, 2);
        geo_pos[4] = w * ppath.los(ihit, 1) + (1 - w) * ppath.los(ihit - 1, 1);
      }
    }
  }

  CREATE_OUT2;
  out2 << "  Sets geo-position to:\n" << geo_pos;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void losAddLosAndDlos(Matrix& new_los,
                      const Vector& ref_los,
                      const Matrix& dlos,
                      const Verbosity&) {
  ARTS_USER_ERROR_IF (ref_los.nelem() != 2,
                      "*ref_los* must have two columns.");
  ARTS_USER_ERROR_IF (dlos.ncols() != 2,
                      "*dlos* must have two columns.");

  const Index nlos = dlos.nrows();

  new_los.resize(nlos, 2);

  for (Index i = 0; i < nlos; i++) {
    add_za_aa(new_los(i, 0),
              new_los(i, 1),
              ref_los[0],
              ref_los[1],
              dlos(i, 0),
              dlos(i, 1));
  }
}

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
void ppathCalcFromAltitude(Workspace& ws,
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
                           const Numeric& altitude,
                           const Numeric& accuracy,
                           const Verbosity& verbosity) {
  ppathCalc(ws,
            ppath,
            ppath_agenda,
            ppath_lmax,
            ppath_lraytrace,
            atmgeom_checked,
            f_grid,
            cloudbox_on,
            cloudbox_checked,
            ppath_inside_cloudbox_do,
            rte_pos,
            rte_los,
            rte_pos2,
            verbosity);

  // Iterate until converging at altitude of interest
  Index pos = first_pos_before_altitude(ppath, altitude);
  Numeric lmax = ppath_lmax;
  while (true) {
    lmax *= 0.5;

    ppathCalc(ws,
              ppath,
              ppath_agenda,
              lmax,
              ppath_lraytrace,
              atmgeom_checked,
              f_grid,
              cloudbox_on,
              cloudbox_checked,
              ppath_inside_cloudbox_do,
              ppath.dim == 1 ? ppath.pos(pos, 0) : ppath.pos(pos, joker),
              ppath.dim == 1 ? ppath.los(pos, 0) : ppath.los(pos, joker),
              rte_pos2,
              verbosity);
    pos = first_pos_before_altitude(ppath, altitude);

    if (std::abs(ppath.pos(pos, 0) - altitude) < accuracy) {
      ppathCalc(ws,
                ppath,
                ppath_agenda,
                ppath_lmax,
                ppath_lraytrace,
                atmgeom_checked,
                f_grid,
                cloudbox_on,
                cloudbox_checked,
                ppath_inside_cloudbox_do,
                ppath.dim == 1 ? ppath.pos(pos, 0) : ppath.pos(pos, joker),
                ppath.dim == 1 ? ppath.los(pos, 0) : ppath.los(pos, joker),
                rte_pos2,
                verbosity);
      break;
    }
  }
}


/* A help function for *ppathFixedLstep*. Determines the radius of the path,
   raius of ellipsoid and surface altitiude for (lat,lon) at a distance l from
   the sensor.
 */
void surf_radius_at_l(Numeric& r,
                      Numeric& r_e,
                      Numeric& z_surf,
                      const Index& atmosphere_dim,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Vector& refellipsoid,
                      const Matrix& z_surface,
                      const Numeric& x0,
                      const Numeric& y0,
                      const Numeric& z0,
                      const Numeric& dx,
                      const Numeric& dy,
                      const Numeric& dz, 
                      const Numeric& l,
                      const Numeric& lat0,
                      const Numeric& lon0,
                      const Numeric& za0,
                      const Numeric& aa0) {
  if (atmosphere_dim == 1) {   // 1D
    const Numeric x = x0+l*dx, y = y0+l*dy, z = z0+l*dz;
    r = sqrt( x*x + y*y + z*z);
    r_e = refellipsoid[0];
    z_surf = z_surface(0,0);
    
  } else {  // 2D and 3D
    
    // Go to spherical coordinates
    Numeric lat, lon;
    cart2sph(r, lat, lon, x0+l*dx, y0+l*dy, z0+l*dz, lat0, lon0, za0, aa0);

    // Latitude grid position
    ARTS_USER_ERROR_IF (lat < lat_grid[0],
        "Search of surface intersection ended up ", lat_grid[0]-lat,
        " degrees below start of *lat_grid*. You need to expand the grid.")
    ARTS_USER_ERROR_IF (lat > last(lat_grid),
        "Search of surface intersection ended up ", lat-last(lat_grid),
        " degrees above end of *lat_grid*. You need to expand the grid.")
    
    GridPos gp_lat, gp_lon;
    gridpos(gp_lat, lat_grid, lat);

    // Interpolate to get ellipsoid radius and surface altitude
    // If 3D we also need to calculate lon grid position
    r_e = refell2d(refellipsoid, lat_grid, gp_lat);
    if (atmosphere_dim==2) {
      Vector itw(2);
      interpweights(itw, gp_lat);
      z_surf = interp(itw, z_surface(joker, 0), gp_lat);
    } else {
      const Numeric lonmax = last(lon_grid);
      resolve_lon(lon, lon_grid[0], lonmax);
      ARTS_USER_ERROR_IF (lon < lon_grid[0],
          "Search of surface intersection ended up ", lon_grid[0]-lon,
          " degrees below start of *lon_grid*. You need to expand the grid.")
      ARTS_USER_ERROR_IF (lon > lonmax,
          "Search of surface intersection ended up ", lon-lonmax,
           " degrees above end of *lon_grid*. You need to expand the grid.")
      
      gridpos(gp_lon, lon_grid, lon);
      Vector itw(4);
      interpweights(itw, gp_lat, gp_lon);
      z_surf = interp(itw, z_surface, gp_lat, gp_lon);
    }
  }
}
                     

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathFixedLstep(Ppath& ppath,
                     const Index& atmfields_checked,
                     const Index& atmgeom_checked,
                     const Index& atmosphere_dim,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& z_field,
                     const Vector& refellipsoid,
                     const Matrix& z_surface,
                     const Index& cloudbox_on,
                     const Vector& rte_pos,
                     const Vector& rte_los,
                     const Numeric& ppath_lmax,
                     const Index& za_scale,
                     const Numeric& z_coarse,
                     const Numeric& l_coarse,
                     const Verbosity&) {
  // Basics checks of input
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  chk_rte_pos(atmosphere_dim, rte_pos);
  chk_rte_los(atmosphere_dim, rte_los);
  ARTS_USER_ERROR_IF (cloudbox_on,
        "This method does not yet handle an active cloudbox.");
  ARTS_USER_ERROR_IF (ppath_lmax <= 0,
        "This method requires that *ppath_lmax* > 0.");
  ARTS_USER_ERROR_IF (atmosphere_dim == 1 && refellipsoid[1] > 1e-7,
    "For 1D, second element of *refellipsoid* must be 0.");

  // Set lat etc according to atmosphere_dim
  const Numeric lon0 = atmosphere_dim == 3 ? rte_pos[2] : 0;
  const Numeric lat0 = atmosphere_dim >= 2 ? rte_pos[1] : 0;
  const Numeric r0 = refell2r(refellipsoid, lat0) + rte_pos[0];
  const Numeric za0 = abs(rte_los[0]);
        Numeric aa0 = atmosphere_dim == 3 ? rte_los[1] : 0;
  if (atmosphere_dim == 2 && rte_los[0]<0) { aa0 = 180; } // 2D special

  // Obtain a top of the atmosphere altitude
  // To make sure that it represents a point inside the atmosphere, we take the lowest
  // altitude of top pressure level
  const Numeric z_toa = min( z_field(z_field.npages()-1,joker,joker) );

  // Set step lengths to use
  const Numeric lstep = ppath_lmax * (za_scale ? 1/abs(cos(DEG2RAD*za0)) : 1); 
  const Numeric lcoarse = l_coarse * (za_scale ? 1/abs(cos(DEG2RAD*za0)) : 1); 

  // ECEF pos and los
  Numeric x0, y0, z0, dx, dy, dz;
  poslos2cart(x0, y0, z0, dx, dy, dz, r0, lat0, lon0, za0, aa0);
  
  // Length to z_coarse
  //
  // We do this, here and below, by adding the search altitude to the major
  // axis of the ellipsoid. This is approxinative. In theory, the eccentricity
  // should be adopted. More important is that in ARTS the radius varies
  // lineraly between points of the latitude grid, while
  // line_refellipsoid_intersect operates with a fully analytical description
  // of the ellipsoid and this causes some inconsistency. That is, the found
  // length is not exact from ARTS's perspective. 
  //
  Vector rell = refellipsoid;  
  Numeric l2coarse = -1;      
  if (z_coarse >= 0) {         
    rell[0] += z_coarse;      
    line_refellipsoid_intersect(l2coarse, rell, x0, y0, z0, dx, dy, dz);
  }
  
  // Create a vector with the distance from the sensor for each ppath point
  //
  Vector lvec(1); lvec[0] = 0;
  Index background = 1;   // Index of radiative background. 1 = space
  //
  // Upward
  if (za0 < 90) {

    if (rte_pos[0] >= z_toa) {
      // If looking up we get an empty path
      // lvec initiated to represent this

    } else {
      // We are inside looking up. Lengths go from 0 to distance to TOA
      Numeric l2toa;
      rell[0] = refellipsoid[0] + z_toa;   
      line_refellipsoid_intersect(l2toa, rell, x0, y0, z0, dx, dy, dz);

      // Check that sensor actually is above the surface
      Numeric r, r_e, z_surf;
      surf_radius_at_l(r, r_e, z_surf, atmosphere_dim, lat_grid, lon_grid,
                       refellipsoid, z_surface, x0, y0, z0, dx, dy, dz,
                       0, lat0, lon0, za0, aa0);
      ARTS_USER_ERROR_IF (r < r_e+z_surf,
        "The sensor is ", r_e+z_surf - r, " m below the surface")

      // Create vector with lengths, considering z_coarse
      if (z_coarse < 0) {
        linspace(lvec, 0, l2toa, lstep);
      } else {
        if (rte_pos[0] > z_coarse) {
          linspace(lvec, 0, l2toa, lcoarse);
        } else {
          Vector l1, l2;
          linspace(l1, 0, l2coarse, lstep);
          linspace(l2, last(l1)+lstep, l2toa, lcoarse);
          const Index n1=l1.nelem(), n2=l2.nelem();
          lvec.resize(n1+n2);
          lvec[Range(0,n1)] = l1;
          lvec[Range(n1,n2)] = l2;
        }
      }
    }
  }

  // Downward
  else {
    // Determine lowest surface altitude
    // We can directly find length to the surface if its altitude is constant
    // so we also look for this. Always true for 1D!
    Numeric z_surf_min = atmosphere_dim == 1 ? z_surface(0,0) : 1e99;
    bool surface_found = true;
    if (atmosphere_dim == 2) {
      for (Index i=1; i<lat_grid.nelem(); i++) {
        if (z_surface(i,0) < z_surf_min)
          z_surf_min = z_surface(i,0);
        if (z_surface(i,0) != z_surface(0,0))
          surface_found = false;
      }
    } else if (atmosphere_dim == 3) {
      for (Index i=1; i<lat_grid.nelem(); i++) {
        for (Index j=1; j<lon_grid.nelem(); j++) {
          if (z_surface(i,j) < z_surf_min)
            z_surf_min = z_surface(i,j);
          if (z_surface(i,j) != z_surface(0,0))
            surface_found = false;
        }
      }
    }

    // Determine length to z_surf_min (l2s)
    Numeric l2s;
    rell[0] = refellipsoid[0] + z_surf_min;
    line_refellipsoid_intersect(l2s, rell, x0, y0, z0, dx, dy, dz);

    // No intersection with surface
    if (l2s < 0) {
      if (!surface_found) {
        ARTS_USER_ERROR (
                            "Cases of limb sounding type are only allowed "
                            "together with constant surface altitude.");
      } else {
        // Make function to find proper (geodetic) tangent point to finish this
        ARTS_USER_ERROR ( "Limb sounding not yet handled.");
      }

    // Ensured intersection with surface
    } else {

      background = 2;

      // If we need to search for the surface, we do this in two steps
      // 1. Move upward until we have a point at or above surface
      // 2. Use intersection to find length to surface
      bool inphase1 = true;
      bool underground = false;  // If first l2s happen to be above surface 
      const Numeric lmoveup = 1e3;
      const Numeric dr = 0.01;    // When we stop
      Numeric llong = l2s;
      Numeric lshort = -1;
      while (!surface_found) {
        // Determine where we have the surface at l2s
        Numeric r, r_e, z_surf;
        surf_radius_at_l(r, r_e, z_surf, atmosphere_dim, lat_grid, lon_grid,
                         refellipsoid, z_surface, x0, y0, z0, dx, dy, dz,
                         l2s, lat0, lon0, za0, aa0);
        
        // Compare radii and update iteration variables
        const Numeric r_s = r_e + z_surf;
        if (abs(r-r_s) <= dr) {
          surface_found = true;
        } else if (inphase1) {
          if (r < r_s) {
            underground = true;
            l2s -= lmoveup;
          } else {
            if (underground) {
              inphase1 = false;
              lshort = l2s;
              l2s = (lshort+llong)/2;
            } else {
              l2s += 0.1*lmoveup;  // If we end up here the start l2s was too
              llong = l2s;         // short. Can happen due to inconsistencies, 
            }                      // as explained above.
          }          
        } else {
          if (r > r_s) {
            lshort = l2s;
          } else {
            llong = l2s;
          }          
          l2s = (lshort+llong)/2;
        }
      }  // while

      // Distance to TOA (same code as above for upward)
      Numeric l2toa;
      if (rte_pos[0] <= z_toa) {
        l2toa = 0;
      } else {
        rell[0] = refellipsoid[0] + z_toa;   
        line_refellipsoid_intersect(l2toa, rell, x0, y0, z0, dx, dy, dz);
      }
      
      // Create vector with lengths, considering z_coarse
      // lvec shall start exactly at l2s
      if (l2s < lstep) {
        linspace(lvec, 0, l2s, l2s);
      } else if (z_coarse < 0) {
        linspace(lvec, l2toa, l2s, lstep);
        lvec += l2s - last(lvec);
      } else {
        if (rte_pos[0] < z_coarse) {
          linspace(lvec, l2toa, l2s, lstep);
          lvec += l2s - last(lvec);
        } else {
          Vector l1, l2;
          linspace(l1, l2coarse, l2s, lstep);
          l1 += l2s - last(l1);
          linspace(l2, l2toa, l1[0]-lstep, lcoarse);
          l2 += l1[0]-lstep - last(l2);
          const Index n1=l1.nelem(), n2=l2.nelem();
          lvec.resize(n1+n2);
          lvec[Range(0,n2)] = l2;
          lvec[Range(n2,n1)] = l1;
        }
      }
      
      
    } // Surface intersection
  }  // Up/down

  // Create ppath structure
  ppath_init_structure(ppath, atmosphere_dim, lvec.nelem());  
  ppath.constant = geometrical_ppc(r0, abs(rte_los[0]));
  ppath_set_background(ppath, background);
  ppath.end_los[joker] = rte_los[joker];   // end_pos done below
  ppath.end_lstep = lvec[0];
  ppath.nreal = 1;
  ppath.ngroup = 1;
  
  // Empty ppath
  if (ppath.np == 1) {
    ppath.r = r0;
    ppath.pos(0,0) = rte_pos[0];
    if (atmosphere_dim == 1)
      ppath.end_pos[1] = 0;
    ppath.los(0,joker) = rte_los[joker];
  }
  // Otherwise loop lvec (split in atm. dim. to make code more efficient)
  else {
    // 1D
    if (atmosphere_dim == 1) {
      ppath.end_pos[0] = rte_pos[0];
      ppath.end_pos[1] = 0;
      for (Index i=0; i<ppath.np; i++) {
        if (i > 0)
          ppath.lstep[i-1] = lvec[i]-lvec[i-1];;
        const Numeric l = lvec[i], x = x0+l*dx, y = y0+l*dy, z = z0+l*dz;
        ppath.r[i] = sqrt( x*x + y*y + z*z);;
        ppath.pos(i,0) = ppath.r[i] - refellipsoid[0];
        ppath.los(i,0) = geompath_za_at_r(ppath.constant,
                                          za0, // Will not work for limb sounding!
                                          ppath.r[i]); 
        ppath.pos(i,1) = geompath_lat_at_za(za0, 0, ppath.los(i,0));
      }
      gridpos(ppath.gp_p, z_field(joker,0,0), ppath.pos(joker,0));

    // 2D
    } else if (atmosphere_dim == 2) {
      ppath.end_pos[joker] = rte_pos[joker];
      const Index nz = z_field.npages(); 
      ArrayOfGridPos gp_z, gp_lat(1), gp_lon(0);
      gridpos_1to1(gp_z, Vector(nz));
      Tensor3 z_grid(nz,1,1);   // The altitudes at one (lat,lon)
      const Numeric lat1=lat_grid[0], lat2=last(lat_grid);
      ARTS_ASSERT( abs(dy) < 1e-9 );    // 2D happens strictly inside plane y=0
      for (Index i=0; i<ppath.np; i++) {
        if (i > 0)
          ppath.lstep[i-1] = lvec[i]-lvec[i-1];;
        const Numeric l = lvec[i], x = x0+l*dx, z = z0+l*dz;
        cart2poslos(ppath.r[i], ppath.pos(i,1), ppath.los(i,0),
                    x, z, dx, dz, ppath.constant, rte_pos[1], rte_los[0]);
        ARTS_USER_ERROR_IF (ppath.pos(i,1) < lat1,
            "The latitude grid must be extended downwards with at "
            "least ", lat1-ppath.pos(i,1), " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        ARTS_USER_ERROR_IF (ppath.pos(i,1) > lat2,
            "The latitude grid must be extended upwards with at "
            "least ", ppath.pos(i,1)-lat2, " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        gridpos(ppath.gp_lat[i], lat_grid, ppath.pos(i,1));
        //
        gridpos_copy(gp_lat[0], ppath.gp_lat[i]);
        regrid_atmfield_by_gp(z_grid, atmosphere_dim, z_field, gp_z, gp_lat, gp_lon);
        const Numeric r_e = refell2d(refellipsoid, lat_grid, ppath.gp_lat[i]);
        ppath.pos(i,0) = ppath.r[i] - r_e;
        gridpos(ppath.gp_p[i], z_grid(joker,0,0), ppath.pos(i,0));
      }

    // 3D
    } else {
      ppath.end_pos[joker] = rte_pos[joker];
      const Index nz = z_field.npages(); 
      ArrayOfGridPos gp_z, gp_lat(1), gp_lon(1);
      gridpos_1to1(gp_z, Vector(nz));
      Tensor3 z_grid(nz,1,1);   // The altitudes at one (lat,lon)
      const Numeric lat1=lat_grid[0], lat2=last(lat_grid);
      const Numeric lon1=lon_grid[0], lon2=last(lon_grid);
      for (Index i=0; i<ppath.np; i++) {
        if (i > 0)
          ppath.lstep[i-1] = lvec[i]-lvec[i-1];;
        const Numeric l = lvec[i], x = x0+l*dx, y = y0+l*dy, z = z0+l*dz;
        cart2poslos(ppath.r[i], ppath.pos(i,1), ppath.pos(i,2),
                    ppath.los(i,0), ppath.los(i,1), x, y, z, dx, dy, dz,
                    ppath.constant, x0, y0, z0, rte_pos[1], rte_pos[2],
                    rte_los[0], rte_los[1]);
        ARTS_USER_ERROR_IF (ppath.pos(i,1) < lat1,
            "The latitude grid must be extended downwards with at "
            "least ", lat1-ppath.pos(i,1), " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        ARTS_USER_ERROR_IF (ppath.pos(i,1) > lat2,
            "The latitude grid must be extended upwards with at "
            "least ", ppath.pos(i,1)-lat2, " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        ARTS_USER_ERROR_IF (ppath.pos(i,2) < lon1,
            "The latitude grid must be extended downwards with at "
            "least ", lon1-ppath.pos(i,2), " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        ARTS_USER_ERROR_IF (ppath.pos(i,2) > lon2,
            "The latitude grid must be extended upwards with at "
            "least ", ppath.pos(i,2)-lon2, " degrees to allow "
            "the ppath to fully be inside of the model atmosphere.")
        gridpos(ppath.gp_lat[i], lat_grid, ppath.pos(i,1));
        gridpos(ppath.gp_lon[i], lon_grid, ppath.pos(i,2));
        //
        gridpos_copy(gp_lat[0], ppath.gp_lat[i]);
        gridpos_copy(gp_lon[0], ppath.gp_lon[i]);
        regrid_atmfield_by_gp(z_grid, atmosphere_dim, z_field, gp_z, gp_lat, gp_lon);
        const Numeric r_e = refell2d(refellipsoid, lat_grid, ppath.gp_lat[i]);
        ppath.pos(i,0) = ppath.r[i] - r_e;
        gridpos(ppath.gp_p[i], z_grid(joker,0,0), ppath.pos(i,0));
      }
    }
  }
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

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathWriteXMLPartial(  //WS Input:
    const String& file_format,
    const Ppath& ppath,
    // WS Generic Input:
    const String& f,
    const Index& file_index,
    const Verbosity& verbosity) {
  String filename = f;
  Ppath ppath_partial = ppath;
  ArrayOfGridPos empty_gp;
  //Vector empty_v;

  ppath_partial.gp_p = empty_gp;
  ppath_partial.gp_lat = empty_gp;
  ppath_partial.gp_lon = empty_gp;
  //ppath_partial.nreal = empty_v;
  //ppath_partial.ngroup = empty_v;

  if (file_index >= 0) {
    // Create default filename if empty
    filename_xml_with_index(filename, file_index, "ppath");
  }

  WriteXML(file_format, ppath_partial, filename, 0, "ppath", "", "", verbosity);
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
void ppath_fieldCalc(Workspace& ws,
                     ArrayOfPpath& ppath_field,
                     const Agenda& ppath_agenda,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Index& atmgeom_checked,
                     const Vector& f_grid,
                     const Index& cloudbox_on,
                     const Index& cloudbox_checked,
                     const Index& ppath_inside_cloudbox_do,
                     const Matrix& sensor_pos,
                     const Matrix& sensor_los,
                     const Vector& rte_pos2,
                     const Verbosity& verbosity) {
  auto n = sensor_pos.nrows();
  ppath_field.resize(n);

  ARTS_USER_ERROR_IF (sensor_los.nrows() not_eq n,
        "Your sensor position matrix and sensor line of sight matrix do not match in size.\n");

  for (auto i = 0; i < n; i++)
    ppathCalc(ws,
              ppath_field[i],
              ppath_agenda,
              ppath_lmax,
              ppath_lraytrace,
              atmgeom_checked,
              f_grid,
              cloudbox_on,
              cloudbox_checked,
              ppath_inside_cloudbox_do,
              sensor_pos(i, joker),
              sensor_los(i, joker),
              rte_pos2,
              verbosity);
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
void rte_losReverse(
    Vector& rte_los,
    const Index& atmosphere_dim,
    const Verbosity&) {

  Vector los;
  Index l = rte_los.nelem();
  mirror_los(los, rte_los, atmosphere_dim);
  rte_los = los[Range(0,l)];
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(Vector& rte_los,
                const Index& atmosphere_dim,
                const Numeric& za,
                const Numeric& aa,
                const Verbosity&) {
  // Check input
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);

  if (atmosphere_dim == 1) {
    rte_los.resize(1);
  } else {
    rte_los.resize(2);
    rte_los[1] = aa;
  }
  rte_los[0] = za;
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

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posSet(Vector& rte_pos,
                const Index& atmosphere_dim,
                const Numeric& z,
                const Numeric& lat,
                const Numeric& lon,
                const Verbosity&) {
  // Check input
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);

  rte_pos.resize(atmosphere_dim);
  rte_pos[0] = z;
  if (atmosphere_dim >= 2) {
    rte_pos[1] = lat;
  }
  if (atmosphere_dim == 3) {
    rte_pos[2] = lon;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losBackwardToAltitude(Vector& rte_pos,
                                   Vector& rte_los,
                                   const Index& atmosphere_dim,
                                   const Vector& lat_grid,
                                   const Vector& lon_grid,
                                   const Vector& refellipsoid,
                                   const Numeric& altitude,
                                   const Index& los_is_reversed,
                                   const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(atmosphere_dim != 3, "This method only works for 3D.");

  // Find los to apply in next step
  Vector los2use;
  if (los_is_reversed) {
    los2use = rte_los;
  } else {
    mirror_los(los2use, rte_los, atmosphere_dim);
  }

  // Move in altitude
  Matrix start_pos(1,3), start_los(1,2), end_pos, end_los;
  start_pos(0, joker) = rte_pos;
  start_los(0, joker) = los2use;
  IntersectionGeometricalWithAltitude(end_pos,
                                      end_los,
                                      start_pos,
                                      start_los,
                                      refellipsoid,
                                      lat_grid,
                                      lon_grid,
                                      altitude,
                                      verbosity);

  // Extract final values
  rte_pos = end_pos(0, joker);
  mirror_los(rte_los, end_los(0, joker), atmosphere_dim);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losForwardToAltitude(Vector& rte_pos,
                                   Vector& rte_los,
                                   const Index& atmosphere_dim,
                                   const Vector& lat_grid,
                                   const Vector& lon_grid,
                                   const Vector& refellipsoid,
                                   const Numeric& altitude,
                                   const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(atmosphere_dim != 3, "This method only works for 3D.");

  // Move in altitude
  Matrix start_pos(1,3), start_los(1,2), end_pos, end_los;
  start_pos(0, joker) = rte_pos;
  start_los(0, joker) = rte_los;
  IntersectionGeometricalWithAltitude(end_pos,
                                      end_los,
                                      start_pos,
                                      start_los,
                                      refellipsoid,
                                      lat_grid,
                                      lon_grid,
                                      altitude,
                                      verbosity);

  // Extract final values
  rte_pos = end_pos(0, joker);
  rte_los = end_los(0, joker);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losStartOfPpath(Vector& rte_pos,
                             Vector& rte_los,
                             const Index& atmosphere_dim,
                             const Ppath& ppath,
                             const Verbosity&) {
  const Index np = ppath.np;

  // Check input
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  ARTS_USER_ERROR_IF (np == 0, "The input *ppath* is empty.");
  ARTS_USER_ERROR_IF (ppath.pos.nrows() != np,
        "Internal inconsistency in *ppath* (size of data "
        "does not match np).");

  rte_pos = ppath.pos(np - 1, Range(0, atmosphere_dim));
  if (atmosphere_dim < 3) {
    rte_los = ppath.los(np - 1, Range(0, 1));
  } else {
    rte_los = ppath.los(np - 1, Range(0, 2));
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricFromSensorPosToOtherPositions(
    Matrix& sensor_los,
    const Index& atmosphere_dim,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Vector& refellipsoid,
    const Matrix& sensor_pos,
    const Matrix& target_pos,
    const Verbosity& verbosity) {
  const Index n = sensor_pos.nrows();

  ARTS_USER_ERROR_IF (sensor_pos.ncols() != atmosphere_dim,
        "The number of columns of sensor_pos must be "
        "equal to the atmospheric dimensionality.");
  ARTS_USER_ERROR_IF ((atmosphere_dim == 1 && target_pos.ncols() != 2) ||
      (atmosphere_dim >= 2 && target_pos.ncols() != atmosphere_dim),
        "The number of columns of targe_pos must be equal to "
        "the atmospheric dimensionality, except for 1D where "
        "two columns are demended (as for *rte_pos2*).");
  ARTS_USER_ERROR_IF (target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same "
        "number of rows.");

  atmosphere_dim < 3 ? sensor_los.resize(n, 1) : sensor_los.resize(n, 2);
  Vector rte_los;
  for (Index i = 0; i < n; i++) {
    rte_losGeometricFromRtePosToRtePos2(rte_los,
                                        atmosphere_dim,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        sensor_pos(i, joker),
                                        target_pos(i, joker),
                                        verbosity);
    sensor_los(i, joker) = rte_los;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losReverse(
    Matrix& sensor_los,
    const Index& atmosphere_dim,
    const Verbosity&) {

  Vector los;
  Index l = sensor_los.ncols();
  for (Index i = 0; i < sensor_los.nrows(); i++) {
    mirror_los(los, sensor_los(i, joker), atmosphere_dim);
    sensor_los(i, joker) = los[Range(0,l)];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_pos_losBackwardToAltitude(Matrix& sensor_pos,
                                      Matrix& sensor_los,
                                      const Index& atmosphere_dim,
                                      const Vector& lat_grid,
                                      const Vector& lon_grid,
                                      const Vector& refellipsoid,
                                      const Numeric& altitude,
                                      const Index& los_is_reversed,
                                      const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(atmosphere_dim != 3, "This method only works for 3D.");

  // Find los to apply in next step
  Matrix los2use = sensor_los;
  if (!los_is_reversed) {
    sensor_losReverse(los2use, atmosphere_dim, verbosity);
  }

  // Move in altitude
  Matrix end_pos, end_los;
  IntersectionGeometricalWithAltitude(end_pos,
                                      end_los,
                                      sensor_pos,
                                      los2use,
                                      refellipsoid,
                                      lat_grid,
                                      lon_grid,
                                      altitude,
                                      verbosity);

  // Extract final values
  sensor_pos = end_pos;
  sensor_los = end_los;
  sensor_losReverse(sensor_los, atmosphere_dim, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_pos_losForwardToAltitude(Matrix& sensor_pos,
                                     Matrix& sensor_los,
                                     const Index& atmosphere_dim,
                                     const Vector& lat_grid,
                                     const Vector& lon_grid,
                                     const Vector& refellipsoid,
                                     const Numeric& altitude,
                                     const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(atmosphere_dim != 3, "This method only works for 3D.");

  // Move in altitude
  Matrix end_pos, end_los;
  IntersectionGeometricalWithAltitude(end_pos,
                                      end_los,
                                      sensor_pos,
                                      sensor_los,
                                      refellipsoid,
                                      lat_grid,
                                      lon_grid,
                                      altitude,
                                      verbosity);

  // Extract final values
  sensor_pos = end_pos;
  sensor_los = end_los;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void TangentPointExtract(Vector& tan_pos,
                         const Ppath& ppath,
                         const Verbosity&) {
  Index it;
  find_tanpoint(it, ppath);

  tan_pos.resize(ppath.pos.ncols());

  if (it < 0) {
    tan_pos = std::numeric_limits<Numeric>::quiet_NaN();
  } else {
    tan_pos[0] = ppath.pos(it, 0);
    tan_pos[1] = ppath.pos(it, 1);
    if (ppath.pos.ncols() == 3) {
      tan_pos[2] = ppath.pos(it, 2);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void TangentPointPrint(const Ppath& ppath,
                       const Index& level,
                       const Verbosity& verbosity) {
  Index it;
  find_tanpoint(it, ppath);

  ostringstream os;

  if (it < 0) {
    os << "Lowest altitude found at the end of the propagation path.\n"
       << "This indicates that the tangent point is either above the\n"
       << "top-of-the-atmosphere or below the planet's surface.";
  } else {
    os << "Tangent point position:\n-----------------------\n"
       << "     z = " << ppath.pos(it, 0) / 1e3 << " km\n"
       << "   lat = " << ppath.pos(it, 1) << " deg";
    if (ppath.pos.ncols() == 3)
      os << "\n   lon: " << ppath.pos(it, 2) << " deg";
  }

  CREATE_OUTS;
  SWITCH_OUTPUT(level, os.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZaRefr1D(Workspace& ws,
                          Vector& za_vector,
                          const Agenda& refr_index_air_agenda,
                          const Matrix& sensor_pos,
                          const Vector& p_grid,
                          const Tensor3& t_field,
                          const Tensor3& z_field,
                          const Tensor4& vmr_field,
                          const Vector& refellipsoid,
                          const Index& atmosphere_dim,
                          const Vector& f_grid,
                          const Vector& ztan_vector,
                          const Verbosity&) {
  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
                      "The function can only be used for 1D atmospheres.");

  ARTS_USER_ERROR_IF (ztan_vector.nelem() != sensor_pos.nrows(),
    "The number of altitudes in true tangent altitude vector must\n"
    "match the number of positions in *sensor_pos*.")

  // Set za_vector's size equal to ztan_vector
  za_vector.resize(ztan_vector.nelem());

  // Define refraction variables
  Numeric refr_index_air, refr_index_air_group;

  // Calculate refractive index for the tangential altitudes
  for (Index i = 0; i < ztan_vector.nelem(); i++) {
    ARTS_USER_ERROR_IF (ztan_vector[i] > sensor_pos(i, 0),
        "Invalid observation geometry: sensor (at z=", sensor_pos(i, 0),
        "m) is located below the requested tangent altitude (tanh=",
        ztan_vector[i], "m)")
    
    get_refr_index_1d(ws,
                      refr_index_air,
                      refr_index_air_group,
                      refr_index_air_agenda,
                      p_grid,
                      refellipsoid[0],
                      z_field,
                      t_field,
                      vmr_field,
                      f_grid,
                      ztan_vector[i] + refellipsoid[0]);

    // Calculate zenith angle
    za_vector[i] = 180 - RAD2DEG * asin(refr_index_air *
                                        (refellipsoid[0] + ztan_vector[i]) /
                                        (refellipsoid[0] + sensor_pos(i, 0)));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZa1D(Vector& za_vector,
                      const Matrix& sensor_pos,
                      const Vector& refellipsoid,
                      const Index& atmosphere_dim,
                      const Vector& ztan_vector,
                      const Verbosity&) {
  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
                      "The function can only be used for 1D atmospheres.");

  const Index npos = sensor_pos.nrows();

  ARTS_USER_ERROR_IF (ztan_vector.nelem() != npos,
      "The number of altitudes in the geometric tangent altitude vector\n"
      "must match the number of positions in *sensor_pos*.")

  za_vector.resize(npos);

  for (Index i = 0; i < npos; i++) {
    ARTS_USER_ERROR_IF (ztan_vector[i] > sensor_pos(i, 0),
        "Invalid observation geometry: sensor (at z=", sensor_pos(i, 0),
        "m) is located below the requested tangent altitude (tanh=",
        ztan_vector[i], "m)")
    
    za_vector[i] = geompath_za_at_r(refellipsoid[0] + ztan_vector[i],
                                    100,
                                    refellipsoid[0] + sensor_pos(i, 0));
  }
}

