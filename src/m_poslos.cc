/* Copyright (C) 2023 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
    @file    m_poslos.cc
    @author  Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date    2023-01-14

    @brief   Workspace methods for setting and extracting positions (pos),
             line-of-sights (los) and relative los (dlos).
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <workspace.h>

#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "ppath.h"
#include "surf.h"
#include "variousZZZ.h"

inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric NAT_LOG_2=Constant::ln_2;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void dlosGauss(Matrix& dlos,
               Vector& dlos_weight_vector,
               const Numeric& fwhm_deg,
               const Index& npoints,
               const Index& include_response_in_weight) {
  // Use FWHM and sigma in radians to get solid angles right
  const Numeric fwhm = DEG2RAD *fwhm_deg;
  const Numeric si = fwhm / (2 * sqrt(2 * NAT_LOG_2));

  // The product between Gauss and radius gives the "density" for the
  // sampling. We place points in radius to cover the cumulative
  // distribution of normalised density with npoints bins. That is,
  // the two first points cover [0,1/npoints[ and [1/npoints,2/npoints[,
  // respectively. To get a correct geo-positioning, we still want a
  // point at (0,0) and the first point is at the end shifted to (0,0).
  
  // Cumulative distribution of Gauss weighted area, as a function of radius x
  Vector xp, cx;
  {
    linspace(xp, 0, 2.0*fwhm, 0.02*fwhm);
    Vector gx;
    VectorGaussian(gx, xp, 0, -1.0, fwhm);
    gx *= xp;  // Weight with radius, no need to include pi
    const Index np = gx.nelem();
    cx.resize(np);
    cumsum(cx, gx);
    cx /= cx[np-1];
  }
  // Flat distribution with bins covering [0, 1]
  Vector cp;
  {
    const Numeric halfbin = 0.5 / (Numeric) npoints;
    nlinspace(cp, halfbin, 1-halfbin, npoints);
  }
  
  // Radii of layers, obtained by interpolating xp(cx) to cp
  Vector r(npoints);
  {
    ArrayOfGridPos gp(npoints);
    gridpos(gp, cx, cp);
    Matrix itw(npoints, 2);
    interpweights(itw, gp);
    interp(r, itw, xp, gp);
  }

  // Factor to rescale Gauss(x) from 1D to 2D (along y=0)
  const Numeric scfac = 1 / (sqrt(2 * PI) * si);
  
  // Prepare for calculating dlos and weights
  dlos.resize(npoints, 2);
  //
  dlos_weight_vector.resize(npoints);
  dlos_weight_vector = 1.0 / (Numeric) npoints;
  // If include_response_in_weight, all weights equal.
  // Otherwise, we need to divide with value of 2D Gauss for radius:
  Vector gv(npoints);
  if (!include_response_in_weight) {
    VectorGaussian(gv, r, 0, -1.0, fwhm);
  }

  // Calculate
  const Numeric dangle = 0.58 * 2.0 * PI;
  for (Index i=0; i<npoints; ++i) {
    const Numeric angle = dangle * (Numeric) i;
    dlos(i, 0) = r[i] * cos(angle);
    dlos(i, 1) = r[i] * sin(angle);
    if (!include_response_in_weight) {
      dlos_weight_vector[i] = dlos_weight_vector[i] / (scfac * gv[i]);
    }
  }
  dlos(0, joker) = 0.0;
  dlos *= RAD2DEG;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dlosUniform(Matrix& dlos,
                 Vector& dlos_weight_vector,
                 const Numeric& width,
                 const Index& npoints,
                 const Index& crop_circular)
{
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
void rte_losGeometricToPosition(Vector& rte_los,
                                const SurfaceField& surface_field,
                                const Vector& rte_pos,
                                const Vector& target_pos)
{
    chk_rte_pos("rte_pos", rte_pos);
    chk_rte_pos("target_pos", target_pos);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    rte_los.resize(2);

    geodetic2ecef(ecef, rte_pos, surface_field.ellipsoid);
    geodetic2ecef(ecef_target, target_pos, surface_field.ellipsoid);
    ecef_vector_distance(decef, ecef, ecef_target);
    decef /= norm2(decef);
    ecef2geodetic_los(dummy, rte_los, ecef, decef, surface_field.ellipsoid);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losReverse(Vector& rte_los)
{
  reverse_los(rte_los, rte_los);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(Vector& rte_los,
                const Numeric& za,
                const Numeric& aa)
{
  rte_los.resize(2);
  rte_los[0] = za;
  rte_los[1] = aa;
  chk_rte_los("rte_los", rte_los);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posSet(Vector& rte_pos,
                const Numeric& z,
                const Numeric& lat,
                const Numeric& lon)
{
  rte_pos.resize(3);
  rte_pos[0] = z;
  rte_pos[1] = lat;
  rte_pos[2] = lon;
  chk_rte_pos("rte_pos", rte_pos);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losEndOfPpath(Vector& rte_pos,
                           Vector& rte_los,
                           const Ppath& ppath)
{
  const Index np = ppath.np;

  // Check input
  ARTS_USER_ERROR_IF(np == 0, "The input *ppath* is empty.");

  rte_pos = ppath.pos(np - 1, joker);
  rte_los = ppath.los(np - 1, joker);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPosition(Matrix& sensor_los,
                                   const SurfaceField& surface_field,
                                   const Matrix& sensor_pos,
                                   const Vector& target_pos)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_rte_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    geodetic2ecef(ecef_target, target_pos, surface_field.ellipsoid);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), surface_field.ellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, surface_field.ellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPositions(Matrix& sensor_los,
                                    const SurfaceField& surface_field,
                                    const Matrix& sensor_pos,
                                    const Matrix& target_pos)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_sensor_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    ARTS_USER_ERROR_IF(target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same number of rows.");
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), surface_field.ellipsoid);
      geodetic2ecef(ecef_target, target_pos(i, joker), surface_field.ellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, surface_field.ellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losReverse(
    Matrix& sensor_los)
{
  for (Index i = 0; i < sensor_los.nrows(); i++)
    reverse_los(sensor_los(i, joker), sensor_los(i, joker));
}
