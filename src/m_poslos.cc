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

#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "ppath.h"
#include "variousZZZ.h"

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
                   const Verbosity&)
{
  chk_rte_los("ref_los", ref_los);
  chk_sensor_los("other_los", other_los);

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
               const Verbosity&)
{
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
                 const Verbosity&)
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
                                const Vector& refellipsoid,
                                const Vector& rte_pos,
                                const Vector& target_pos,
                                const Verbosity&)
{
    chk_rte_pos("rte_pos", rte_pos);
    chk_rte_pos("target_pos", target_pos);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    rte_los.resize(2);

    geodetic2ecef(ecef, rte_pos, refellipsoid);
    geodetic2ecef(ecef_target, target_pos, refellipsoid);
    ecef_vector_distance(decef, ecef, ecef_target);
    decef /= norm2(decef);
    ecef2geodetic_los(dummy, rte_los, ecef, decef, refellipsoid);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losRefractedToPosition(Workspace& ws,
                                Vector& rte_los,
                                Ppath& ppath,
                                const Agenda& refr_index_air_ZZZ_agenda,
                                const Numeric& ppath_lstep,
                                const Numeric& ppath_lraytrace,
                                const Vector& refellipsoid,
                                const GriddedField2& surface_elevation,
                                const Numeric& surface_search_accuracy,
                                const Vector& rte_pos,
                                const Vector& target_pos,
                                const Numeric& target_dl,
                                const String& algorithm,
                                const Index& max_iterations,
                                const Index& robust,
                                const Numeric& z_toa,
                                const Index& do_horizontal_gradients,
                                const Index& do_twosided_perturb,
                                const Verbosity&)
{
    chk_rte_pos("rte_pos", rte_pos);
    chk_rte_pos("target_pos", target_pos);

    if (algorithm == "basic") {
      refracted_link_basic(ws,
                           ppath,
                           refr_index_air_ZZZ_agenda,
                           ppath_lstep,
                           ppath_lraytrace,
                           refellipsoid,
                           surface_elevation,
                           surface_search_accuracy,
                           z_toa,
                           do_horizontal_gradients,
                           do_twosided_perturb,
                           rte_pos,
                           target_pos,
                           target_dl,
                           max_iterations,
                           robust);

    } else {
      ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
    }

    rte_los = ppath.start_los;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losReverse(Vector& rte_los,
                    const Verbosity&)
{
  reverse_los(rte_los, rte_los);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(Vector& rte_los,
                const Numeric& za,
                const Numeric& aa,
                const Verbosity&)
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
                const Numeric& lon,
                const Verbosity&)
{
  rte_pos.resize(3);
  rte_pos[0] = z;
  rte_pos[1] = lat;
  rte_pos[2] = lon;
  chk_rte_pos("rte_pos", rte_pos);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losBackwardToAltitude(Vector& rte_pos,
                                   Vector& rte_los,
                                   const Vector& refellipsoid,
                                   const Numeric& altitude,
                                   const Index& los_is_reversed,
                                   const Verbosity& verbosity)
{
  // Find los to apply in next step
  Vector los2use;
  if (los_is_reversed) {
    los2use = rte_los;
  } else {
    reverse_los(los2use, rte_los);
  }

  // Move in altitude
  Matrix start_pos(1,3), start_los(1,2), end_pos, end_los;
  start_pos(0, joker) = rte_pos;
  start_los(0, joker) = los2use;
  IntersectionGeometricAltitude(end_pos,
                                end_los,
                                start_pos,
                                start_los,
                                refellipsoid,
                                altitude,
                                verbosity);

  // Extract final values
  rte_pos = end_pos(0, joker);
  reverse_los(rte_los, end_los(0, joker));
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losForwardToAltitude(Vector& rte_pos,
                                   Vector& rte_los,
                                   const Vector& refellipsoid,
                                   const Numeric& altitude,
                                   const Verbosity& verbosity)
{
  // Move in altitude
  Matrix start_pos(1,3), start_los(1,2), end_pos, end_los;
  start_pos(0, joker) = rte_pos;
  start_los(0, joker) = rte_los;
  IntersectionGeometricAltitude(end_pos,
                                end_los,
                                start_pos,
                                start_los,
                                refellipsoid,
                                altitude,
                                verbosity);

  // Extract final values
  rte_pos = end_pos(0, joker);
  rte_los = end_los(0, joker);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_losEndOfPpath(Vector& rte_pos,
                           Vector& rte_los,
                           const Ppath& ppath,
                           const Verbosity&)
{
  const Index np = ppath.np;

  // Check input
  ARTS_USER_ERROR_IF(np == 0, "The input *ppath* is empty.");

  rte_pos = ppath.pos(np - 1, joker);
  rte_los = ppath.los(np - 1, joker);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losAddLosAndDlos(Matrix& sensor_los,
                             const Vector& ref_los,
                             const Matrix& dlos,
                             const Verbosity&)
{
  chk_rte_los("ref_los", ref_los);
  ARTS_USER_ERROR_IF (dlos.ncols() != 2, "*dlos* must have two columns.");

  const Index nlos = dlos.nrows();
  sensor_los.resize(nlos, 2);

  for (Index i = 0; i < nlos; i++)
    add_za_aa(sensor_los(i, 0),
              sensor_los(i, 1),
              ref_los[0],
              ref_los[1],
              dlos(i, 0),
              dlos(i, 1));
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPosition(Matrix& sensor_los,
                                   const Vector& refellipsoid,
                                   const Matrix& sensor_pos,
                                   const Vector& target_pos,
                                   const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_rte_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    geodetic2ecef(ecef_target, target_pos, refellipsoid);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), refellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, refellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPositions(Matrix& sensor_los,
                                    const Vector& refellipsoid,
                                    const Matrix& sensor_pos,
                                    const Matrix& target_pos,
                                    const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_sensor_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    ARTS_USER_ERROR_IF(target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same number of rows.");
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), refellipsoid);
      geodetic2ecef(ecef_target, target_pos(i, joker), refellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, refellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losRefractedToPosition(Workspace& ws,
                                   Matrix& sensor_los,
                                   const Agenda& refr_index_air_ZZZ_agenda,
                                   const Numeric& ppath_lstep,
                                   const Numeric& ppath_lraytrace,
                                   const Vector& refellipsoid,
                                   const GriddedField2& surface_elevation,
                                   const Numeric& surface_search_accuracy,
                                   const Matrix& sensor_pos,
                                   const Vector& target_pos,
                                   const Numeric& target_dl,
                                   const String& algorithm,
                                   const Index& max_iterations,
                                   const Index& robust,
                                   const Numeric& z_toa,
                                   const Index& do_horizontal_gradients,
                                   const Index& do_twosided_perturb,
                                   const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_rte_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    sensor_los.resize(n, 2);

    if (algorithm == "basic") {
      for (Index i=0; i<n; ++i) {
        Ppath ppath;
        refracted_link_basic(ws,
                             ppath,
                             refr_index_air_ZZZ_agenda,
                             ppath_lstep,
                             ppath_lraytrace,
                             refellipsoid,
                             surface_elevation,
                             surface_search_accuracy,
                             z_toa,
                             do_horizontal_gradients,
                             do_twosided_perturb,
                             sensor_pos(i, joker),
                             target_pos,
                             target_dl,
                             max_iterations,
                             robust);
        sensor_los(i, joker) = ppath.start_los;
      }

  } else {
    ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losRefractedToPositions(Workspace& ws,
                                    Matrix& sensor_los,
                                    const Agenda& refr_index_air_ZZZ_agenda,
                                    const Numeric& ppath_lstep,
                                    const Numeric& ppath_lraytrace,
                                    const Vector& refellipsoid,
                                    const GriddedField2& surface_elevation,
                                    const Numeric& surface_search_accuracy,
                                    const Matrix& sensor_pos,
                                    const Matrix& target_pos,
                                    const Numeric& target_dl,
                                    const String& algorithm,
                                    const Index& max_iterations,
                                    const Index& robust,
                                    const Numeric& z_toa,
                                    const Index& do_horizontal_gradients,
                                    const Index& do_twosided_perturb,
                                    const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_sensor_pos("target_pos", target_pos);

    const Index n = sensor_pos.nrows();
    ARTS_USER_ERROR_IF(target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same number of rows.");
    sensor_los.resize(n, 2);

    if (algorithm == "basic") {
      for (Index i=0; i<n; ++i) {
        Ppath ppath;
        refracted_link_basic(ws,
                             ppath,
                             refr_index_air_ZZZ_agenda,
                             ppath_lstep,
                             ppath_lraytrace,
                             refellipsoid,
                             surface_elevation,
                             surface_search_accuracy,
                             z_toa,
                             do_horizontal_gradients,
                             do_twosided_perturb,
                             sensor_pos(i, joker),
                             target_pos(i, joker),
                             target_dl,
                             max_iterations,
                             robust);
        sensor_los(i, joker) = ppath.start_los;
      }

  } else {
    ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losReverse(
    Matrix& sensor_los,
    const Verbosity&)
{
  for (Index i = 0; i < sensor_los.nrows(); i++)
    reverse_los(sensor_los(i, joker), sensor_los(i, joker));
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_pos_losBackwardToAltitude(Matrix& sensor_pos,
                                      Matrix& sensor_los,
                                      const Vector& refellipsoid,
                                      const Numeric& altitude,
                                      const Index& los_is_reversed,
                                      const Verbosity& verbosity)
{
  // Find los to apply in next step
  Matrix los2use = sensor_los;
  if (!los_is_reversed) {
    sensor_losReverse(los2use, verbosity);
  }

  // Move in altitude
  Matrix end_pos, end_los;
  IntersectionGeometricAltitude(end_pos,
                                end_los,
                                sensor_pos,
                                los2use,
                                refellipsoid,
                                altitude,
                                verbosity);

  // Extract final values
  sensor_pos = end_pos;
  sensor_los = end_los;
  sensor_losReverse(sensor_los, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_pos_losForwardToAltitude(Matrix& sensor_pos,
                                     Matrix& sensor_los,
                                     const Vector& refellipsoid,
                                     const Numeric& altitude,
                                     const Verbosity& verbosity)
{
  // Move in altitude
  Matrix end_pos, end_los;
  IntersectionGeometricAltitude(end_pos,
                                end_los,
                                sensor_pos,
                                sensor_los,
                                refellipsoid,
                                altitude,
                                verbosity);

  // Extract final values
  sensor_pos = end_pos;
  sensor_los = end_los;
}
