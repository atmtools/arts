/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/**
   @file   surface.cc
   @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   @date   2012-02-06 

   This file contains internal functions associated with the surface.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "surface.h"
#include <cmath>
#include "auto_md.h"
#include "check_input.h"
#include "matpack_complex.h"
#include "geodetic.h"
#include "geodetic_OLD.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "variousZZZ.h"
#include "workspace_ng.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

// Expression double-checked 210330 (PE)
Numeric calc_incang(ConstVectorView rte_los, ConstVectorView specular_los) {
  return (180 - abs(rte_los[0]) + abs(specular_los[0])) / 2;
}

Index index_of_zsurface(const Numeric& z_surface,
                        ConstVectorView z_profile) {
  Index ip = 0;
  while (z_surface >= z_profile[ip+1]) {
    ip++;
  }
  return ip;
}

void surface_calc(Matrix& iy,
                  ConstTensor3View I,
                  ConstMatrixView surface_los,
                  ConstTensor4View surface_rmatrix,
                  ConstMatrixView surface_emission) {
  // Some sizes
  const Index nf = I.nrows();
  const Index stokes_dim = I.ncols();
  const Index nlos = surface_los.nrows();

  iy = surface_emission;

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector rtmp(stokes_dim);  // Reflected Stokes vector for 1 frequency

      for (Index iv = 0; iv < nf; iv++) {
        mult(rtmp, surface_rmatrix(ilos, iv, joker, joker), I(ilos, iv, joker));
        iy(iv, joker) += rtmp;
      }
    }
  }
}

void surface_specular_R_and_b(MatrixView surface_rmatrix,
                              VectorView surface_emission,
                              const Complex& Rv,
                              const Complex& Rh,
                              const Numeric& f,
                              const Index& stokes_dim,
                              const Numeric& surface_skin_t) {
  ARTS_ASSERT(surface_rmatrix.nrows() == stokes_dim);
  ARTS_ASSERT(surface_rmatrix.ncols() == stokes_dim);
  ARTS_ASSERT(surface_emission.nelem() == stokes_dim);

  // Expressions are derived in the surface chapter in the user guide

  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  Numeric B = planck(f, surface_skin_t);

  const Numeric rv = pow(abs(Rv), 2.0);
  const Numeric rh = pow(abs(Rh), 2.0);
  const Numeric rmean = (rv + rh) / 2;

  surface_rmatrix(0, 0) = rmean;
  surface_emission[0] = B * (1 - rmean);

  if (stokes_dim > 1) {
    const Numeric rdiff = (rv - rh) / 2;

    surface_rmatrix(1, 0) = rdiff;
    surface_rmatrix(0, 1) = rdiff;
    surface_rmatrix(1, 1) = rmean;
    surface_emission[1] = -B * rdiff;

    if (stokes_dim > 2) {
      const Complex a = Rh * conj(Rv);
      const Complex b = Rv * conj(Rh);
      const Numeric c = real(a + b) / 2.0;

      surface_rmatrix(2, 2) = c;

      if (stokes_dim > 3) {
        const Numeric d = imag(a - b) / 2.0;

        surface_rmatrix(2, 3) = d;
        surface_rmatrix(3, 2) = -d;
        surface_rmatrix(3, 3) = c;
      }
    }
  }
}

void surface_props_check(const Index& atmosphere_dim,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& surface_props_data,
                         const ArrayOfString& surface_props_names) {
  // Check sizes
  ARTS_USER_ERROR_IF (surface_props_data.npages() != surface_props_names.nelem(),
        "The number of pages in *surface_props_data* and "
        "length of *surface_props_names* differ.");
  // If no surface properties, then we are ready
  if (surface_props_names.nelem() == 0) {
    return;
  }
  ARTS_USER_ERROR_IF (surface_props_data.nrows() !=
      (atmosphere_dim == 1 ? 1 : lat_grid.nelem()),
                      "Row-size of *surface_props_data* not as expected.");
  ARTS_USER_ERROR_IF (surface_props_data.ncols() !=
      (atmosphere_dim <= 2 ? 1 : lon_grid.nelem()),
                      "Column-size of *surface_props_data* not as expected.");

  for (Index i = 0; i < surface_props_names.nelem(); i++) {
    ARTS_USER_ERROR_IF (surface_props_names[i].nelem() == 0,
      "Element ", i, " (0-based) of *surface_props_names* is empty.")
    for (Index j = i + 1; j < surface_props_names.nelem(); j++) {
      ARTS_USER_ERROR_IF (surface_props_names[j] == surface_props_names[i],
        "Two surface properties with same name found!\n"
        "This found for these two properties\n"
        "   index: ", i, '\n',
        "   index: ", j, '\n',
        "    name: ", surface_props_names[i])
    }
  }
}

void surface_props_interp(Vector& v,
                          const String& vname,
                          const Index& atmosphere_dim,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const Matrix& itw,
                          const Tensor3& surface_props_data,
                          const ArrayOfString& surface_props_names) {
  ARTS_ASSERT(v.nelem() == 1);
  ARTS_ASSERT(surface_props_data.npages() == surface_props_names.nelem());

  for (Index i = 0; i < surface_props_names.nelem(); i++) {
    if (surface_props_names[i] == vname) {
      interp_atmsurface_by_itw(v,
                               atmosphere_dim,
                               surface_props_data(i, joker, joker),
                               gp_lat,
                               gp_lon,
                               itw);
      return;
    }
  }

  ARTS_USER_ERROR (
                      "The following property was requested\n"
                      "   ", vname, '\n',
                      "but it could not be found in *surface_props_names*.")
}

void dsurface_check(const ArrayOfString& surface_props_names,
                    const ArrayOfString& dsurface_names,
                    const ArrayOfTensor4 dsurface_rmatrix_dx,
                    const ArrayOfMatrix& dsurface_emission_dx) {
  const Index nq = dsurface_names.nelem();

  ARTS_USER_ERROR_IF (dsurface_rmatrix_dx.nelem() != nq,
        "The lengths of *dsurface_names* and *dsurface_rmatrix_dx* differ.");
  ARTS_USER_ERROR_IF (dsurface_emission_dx.nelem() != nq,
        "The lengths of *dsurface_names* and *dsurface_emission_dx* differ.");

  for (Index i = 0; i < nq; i++) {
    bool found = false;
    for (Index j = 0; j < surface_props_names.nelem() && !found; j++) {
      if (dsurface_names[i] == surface_props_names[j]) {
        found = true;
      }
    }
    ARTS_USER_ERROR_IF (!found,
        "String ", i, " (0-based) of *dsurface_names* is \"",
        dsurface_names[i], "\"\n"
        "but this string could not be found in *surface_props_names*.\n"
        "This is likely due to incorrect choice of quantity when\n"
        " calling *jacobianAddSurfaceQuantity*.")
  }
}


void surface_get_incoming_direct(
    Workspace& ws,
    Matrix& iy_incoming,
    Index& stars_visible,
    Vector& specular_los,
    const Vector& rtp_pos,
    const Vector& rtp_los,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    const Tensor3& t_field,
    const EnergyLevelMap& nlte_field,
    const Tensor4& vmr_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Tensor3& wind_u_field,
    const Tensor3& wind_v_field,
    const Tensor3& wind_w_field,
    const Tensor3& mag_u_field,
    const Tensor3& mag_v_field,
    const Tensor3& mag_w_field,
    const Matrix& z_surface,
    const Vector& refellipsoid,
    const Tensor4& pnd_field,
    const ArrayOfTensor4& dpnd_field_dx,
    const ArrayOfString& scat_species,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Index& ppath_inside_cloudbox_do,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Index& gas_scattering_do,
    const Index& jacobian_do,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfSun& suns,
    const Numeric& rte_alonglos_v,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const Agenda& gas_scattering_agenda,
    const Agenda& ppath_step_agenda,
    const Verbosity& verbosity){

  //Allocate
  Vector surface_normal;
  Matrix iy_sun_toa;

  //get specular line of sight
  specular_losCalcOld(specular_los,
                   surface_normal,
                   rtp_pos,
                   rtp_los,
                   atmosphere_dim,
                   lat_grid,
                   lon_grid,
                   refellipsoid,
                   z_surface,
                   0,
                   verbosity);

  //calculate propagation path from the surface to the space in line of sight
  Ppath ppath;
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
             rtp_pos,
             specular_los,
             ppath_lmax,
             ppath_lraytrace,
             ppath_inside_cloudbox_do,
             verbosity);


  //get the incoming spectral radiance of the sun at toa. If there is no in
  //line of sight, then iy_sun_toa is simply zero and we are finished. No further
  //calculations needed.
  stars_visible=0;
  get_sun_background(iy_sun_toa,
                      stars_visible,
                      suns,
                      ppath,
                      f_grid,
                      stokes_dim,
                      atmosphere_dim,
                      refellipsoid);

  if (stars_visible){

    //dummy variables needed for the output and input of iyTransmission
    ArrayOfMatrix iy_aux_dummy;
    ArrayOfString iy_aux_vars_dummy;
    Vector ppvar_p_dummy;
    Vector ppvar_t_dummy;
    EnergyLevelMap ppvar_nlte_dummy;
    Matrix ppvar_vmr_dummy;
    Matrix ppvar_wind_dummy;
    Matrix ppvar_mag_dummy;
    Matrix ppvar_pnd_dummy;
    Matrix ppvar_f_dummy;
    Tensor3 ppvar_iy_dummy;
    Tensor4 ppvar_trans_cumulat_dummy;
    Tensor4 ppvar_trans_partial_dummy;
    ArrayOfTensor3 diy_incoming_dummy;

    //Calculate the transmitted radiation from toa to the surface
    iyTransmissionStandard(ws,
                           iy_incoming,
                           iy_aux_dummy,
                           diy_incoming_dummy,
                           ppvar_p_dummy,
                           ppvar_t_dummy,
                           ppvar_nlte_dummy,
                           ppvar_vmr_dummy,
                           ppvar_wind_dummy,
                           ppvar_mag_dummy,
                           ppvar_pnd_dummy,
                           ppvar_f_dummy,
                           ppvar_iy_dummy,
                           ppvar_trans_cumulat_dummy,
                           ppvar_trans_partial_dummy,
                           stokes_dim,
                           f_grid,
                           atmosphere_dim,
                           p_grid,
                           t_field,
                           nlte_field,
                           vmr_field,
                           abs_species,
                           wind_u_field,
                           wind_v_field,
                           wind_w_field,
                           mag_u_field,
                           mag_v_field,
                           mag_w_field,
                           cloudbox_on,
                           cloudbox_limits,
                           gas_scattering_do,
                           pnd_field,
                           dpnd_field_dx,
                           scat_species,
                           scat_data,
                           iy_aux_vars_dummy,
                           jacobian_do,
                           jacobian_quantities,
                           ppath,
                           iy_sun_toa,
                           propmat_clearsky_agenda,
                           water_p_eq_agenda,
                           gas_scattering_agenda,
                           1,
                           Tensor3(),
                           rte_alonglos_v,
                           verbosity);

  } else {
    iy_incoming.resize(f_grid.nelem(),stokes_dim);
    iy_incoming=0;
  }
}


void surface_normal_calc(VectorView pos,
                         VectorView ecef,
                         VectorView decef,
                         const Vector& refellipsoid,
                         const GriddedField2& surface_elevation,
                         ConstVectorView pos2D)
{
  ARTS_ASSERT(pos.nelem() == 3); 
  ARTS_ASSERT(ecef.nelem() == 3); 
  ARTS_ASSERT(decef.nelem() == 3); 
  ARTS_ASSERT(pos2D.nelem() == 2);
  
  // We need two orthogonal vectors inside the surface plane. We
  // follow ENU and the first one should be towards E and the second
  // towards N. We obtain the vectors quite easily by dl shifts,
  // except when we are are very close to the North or South pole. To
  // be sure that a dl shift does not pass any of the poles, we
  // consider here that all positions inside a distance 5*dl on the
  // side. These points are shifted to be at the pole.
  //
  const Numeric dl = 1.0;  
  const Numeric lat_limit = 90.0 - RAD2DEG * 5 * dl / refellipsoid[1];

  // Determine 3D pos at pos
  Numeric lat = pos2D[0];
  if (lat > lat_limit)
    lat = 90.0;
  else if (lat < -lat_limit)
    lat = -90.0;
  //
  pos[1] = lat;
  pos[2] = pos2D[1];
  pos[0] = interp_gfield2(surface_elevation, pos[Range(1, 2)]);

  // Radius at pos0
  const Numeric r = pos[0] + prime_vertical_radius(refellipsoid, lat);

  // Shifted positions
  Vector posWE = pos, posSN = pos;
  //
  // North pole case
  if (lat > lat_limit) {
    posSN[1] -= RAD2DEG * dl / r;  
    posSN[2] = 90; 
    posWE[1] = posSN[1];
    posWE[2] = 0; 
  // South pole case
  } else if (lat < -lat_limit) {
    posSN[1] += RAD2DEG * dl / r;  
    posSN[2] = 0; 
    posWE[1] = posSN[1];
    posWE[2] = 90; 
  // A general case 
  } else {
    posSN[1] += RAD2DEG * dl / r;  
    posWE[2] += RAD2DEG * dl / (r * cos(DEG2RAD * posWE[1]));
    if (posWE[2] >= 180)
      posWE[2] -= 360;
  }
  //
  posSN[0] = interp_gfield2(surface_elevation, posSN[Range(1, 2)]);
  posWE[0] = interp_gfield2(surface_elevation, posWE[Range(1, 2)]);
  
  // Convert all three positions to ECEF
  Vector ecefSN(3), ecefWE(3);
  geodetic2ecef(ecef, pos, refellipsoid);
  geodetic2ecef(ecefSN, posSN, refellipsoid);
  geodetic2ecef(ecefWE, posWE, refellipsoid);

  // Directional vectors to shifted positions
  Vector decefSN(3), decefWE(3);
  ecef_vector_distance(decefSN, ecef, ecefSN);
  ecef_vector_distance(decefWE, ecef, ecefWE);

  // Normal is cross product of the two decef (
  cross3(decef, decefWE, decefSN);
  decef /= norm2(decef);
}
