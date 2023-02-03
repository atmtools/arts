/* Copyright (C) 2002-2012 Patrick Eriksson <patrick.eriksson@chalmers.se>

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

/*!
  \file   refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-17
  
  \brief  Functions releated to calculation of refractive index.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "refraction.h"
#include <cmath>
#include "arts_constants.h"
#include "auto_md.h"
#include "matpack_complex.h"
#include "geodetic.h"
#include "interpolation.h"
#include "special_interp.h"
#include "check_input.h"
#include "arts_conversions.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! complex_n_water_liebe93
/*! 
  Complex refractive index of liquid water according to Liebe 1993.

  The method treats liquid water without salt. Thus, not valid below 10 GHz.
  Upper frequency limit not known, here set to 1000 GHz. Model parameters taken
  from Atmlab function epswater93 (by C. Maetzler), which refer to Liebe 1993
  without closer specifications.
 
  Temperature must be between 0 and 100 degrees Celsius.

  The output matrix has two columns, where column 0 is real part and column 1
  is imaginary part. And rows matches f_grid.
   
   \param   complex_n   Out: Complex refractive index.        
   \param   f_grid      As the WSV with the same name.
   \param   t           Temperature

   \author Patrick Eriksson
   \date   2003-08-15
*/
void complex_n_water_liebe93(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t) {
  chk_if_in_range("t", t, TEMP_0_C - 40, TEMP_0_C + 100);
  chk_if_in_range("min of f_grid", min(f_grid), 10e9, 1000e9);
  chk_if_in_range("max of f_grid", max(f_grid), 10e9, 1000e9);

  const Index nf = f_grid.nelem();

  complex_n.resize(nf, 2);

  // Implementation following epswater93.m (by C. Mätzler), part of Atmlab,
  // but numeric values strictly following the paper version (146, not 146.4)
  const Numeric theta = 1 - 300 / t;
  const Numeric e0 = 77.66 - 103.3 * theta;
  const Numeric e1 = 0.0671 * e0;
  const Numeric f1 = 20.2 + 146 * theta + 316 * theta * theta;
  const Numeric e2 = 3.52;
  const Numeric f2 = 39.8 * f1;

  for (Index iv = 0; iv < nf; iv++) {
    const Complex ifGHz(0.0, f_grid[iv] / 1e9);

    Complex n = sqrt(e2 + (e1 - e2) / (Numeric(1.0) - ifGHz / f2) +
                     (e0 - e1) / (Numeric(1.0) - ifGHz / f1));

    complex_n(iv, 0) = n.real();
    complex_n(iv, 1) = n.imag();
  }
}

//! complex_n_ice_matzler06
/*! 
  Complex refractive index of water ice according to Matzler 2006 (equivalent to
  Warren 2008).

  The method treats pure water ice (no impurities like salt). Valid from 10 MHz
  up to 3 THz. Thus, not valid below 10 GHz. Follows the atmlab implementation,
  including some relaxation of upper temperature limit to 280K.

  The output matrix has two columns, where column 0 is real part and column 1
  is imaginary part; rows match f_grid.
   
   \param   complex_n   Out: Complex refractive index.        
   \param   f_grid      As the WSV with the same name.
   \param   t           Temperature

   \author Jana Mendrok
   \date   2016-03-21
*/
void complex_n_ice_matzler06(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t) {
  chk_if_in_range("t", t, 20., 280.);
  chk_if_in_range("min of f_grid", min(f_grid), 10e6, 3000e9);
  chk_if_in_range("max of f_grid", max(f_grid), 10e6, 3000e9);

  const Index nf = f_grid.nelem();

  complex_n.resize(nf, 2);

  // some parametrization constants
  const Numeric B1 = 0.0207;
  const Numeric B2 = 1.16e-11;
  const Numeric b = 335.;

  const Numeric deltabeta = exp(-9.963 + 0.0372 * (t - 273));
  const Numeric ebdt = exp(b / t);
  const Numeric betam = (B1 / t) * ebdt / ((ebdt - 1.) * (ebdt - 1.));

  const Numeric theta = 300. / t - 1;
  const Numeric alfa = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);
  const Numeric reps = 3.1884 + 9.1e-4 * (t - 273);

  for (Index iv = 0; iv < nf; iv++) {
    Numeric f = f_grid[iv] / 1e9;
    Numeric beta = betam + B2 * f * f + deltabeta;
    Numeric ieps = alfa / f + beta * f;

    Complex eps(reps, ieps);
    Complex n = sqrt(eps);
    complex_n(iv, 0) = n.real();
    complex_n(iv, 1) = n.imag();
  }
}

//! get_refr_index_1d
/*! 
   Extracts the refractive index for 1D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_air_agenda* to determine the refractive index for the
   given point.

   The atmosphere is given by its 1D view. That is, the latitude and
   longitude dimensions are removed from the atmospheric fields. For
   example, the temperature is given as a vector (the vertical profile).

   \param   ws                  Current Workspace
   \param   refr_index_air          Output: As the WSV with the same name.
   \param   refr_index_air_group    Output: As the WSV with the same name.
   \param   refr_index_air_agenda   As the WSV with the same name.
   \param   p_grid              As the WSV with the same name.   
   \param   refellipsoid        As the WSV with the same name.
   \param   z_field             As the WSV with the same name.
   \param   t_field             As the WSV with the same name.
   \param   vmr_field           As the WSV with the same name.
   \param   f_grid              As the WSV with the same name.
   \param   r                   The radius of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-16
*/
void get_refr_index_1d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r) {
  Numeric rtp_pressure, rtp_temperature;
  Vector rtp_vmr;

  // Pressure grid position
  ArrayOfGridPos gp(1);
  gridpos(gp, z_field(joker, 0, 0), Vector(1, r - refellipsoid[0]));

  // Altitude interpolation weights
  Matrix itw(1, 2);
  interpweights(itw, gp);

  // Pressure
  Vector dummy(1);
  itw2p(dummy, p_grid, gp, itw);
  rtp_pressure = dummy[0];

  // Temperature
  interp(dummy, itw, t_field(joker, 0, 0), gp);
  rtp_temperature = dummy[0];

  // VMR
  const Index ns = vmr_field.nbooks();
  //
  rtp_vmr.resize(ns);
  //
  for (Index is = 0; is < ns; is++) {
    interp(dummy, itw, vmr_field(is, joker, 0, 0), gp);
    rtp_vmr[is] = dummy[0];
  }

  refr_index_air_agendaExecute(ws,
                               refr_index_air,
                               refr_index_air_group,
                               rtp_pressure,
                               rtp_temperature,
                               rtp_vmr,
                               Vector{f_grid},
                               refr_index_air_agenda);
}

//! get_refr_index_2d
/*! 
   Extracts the refractive index for 2D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_air_agenda* to determine the refractive index for the
   given point.

   The atmosphere is given by its 2D view. That is, the longitude
   dimension is removed from the atmospheric fields. For example,
   the temperature is given as a matrix.

   \param   ws                      Current Workspace
   \param   refr_index_air          Output: As the WSV with the same name.
   \param   refr_index_air_group    Output: As the WSV with the same name.
   \param   refr_index_air_agenda   As the WSV with the same name.
   \param   p_grid                  As the WSV with the same name.
   \param   lat_grid                As the WSV with the same name.
   \param   refellipsoid            As the WSV with the same name.
   \param   z_field                 As the WSV with the same name.
   \param   t_field                 As the WSV with the same name.
   \param   vmr_field               As the WSV with the same name.
   \param   f_grid                  As the WSV with the same name.
   \param   r                       The radius of the position of interest.
   \param   lat                     The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-14
*/
void get_refr_index_2d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat) {
  Numeric rtp_pressure, rtp_temperature;
  Vector rtp_vmr;

  // Determine the geometric altitudes at *lat*
  const Index np = p_grid.nelem();
  Vector z_grid(np);
  ArrayOfGridPos gp_lat(1);
  //
  gridpos(gp_lat, lat_grid, ExhaustiveConstVectorView{lat});
  z_at_lat_2d(z_grid, p_grid, lat_grid, z_field(joker, joker, 0), gp_lat[0]);

  // Determine the ellipsoid radius at *lat*
  const Numeric rellips = refell2d(refellipsoid, lat_grid, gp_lat[0]);

  // Altitude (equal to pressure) grid position
  ArrayOfGridPos gp_p(1);
  gridpos(gp_p, z_grid, Vector(1, r - rellips));

  // Altitude interpolation weights
  Matrix itw(1, 2);
  Vector dummy(1);
  interpweights(itw, gp_p);

  // Pressure
  itw2p(dummy, p_grid, gp_p, itw);
  rtp_pressure = dummy[0];

  // Temperature
  itw.resize(1, 4);
  interpweights(itw, gp_p, gp_lat);
  interp(dummy, itw, t_field(joker, joker, 0), gp_p, gp_lat);
  rtp_temperature = dummy[0];

  // VMR
  const Index ns = vmr_field.nbooks();
  //
  rtp_vmr.resize(ns);
  //
  for (Index is = 0; is < ns; is++) {
    interp(dummy, itw, vmr_field(is, joker, joker, 0), gp_p, gp_lat);
    rtp_vmr[is] = dummy[0];
  }

  refr_index_air_agendaExecute(ws,
                               refr_index_air,
                               refr_index_air_group,
                               rtp_pressure,
                               rtp_temperature,
                               rtp_vmr,
                               Vector{f_grid},
                               refr_index_air_agenda);
}

/*! get_refr_index_3d

   Extracts the refractive index for 3D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_air_agenda* to determine the refractive index for the
   given point.

   \param   ws                      Current Workspace
   \param   refr_index_air          Output: As the WSV with the same name.
   \param   refr_index_air_group    Output: As the WSV with the same name.
   \param   refr_index_air_agenda   As the WSV with the same name.
   \param   p_grid                  As the WSV with the same name.
   \param   lat_grid                As the WSV with the same name.
   \param   lon_grid                As the WSV with the same name.
   \param   refellipsoid            As the WSV with the same name.
   \param   z_field                 As the WSV with the same name.
   \param   t_field                 As the WSV with the same name.
   \param   vmr_field               As the WSV with the same name.
   \param   f_grid                  As the WSV with the same name.
   \param   r                       The radius of the position of interest.
   \param   lat                     The latitude of the position of interest.
   \param   lon                     The longitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-17
*/
void get_refr_index_3d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView lon_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat,
                       const Numeric& lon) {
  Numeric rtp_pressure, rtp_temperature;
  Vector rtp_vmr;

  // Determine the geometric altitudes at *lat* and *lon*
  const Index np = p_grid.nelem();
  Vector z_grid(np);
  ArrayOfGridPos gp_lat(1), gp_lon(1);
  //
  gridpos(gp_lat, lat_grid, ExhaustiveConstVectorView{lat});
  gridpos(gp_lon, lon_grid, ExhaustiveConstVectorView{lon});
  z_at_latlon(
      z_grid, p_grid, lat_grid, lon_grid, z_field, gp_lat[0], gp_lon[0]);

  // Determine the elipsoid radius at *lat*
  const Numeric rellips = refell2d(refellipsoid, lat_grid, gp_lat[0]);

  // Altitude (equal to pressure) grid position
  ArrayOfGridPos gp_p(1);
  gridpos(gp_p, z_grid, Vector(1, r - rellips));

  // Altitude interpolation weights
  Matrix itw(1, 2);
  Vector dummy(1);
  interpweights(itw, gp_p);

  // Pressure
  itw2p(dummy, p_grid, gp_p, itw);
  rtp_pressure = dummy[0];

  // Temperature
  itw.resize(1, 8);
  interpweights(itw, gp_p, gp_lat, gp_lon);
  interp(dummy, itw, t_field, gp_p, gp_lat, gp_lon);
  rtp_temperature = dummy[0];

  // VMR
  const Index ns = vmr_field.nbooks();
  //
  rtp_vmr.resize(ns);
  //
  for (Index is = 0; is < ns; is++) {
    interp(
        dummy, itw, vmr_field(is, joker, joker, joker), gp_p, gp_lat, gp_lon);
    rtp_vmr[is] = dummy[0];
  }

  refr_index_air_agendaExecute(ws,
                               refr_index_air,
                               refr_index_air_group,
                               rtp_pressure,
                               rtp_temperature,
                               rtp_vmr,
                               Vector{f_grid},
                               refr_index_air_agenda);
}

//! refr_gradients_1d
/*! 
   Determines the refractive index, and its gradients, for the given position.

   The gradients are calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius or
   latitude and the difference to the refractive index at the given
   point determines the gradient.

   \param   ws                    Current Workspace
   \param   refr_index_air        Output: As the WSV with the same name.
   \param   refr_index_air_group  Output: As the WSV with the same name.
   \param   dndr                  Output: Radial gradient of refractive index.
   \param   refr_index_air_agenda As the WSV with the same name.
   \param   p_grid                As the WSV with the same name.
   \param   refellipsoid          As the WSV with the same name.
   \param   z_field               As the WSV with the same name.
   \param   t_field               As the WSV with the same name.
   \param   vmr_field             As the WSV with the same name.
   \param   f_grid                As the WSV with the same name.
   \param   r                     The radius of the position of interest.
   \param   lat                   The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-14
*/
void refr_gradients_1d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       Numeric& dndr,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r) {
  get_refr_index_1d(ws,
                    refr_index_air,
                    refr_index_air_group,
                    refr_index_air_agenda,
                    p_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r);

  const Numeric n0 = refr_index_air;
  Numeric dummy;

  get_refr_index_1d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r + 1);

  dndr = refr_index_air - n0;

  refr_index_air = n0;
}

//! refr_gradients_2d
/*! 
   Determines the refractive index, and its gradients, for the given position.

   The gradients are calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius or
   latitude and the difference to the refractive index at the given
   point determines the gradient.

   The latitude gradient is scaled with the radius to obtain the same
   unit ([1/m]) for both gradients. That is, the returned value is the
   change of the refractive index for a movement of 1m in the latitude
   direction.

   \param   ws                    Current Workspace
   \param   refr_index_air        Output: As the WSV with the same name.
   \param   refr_index_air_group  Output: As the WSV with the same name.
   \param   dndr                  Output: Radial gradient of refractive index.
   \param   dndlat                Output: Latitude gradient of refractive index.
   \param   refr_index_air_agenda As the WSV with the same name.
   \param   p_grid                As the WSV with the same name.
   \param   lat_grid              As the WSV with the same name.
   \param   refellipsoid          As the WSV with the same name.
   \param   z_field               As the WSV with the same name.
   \param   t_field               As the WSV with the same name.
   \param   vmr_field             As the WSV with the same name.
   \param   f_grid                As the WSV with the same name.
   \param   r                     The radius of the position of interest.
   \param   lat                   The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-14
*/
void refr_gradients_2d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       Numeric& dndr,
                       Numeric& dndlat,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat) {
  get_refr_index_2d(ws,
                    refr_index_air,
                    refr_index_air_group,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r,
                    lat);

  const Numeric n0 = refr_index_air;
  Numeric dummy;

  get_refr_index_2d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r + 1,
                    lat);

  dndr = refr_index_air - n0;

  const Numeric dlat = 1e-4;

  get_refr_index_2d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r,
                    lat + dlat);

  dndlat = (refr_index_air - n0) / (DEG2RAD * dlat * r);

  refr_index_air = n0;
}

//! refr_gradients_3d
/*! 
   Determines the refractive index, and its gradients, for the given position.

   The gradients are calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius,
   latitude or longitude and the difference to the refractive index at
   the given point determines the gradient.

   The latitude and longitude gradients are scaled with the
   (effective) radius to obtain the same unit ([1/m]) for all
   gradients. That is, the returned values are the change of the
   refractive index for a movement of 1m in the latitude or longitude
   direction.

   \param   ws                   Current Workspace
   \param   refr_index_air       Output: As the WSV with the same name.
   \param   refr_index_air_group Output: As the WSV with the same name.
   \param   dndr                 Output: Radial gradient of refractive index.
   \param   dndlat               Output: Latitude gradient of refractive index.
   \param   dndlon               Output: Longitude gradient of refractive index.
   \param   refr_index_air_agenda As the WSV with the same name.
   \param   p_grid               As the WSV with the same name.
   \param   lat_grid             As the WSV with the same name.
   \param   lon_grid             As the WSV with the same name.
   \param   refellipsoid         As the WSV with the same name.
   \param   z_field              As the WSV with the same name.
   \param   t_field              As the WSV with the same name.
   \param   vmr_field            As the WSV with the same name.
   \param   f_grid               As the WSV with the same name.
   \param   r                    The radius of the position of interest.
   \param   lat                  The latitude of the position of interest.
   \param   lon                  The longitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-17
*/
void refr_gradients_3d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       Numeric& dndr,
                       Numeric& dndlat,
                       Numeric& dndlon,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView lon_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat,
                       const Numeric& lon) {
  get_refr_index_3d(ws,
                    refr_index_air,
                    refr_index_air_group,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r,
                    lat,
                    lon);

  const Numeric n0 = refr_index_air;
  Numeric dummy;

  get_refr_index_3d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r + 1,
                    lat,
                    lon);

  dndr = refr_index_air - n0;

  const Numeric dlat = 1e-4;

  get_refr_index_3d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r,
                    lat + dlat,
                    lon);

  dndlat = (refr_index_air - n0) / (DEG2RAD * dlat * r);

  const Numeric dlon = 1e-4;

  get_refr_index_3d(ws,
                    refr_index_air,
                    dummy,
                    refr_index_air_agenda,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    refellipsoid,
                    z_field,
                    t_field,
                    vmr_field,
                    f_grid,
                    r,
                    lat,
                    lon + dlon);

  dndlon = (refr_index_air - n0) / (DEG2RAD * dlon * r * cos(DEG2RAD * lat));

  refr_index_air = n0;
}

void refractive_index_water_and_steam_VisNIR(Numeric& n,
                                             const Index& only_valid_range,
                                             const Numeric& frequency,
                                             const Numeric& temperature,
                                             const Numeric& density) {
  //convert frequency to wavelength
  const Numeric wavelength = Conversion::freq2wavelen(frequency) * 1e6;  // [µm]

  //Reference values
  const Numeric T_star = 273.15;      // [K]
  const Numeric rho_star = 1000;      // [kg/m^3]
  const Numeric lambda_star = 0.589;  // [µm]

  //check input
  bool T_ok = (temperature > T_star - 12) && (temperature < T_star + 500);
  bool rho_ok = (density > 0) && (density < 1060);
  bool wvl_ok = (wavelength > 0.2) && (wavelength < 1.9);

  if (only_valid_range) {
    ARTS_USER_ERROR_IF(!T_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Temperature must be between ",
                       T_star - 12,
                       "K and ",
                       T_star + 500,
                       "K\n"
                       "Desired temperature: ",
                       temperature,
                       "K \n")

    ARTS_USER_ERROR_IF(!rho_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Density must be between ",
                       0,
                       "kg m^-3 and ",
                       1060,
                       "kg m^-3\n"
                       "Desired density: ",
                       density,
                       "kg m^-3 \n")

    Numeric frq_upper_limit = Conversion::wavelen2freq(0.2 * 1e-6);
    Numeric frq_lower_limit = Conversion::wavelen2freq(1.9 * 1e-6);

    ARTS_USER_ERROR_IF(!wvl_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Frequency must be between ",
                       frq_lower_limit,
                       "Hz and ",
                       frq_upper_limit,
                       "Hz\n"
                       "Desired density: ",
                       frequency,
                       "Hz \n")
  }

  //coefficients
  const Numeric a0 = 0.244257733;
  const Numeric a1 = 9.74634476e-3;
  const Numeric a2 = -3.73234996e-3;
  const Numeric a3 = 2.68678472e-4;
  const Numeric a4 = 1.58920570e-3;
  const Numeric a5 = 2.45934259e-3;
  const Numeric a6 = 0.900704920;
  const Numeric a7 = -1.66626219e-2;

  const Numeric lambda_uv = 0.2292020;
  const Numeric lambda_ir = 5.432937;

  //normalize input variables
  Numeric T_bar = temperature / T_star;
  Numeric rho_bar = density / rho_star;
  Numeric lambda_bar = wavelength / lambda_star;

  //set up right hands side of eq A1
  Numeric rhs = a0;
  rhs += a1 * rho_bar;
  rhs += a2 * T_bar;
  rhs += a3 * lambda_bar * lambda_bar * T_bar;
  rhs += a4 / (lambda_bar * lambda_bar);
  rhs += a5 / (lambda_bar * lambda_bar - lambda_uv * lambda_uv);
  rhs += a6 / (lambda_bar * lambda_bar - lambda_ir * lambda_ir);
  rhs += a7 * rho_bar * rho_bar;

  Complex a;
  a = rho_bar * rhs;

  //solve Eq A1 (see paper in documentation after n
  Complex n_cmplx = sqrt(-2 * a - 1);
  n_cmplx /= sqrt(a - 1);

  //refractive index
  n = - real(n_cmplx);
}