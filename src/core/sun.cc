/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   sun.cc
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Functions needed for radiative transfer with direct sources.
*/

#include "sun.h"

#include "arts_conversions.h"
#include "check_input.h"
#include "debug.h"
#include "interpolation.h"
#include "matpack_data.h"
#include "physics_funcs.h"

using Constant::pi;

/*===========================================================================
  === The functions
  ===========================================================================*/

std::ostream& operator<<(std::ostream& os, const Sun& sun) {
  os << "Sun: " << sun.description;
  os << " Radius: " << sun.radius << "m ";
  os << " Distance: " << sun.distance << "m \n";
  os << " Latitude: " << sun.latitude << "° \n";
  os << " Longitude: " << sun.longitude << "° \n";
  os << " Spectrum [W/m2/Hz]: \n" << sun.spectrum;
  return os;
}

Matrix regrid_sun_spectrum(const GriddedField2& sun_spectrum_raw,
                           const Vector& f_grid,
                           const Numeric& temperature) {
  const Index nf = f_grid.size();
  const Vector& data_f_grid = sun_spectrum_raw.grid<0>();
  const Numeric data_fmin = data_f_grid[0];
  const Numeric data_fmax = data_f_grid[data_f_grid.size() - 1];

  // Result array
  Matrix int_data(f_grid.size(), 4, 0.);

  const Numeric* f_grid_begin = f_grid.unsafe_data_handle();
  const Numeric* f_grid_end = f_grid_begin + f_grid.size();
  const Index i_fstart = std::distance(
      f_grid_begin, std::lower_bound(f_grid_begin, f_grid_end, data_fmin));
  const Index i_fstop =
      std::distance(
          f_grid_begin,
          std::upper_bound(f_grid_begin + i_fstart, f_grid_end, data_fmax)) -
      1;

  const Index f_extent = i_fstop - i_fstart + 1;

  // This is the part of f_grid that lies inside the spectrum data band
  const ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

  // Determine the index of the first grid points in the spectrum band that
  // lies outside or on the boundary of the part of f_grid that lies inside
  // the spectral band.
  const Numeric f_grid_fmin = f_grid[i_fstart];
  const Numeric f_grid_fmax = f_grid[i_fstop];

  const Numeric* data_f_grid_begin = data_f_grid.unsafe_data_handle();
  const Numeric* data_f_grid_end = data_f_grid_begin + data_f_grid.size() - 1;
  const Index i_data_fstart =
      std::distance(
          data_f_grid_begin,
          std::upper_bound(data_f_grid_begin, data_f_grid_end, f_grid_fmin)) -
      1;
  const Index i_data_fstop = std::distance(
      data_f_grid_begin,
      std::upper_bound(
          data_f_grid_begin + i_data_fstart, data_f_grid_end, f_grid_fmax));

  // Extent for active data frequency vector:
  const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

  // For this range the interpolation will find place.
  const Range active_range(i_data_fstart, data_f_extent);
  const ConstVectorView data_f_grid_active = data_f_grid[active_range];

  // Check if frequency is inside the range covered by the data:
  chk_interpolation_grids(
      "Frequency interpolation for cross sections", data_f_grid, f_grid_active);

  {
    // Find frequency grid positions:
    ArrayOfGridPos f_gp(f_grid_active.size());
    gridpos(f_gp, data_f_grid_active, f_grid_active, 0);

    Matrix itw(f_gp.size(), 2);
    interpweights(itw, f_gp);

    for (int i = 0; i < 4; i++) {
      interp(int_data(Range(i_fstart, f_extent), i),
             itw,
             sun_spectrum_raw.data(active_range, i),
             f_gp);
    }
  }

  // Padding
  if (f_extent < nf) {
    if (temperature == -1) {
      ARTS_USER_ERROR(
          "f_grid is (partially) outside the sun spectrum data"
          "Set temperature to zero to have a padding of "
          "0 or a fitting effective temperature"
          "For further information take a look at the "
          "documentation for regrid_sun_spectrum")
    }
    if (temperature > 0) {
      for (int i = 0; i < i_fstart; i++) {
        int_data(i, 0) = planck(f_grid[i], temperature);
      }
      for (Index i = f_extent; i < nf; i++) {
        int_data(i, 0) = planck(f_grid[i], temperature);
      }
    }
  }

  return int_data;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfSun& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

/*!
    Reference ellipsoid radius, directly from *refellipsoid*.

    Gives the distance from the Earth's centre and the reference ellipsoid as a
    function of geoCENTRIC latitude.

    For 1D, extract r directly as refellipsoid[0] to save time.

    This is the basic function to calculate the reference ellipsoid radius.
    However, inside the atmosphere this radius is just used at the positions of
    the lat_grid. A linear interpolation is applied between these points. This
    is handled by other functions. For 2D and 3D and the grid position is
    known, use *refell2d*. The function pos2refell_r handles all this in a
    general way (but not always the fastest option).

    \return                 Ellispoid radius
    \param  refellipsoid    [a, b]
    \param  latitude        A geoecentric latitude.

    \author Patrick Eriksson 
    \date   2012-02-07

    Update for new definition of ellipsoid and interface
    \author Richard Larsson
    \date   2024-05-15
*/
Numeric refell2r(const Vector2 ell, const Numeric lat) {
  using Conversion::cosd;
  using Conversion::sind;

  const auto [a, b] = ell;

  ARTS_ASSERT(a > 0);
  ARTS_ASSERT(b > 0);
  ARTS_ASSERT(a >= b);

  const Numeric e2 = 1 - Math::pow2(b / a);
  const Numeric c = 1 - e2;

  const Numeric ct = cosd(lat);
  const Numeric st = sind(lat);

  return b / std::sqrt(c * ct * ct + st * st);
}

/*! 
   Conversion from spherical to cartesian coordinates.

   The cartesian coordinate system is defined such as the x-axis goes along
   lat=0 and lon=0, the z-axis goes along lat=0 and lon=90, and z-axis goes
   along lat=90. 

   \param   sph [r, lat, lon] Radius, Latitude, Longitude.
   \return [x, y, z] Cartesian coordinates.

   \author Patrick Eriksson
   \date   2002-12-17

   Update for new definition of ellipsoid and interface
   \author Richard Larsson
   \date   2024-05-15
*/
Vector3 sph2cart(const Vector3 sph) {
  using Conversion::cosd;
  using Conversion::sind;
  
  const auto& [r, lat, lon] = sph;

  ARTS_ASSERT(std::abs(lat) <= 90);
  ARTS_ASSERT(std::abs(lon) <= 360);
  ARTS_ASSERT(r > 0);

  return {r * cosd(lat) * cosd(lon), r * cosd(lat) * sind(lon), r * sind(lat)};
}

/*! 
   Conversion from position and LOS to cartesian coordinates

   A position (in geographical coordinates) and LOS are converted to a
   cartesian position and a viewing vector. The viewing direction is the
   the vector [dx,dy,dz]. This vector is normalised (it has length 1).

   See the user guide for definition on the zenith and azimuth angles.

   \param   sph [r, lat, lon] Radius, Latitude, Longitude.
   \param   los [za, aa] Zenith and azimuth angles.
   \return [x, y, z] Cartesian coordinates. [dx, dy, dz] Viewing vector.

   \author Patrick Eriksson
   \date   2002-12-30

   Update for new definition of ellipsoid and interface
   \author Richard Larsson
   \date   2024-05-15
*/
std::pair<Vector3, Vector3> poslos2cart(const Vector3 sph, const Vector2 los) {
  using Conversion::cosd;
  using Conversion::sind;

  const auto [r, lat, lon] = sph;
  const auto [za, aa] = los;

  ARTS_ASSERT(r > 0);
  ARTS_ASSERT(std::abs(lat) <= 90);
  ARTS_ASSERT(za >= 0 && za <= 180);

  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if (std::abs(lat) > 90 - 1e-8) {
    const Numeric s = lat > 0 ? 1 : -1;

    return {{0, 0, s * r},
            {sind(za) * cosd(aa), sind(za) * sind(aa), s * cosd(za)}};
  }

  const Numeric coslat = cosd(lat);
  const Numeric sinlat = sind(lat);
  const Numeric coslon = cosd(lon);
  const Numeric sinlon = sind(lon);
  const Numeric cosza = cosd(za);
  const Numeric sinza = sind(za);
  const Numeric cosaa = cosd(aa);
  const Numeric sinaa = sind(aa);

  const Numeric dr = cosza;
  const Numeric dlat = sinza * cosaa;  // r-term cancel out below
  const Numeric dlon = sinza * sinaa / coslat;

  return {
      {r * coslat * coslon, r * coslat * sinlon, r * sinlat},
      {coslat * coslon * dr - sinlat * coslon * dlat - coslat * sinlon * dlon,
       coslat * sinlon * dr - sinlat * sinlon * dlat + coslat * coslon * dlon,
       sinlat * dr + coslat * dlat}};
}

std::pair<Numeric, bool> hit_sun(const Sun& sun,
                                 const Vector3 pos,
                                 const Vector2 los,
                                 const Vector2 ell) {
  //Calculate earth centric radial component of sun_pos and rtp_pos.
  const Numeric R_sun = sun.distance;
  const Numeric R_rte = refell2r(ell, pos[1]) + pos[0];

  // r_sun
  const auto r_sun = sph2cart({R_sun, sun.latitude, sun.longitude});

  // r_rte, r_los
  const auto [r_rte, r_los] = poslos2cart({R_rte, pos[1], pos[2]}, los);

  const auto& [r_los_x, r_los_y, r_los_z] = r_los;

  //Calculate vector of line of sight and unit vector pointing from
  //ppath point to the sun.
  const auto [r_ps_x, r_ps_y, r_ps_z] = r_sun - r_rte;

  //abs value of r_ps
  const Numeric r_ps = std::hypot(r_ps_x, r_ps_y, r_ps_z);

  //abs value of r_los
  const Numeric r_glos = std::hypot(r_los_x, r_los_y, r_los_z);

  //Calculate angle beta between line of sight and the line between ppath point and the sun
  //using scalar product
  const Numeric cos_beta =
      std::clamp((r_ps_x * r_los_x + r_ps_y * r_los_y + r_ps_z * r_los_z) /
                     (r_ps * r_glos),
                 -1.0,
                 1.0);
  const Numeric beta = std::acos(cos_beta);

  // angular radius of sun
  const Numeric alpha = std::atan2(sun.radius, r_ps);

  return {beta, beta <= alpha};
}

Numeric Sun::sin_alpha_squared(Vector3 pos, Vector2 ell) const try {
  using Math::pow2;

  const auto[alt, lat, lon] = pos;

  // r_sun
  const Numeric R_rte = refell2r(ell, lat) + alt;
  const auto r_sun = sph2cart({distance, latitude, longitude});
  const auto r_rte = sph2cart({R_rte, lat, lon});

  const auto d = r_sun - r_rte;

  //abs value of r_ps
  const Numeric r_ps2 = d * d;

  return pow2(radius) / (pow2(radius) + r_ps2);
} ARTS_METHOD_ERROR_CATCH

//! Returns true if the sun was hit, false otherwise.  Sets the spectral radiance if the sun was hit.
bool set_spectral_radiance_if_sun_intersection(
    StokvecVector& spectral_radiance,
    const Sun& sun,
    const PropagationPathPoint& propagation_path_point,
    const SurfaceField& surface_field) {
  const Index nf = spectral_radiance.size();
  ARTS_ASSERT(
      nf == sun.spectrum.nrows(),
      "Spectral radiance and sun spectrum must have the same frequency size")
  ARTS_ASSERT(
      4 == sun.spectrum.ncols(),
      "Spectral radiance and sun spectrum must have the same stokes size")

  //Check if we see the sun.
  if (propagation_path_point.los_type == PathPositionType::space and
      hit_sun(sun,
              propagation_path_point.pos,
              path::mirror(propagation_path_point.los),
              surface_field.ellipsoid)
          .second) {
    //Here we assume that the sun radiates isotropically.
    //Matrix sun_irradiance = sun.spectrum;

    for (Index iv = 0; iv < nf; ++iv) {
      for (Index is = 0; is < 4; ++is) {
        spectral_radiance[iv][is] = sun.spectrum(iv, is) / Constant::pi;
      }
    }

    return true;
  }

  return false;
}
