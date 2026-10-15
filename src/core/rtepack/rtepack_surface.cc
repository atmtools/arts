#include "rtepack_surface.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <debug.h>

#include <span>

#include "rtepack_multitype.h"

namespace rtepack {

muelmat power_reflectance(Numeric rv, Numeric rh) {
  const Numeric mean       = 0.5 * (rv + rh);
  const Numeric difference = 0.5 * (rv - rh);
  const Numeric product    = -std::sqrt(rv * rh);
  return {mean, difference, 0, 0, difference, mean, 0, 0, 0, 0, product, 0, 0, 0, 0, product};
}

muelmat fresnel_reflectance_specular(Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& n_surface) {
  return fresnel_reflectance_nonspecular(Rv, Rh, k_inc, specular_reflected_direction(k_inc, n_surface), n_surface);
}

muelmat fresnel_reflectance_nonspecular(
    Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& k_out, const Vector3& n_surface) {
  /* The fields of both waves in their (p, s) frames, s across the plane of
     incidence and p = s x k, in which the reflection is
     E_p' = Rv E_p, E_s' = Rh E_s with the amplitudes of fresnel()
     (Rv -> +1 and Rh -> -1 for a perfect conductor).  Both (p, s, k) are
     right-handed, as the meridional bases (v, h, k) are, so the frames are
     related by rotations alone.  At normal incidence the plane is undefined
     and any s across k_inc will do. */
  const Vector3 bisector = k_out - k_inc;
  ARTS_USER_ERROR_IF(hypot(cross(n_surface, bisector)) > 1e-9 * hypot(bisector) or dot(n_surface, bisector) <= 0.0,
                     R"(The surface normal {} does not bisect -k_inc = {} and k_out = {}.

A reflection exists only for the facet whose normal is along k_out - k_inc.
)",
                     n_surface,
                     -1.0 * k_inc,
                     k_out);

  const Vector3            x{1.0, 0.0, 0.0};
  const polarization_basis in  = meridional_basis(k_inc, x);
  const polarization_basis out = meridional_basis(k_out, x);

  Vector3       s      = cross(k_inc, n_surface);
  const Numeric norm_s = hypot(s);
  if (norm_s < 1e-12) {
    s = in.second;
  } else {
    s /= norm_s;
  }

  return basis_rotation({cross(s, k_out), s}, out) * fresnel_reflectance(Rv, Rh) *
         basis_rotation(in, {cross(s, k_inc), s});
}

Vector3 specular_reflected_direction(const Vector3& k_inc, const Vector3& n_surface) {
  const Vector3 out = k_inc - 2.0 * dot(k_inc, n_surface) * n_surface;
  return normalized(out);
}

stokvec specular_radiance(
    const stokvec& I_in, const stokvec& J, Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& n_surface) {
  const muelmat R = fresnel_reflectance_specular(Rv, Rh, k_inc, n_surface);
  return J + R * (I_in - J);
}

stokvec nonspecular_radiance(const stokvec& I_in,
                             const stokvec& J,
                             Complex        Rv,
                             Complex        Rh,
                             const Vector3& k_inc,
                             const Vector3& k_out,
                             const Vector3& n_surface) {
  const muelmat R = fresnel_reflectance_nonspecular(Rv, Rh, k_inc, k_out, n_surface);
  return J + R * (I_in - J);
}

stokvec nonspecular_radiance_from_patches(std::span<const Vector2>  coords,
                                          stokvec_vector_const_view sources,
                                          const stokvec&            J,
                                          Complex                   Rv,
                                          Complex                   Rh,
                                          Vector2                   pos,
                                          Numeric                   h_pos,
                                          const Vector3&            n_surface,
                                          const Vector3&            k_out,
                                          Vector2                   ellipsoid,
                                          const GeodeticField2&     hfield) {
  constexpr Numeric pi = Constant::pi;
  using Conversion::deg2rad;

  assert(coords.size() == sources.size());

  const Vector& lats     = hfield.grid<0>();
  const Vector& lons     = hfield.grid<1>();
  const Numeric dlat_rad = deg2rad(lats[1] - lats[0]);
  const Numeric dlon_rad = deg2rad(lons[1] - lons[0]);
  const Numeric dlat     = lats[1] - lats[0];
  const Numeric dlon     = lons[1] - lons[0];
  const Index   nlat     = static_cast<Index>(lats.size());
  const Index   nlon     = static_cast<Index>(lons.size());
  const Numeric a        = ellipsoid[0];  // semi-major axis as effective radius

  // ECEF position of the scatter point
  const Vector3 pos_P = geodetic2ecef({h_pos, pos[0], pos[1]}, ellipsoid);

  // The Stokes vectors are in the meridional bases of the local vertical at P, so the reflection
  // takes its directions in P's east-north-up frame, whose up is the ellipsoid normal
  const Numeric sin_lat = std::sin(deg2rad(pos[0])), cos_lat = std::cos(deg2rad(pos[0]));
  const Numeric sin_lon = std::sin(deg2rad(pos[1])), cos_lon = std::cos(deg2rad(pos[1]));
  const Vector3 east{-sin_lon, cos_lon, 0.0};
  const Vector3 north{-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat};
  const Vector3 up{cos_lat * cos_lon, cos_lat * sin_lon, sin_lat};
  const auto    local       = [&](const Vector3& v) { return Vector3{dot(v, east), dot(v, north), dot(v, up)}; };
  const Vector3 k_out_local = local(k_out);

  stokvec result = J;

  for (Index idx = 0; idx < static_cast<Index>(coords.size()); ++idx) {
    const Numeric lat_j = coords[idx][0];
    const Numeric lon_j = coords[idx][1];

    // Height at this patch (nearest-grid lookup — coords come from the grid)
    const Index   ilat_j = std::clamp(static_cast<Index>(std::lround((lat_j - lats[0]) / dlat)), Index(0), nlat - 1);
    const Index   ilon_j = std::clamp(static_cast<Index>(std::lround((lon_j - lons[0]) / dlon)), Index(0), nlon - 1);
    const Numeric h_j    = hfield.data[ilat_j, ilon_j];

    // ECEF position of patch j and vector toward scatter point
    const Vector3 pos_j = geodetic2ecef({h_j, lat_j, lon_j}, ellipsoid);
    const Vector3 r_vec = pos_P - pos_j;
    const Numeric r     = hypot(r_vec);
    if (r < 1.0) continue;            // degenerate: same point
    const Vector3 k_inc = r_vec / r;  // propagation direction: from j toward P

    // Outward surface normal at patch j (ellipsoid normal)
    const Numeric latrad_j = deg2rad(lat_j);
    const Numeric lonrad_j = deg2rad(lon_j);
    const Vector3 n_j{
        std::cos(latrad_j) * std::cos(lonrad_j), std::cos(latrad_j) * std::sin(lonrad_j), std::sin(latrad_j)};

    // Emission angle at patch j (angle between j's normal and direction toward P)
    const Numeric cos_alpha = dot(n_j, k_inc);
    if (cos_alpha <= 0.0) continue;  // j not emitting toward P

    // Incidence angle at scatter point P
    // (k_inc points INTO the surface; n_surface points outward)
    const Numeric cos_theta = -dot(n_surface, k_inc);
    if (cos_theta <= 0.0) continue;  // behind the surface at P

    // Cell area on the ellipsoid surface [m^2]
    const Numeric R_j = a + h_j;
    const Numeric A_j = R_j * R_j * std::abs(dlat_rad) * std::abs(dlon_rad) * std::abs(std::cos(latrad_j));

    // Solid angle of patch j as seen from P [sr]
    const Numeric dOmega = A_j * cos_alpha / (r * r);

    // Discrete contribution to the reflectance integral:
    //   dL = (1/π) R(k_inc, k_out) · L_j · cos_theta · dΩ
    // with R the reflection by the facet that turns k_inc into k_out
    const Vector3 k_inc_local = local(k_inc);
    Vector3       facet       = k_out_local - k_inc_local;
    const Numeric norm_facet  = hypot(facet);
    if (norm_facet == 0.0) continue;  // grazing: no facet reflects k_inc into k_out
    facet           /= norm_facet;
    const muelmat R  = fresnel_reflectance_nonspecular(Rv, Rh, k_inc_local, k_out_local, facet);
    result          += (cos_theta * dOmega / pi) * R * sources[idx];
  }

  return result;
}

}  // namespace rtepack
