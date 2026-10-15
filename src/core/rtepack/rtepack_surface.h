#pragma once

#include <geodetic.h>

#include <span>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_polarization.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** Fresnel style reflectance matrix
 *
 *  @param Rv Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh Reflection AMPLITUDE coefficient for horizontal polarisation.
 */
constexpr muelmat fresnel_reflectance(Complex Rv, Complex Rh) {
  const Numeric rv    = std::norm(Rv);
  const Numeric rh    = std::norm(Rh);
  const Numeric rmean = 0.5 * (rv + rh);
  const Numeric rdiff = 0.5 * (rv - rh);
  const Complex a     = Rh * std::conj(Rv);
  const Complex b     = Rv * std::conj(Rh);
  const Numeric c     = 0.5 * std::real(a + b);
  const Numeric d     = 0.5 * std::imag(a - b);

  muelmat out{};

  out[0, 0] = rmean;
  out[1, 0] = rdiff;
  out[0, 1] = rdiff;
  out[1, 1] = rmean;
  out[2, 2] = c;
  out[2, 3] = d;
  out[3, 2] = -d;
  out[3, 3] = c;

  return out;
}

/** Reflectance matrix of power reflectances, without the phase between the amplitudes
 *
 *  rv = |Rv|^2 and rh = |Rh|^2 give the I and Q block.  The phase of Rv Rh* is taken
 *  as pi, that of a loss-free Fresnel surface below the Brewster angle and of any
 *  isotropic surface at normal incidence, so the U and V diagonal is -sqrt(rv rh)
 *  and U and V do not couple.
 *
 *  @param rv Power reflectance for vertical polarisation
 *  @param rh Power reflectance for horizontal polarisation
 */
muelmat power_reflectance(Numeric rv, Numeric rh);

/** power_reflectance(r, r): an isotropic mirror of power reflectance r, diag(r, r, -r, -r) */
constexpr muelmat power_reflectance(Numeric r) { return {r, 0, 0, 0, 0, r, 0, 0, 0, 0, -r, 0, 0, 0, 0, -r}; }

/** Fresnel reflection Mueller matrix for specular reflection off a
 *  surface with arbitrary tilt: fresnel_reflectance_nonspecular into the
 *  specularly reflected direction (specular_reflected_direction).
 *
 *  @param Rv        Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh        Reflection AMPLITUDE coefficient for horizontal polarisation.
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param n_surface Outward unit surface normal
 *  @return          4x4 Mueller matrix from the incident Stokes vector to the
 *                   reflected one, each in its meridional basis
 *
 *  For a flat horizontal surface (n_surface = (0,0,1)) this is
 *  fresnel_reflectance(Rv, Rh).
 */
muelmat fresnel_reflectance_specular(Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& n_surface);

/** Fresnel reflection Mueller matrix of the facet that reflects k_inc into
 *  k_out, for incident and outgoing directions that are independent of the
 *  mean surface (e.g. the facets of a rough surface).
 *
 *  The directions are those of propagation, so k_inc points into the
 *  surface, in a local frame whose z is the local up, e.g. east-north-up.
 *  Each Stokes vector is in the meridional basis of its direction,
 *  meridional_basis(k, x): h = k x z / |k x z| and v = h x k, and along +-z
 *  h = x, the meridian of azimuth 0 in east-north-up.  Rv and Rh are those
 *  of fresnel(), Rv -> +1 and Rh -> -1 for a perfect conductor, i.e. for the
 *  fields in (p, s) with s across the plane of incidence (k_inc x n_surface)
 *  and p = s x k for both waves.  The reflection is a rotation into that
 *  frame, the Fresnel matrix, and a rotation out of it.  It is exactly the
 *  reflection of the fields by the facet, and it reverses circular
 *  polarization at normal incidence, as a mirror does.
 *
 *  @param Rv        Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh        Reflection AMPLITUDE coefficient for horizontal polarisation.
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param k_out     Unit propagation vector of outgoing beam (toward observer)
 *  @param n_surface Outward unit normal of the facet, along k_out - k_inc;
 *                   any other normal reflects k_inc elsewhere and throws
 *  @return          4x4 Mueller matrix mapping incident Stokes to outgoing Stokes
 */
muelmat fresnel_reflectance_nonspecular(
    Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& k_out, const Vector3& n_surface);

/** Reflect and emit a Stokes vector
 *
 * To maintain R + E = 1, (id - R) * B is added as emission.  Named B
 * because this is often the Planck emission
 *
 * R takes the incident Stokes vector, in the meridional basis of its
 * propagation direction, to the reflected one in the meridional basis of
 * its propagation direction.  Any change of circular polarization is
 * part of R, e.g., fresnel_reflectance reverses it at normal incidence.
 *
 * @param I Incoming radiation
 * @param R Matrix reflection
 * @param B Vector emission
 * @return  Outgoing radiation - R I + (1 - R) B
 */
constexpr stokvec reflection(stokvec I, const muelmat R, const stokvec B) {
  I  = R * I;
  I += (muelmat::id() - R) * B;
  return I;
}

constexpr stokvec dreflection(stokvec I, const muelmat dR, const stokvec B) {
  I  = dR * I;
  I -= dR * B;
  return I;
}

/** Return the outgoing propagation direction for specular reflection.
     *
     *  Both `k_inc` and `n_surface` are assumed to be unit vectors with
     *  `k_inc` pointing toward the surface and `n_surface` pointing outwards.
     */
Vector3 specular_reflected_direction(const Vector3& k_inc, const Vector3& n_surface);

/** Specular surface emission and reflection.
 *
 *  Implements I_out = J + R_spec * (I_in - J), where R_spec is the
 *  polarised Fresnel Mueller matrix for the specular incident direction,
 *  J is the subsurface thermal emission, and I_in is the incident field
 *  arriving from the specular direction (e.g. CMB or sky radiance).
 *
 *  @param I_in      Incident Stokes vector from the specular direction
 *  @param J         Subsurface thermal emission Stokes vector
 *  @param Rv        Complex Fresnel amplitude for vertical polarisation
 *  @param Rh        Complex Fresnel amplitude for horizontal polarisation
 *  @param k_inc     Unit propagation vector of the incident beam
 *  @param n_surface Outward unit surface normal
 *  @return          Outgoing Stokes vector
 */
stokvec specular_radiance(
    const stokvec& I_in, const stokvec& J, Complex Rv, Complex Rh, const Vector3& k_inc, const Vector3& n_surface);

/** Non-specular surface emission and reflection for one incident direction.
 *
 *  Implements I_out = J + R_nonspec * (I_in - J), where R_nonspec is the
 *  polarised Fresnel Mueller matrix for an arbitrary incident/outgoing pair,
 *  J is the subsurface thermal emission, and I_in is the incident field
 *  arriving from direction k_inc (e.g. radiance from a visible surface patch).
 *  Callers accumulate weighted contributions over all visible incoming
 *  directions to build the full scattered field.
 *
 *  @param I_in      Incident Stokes vector from direction k_inc
 *  @param J         Subsurface thermal emission Stokes vector
 *  @param Rv        Complex Fresnel amplitude for vertical polarisation
 *  @param Rh        Complex Fresnel amplitude for horizontal polarisation
 *  @param k_inc     Unit propagation vector of the incident beam (toward surface)
 *  @param k_out     Unit propagation vector of the outgoing beam (toward sensor)
 *  @param n_surface Outward unit surface normal
 *  @return          Outgoing Stokes vector
 */
stokvec nonspecular_radiance(const stokvec& I_in,
                             const stokvec& J,
                             Complex        Rv,
                             Complex        Rh,
                             const Vector3& k_inc,
                             const Vector3& k_out,
                             const Vector3& n_surface);

/** Accumulate non-specular scattered radiance from all visible surface patches.
 *
 *  Discretises the reflectance integral
 *    L_out = J + (1/π) ∑_j R(k̂_j, k̂_out) · L_j · cosθ_P · dΩ_j
 *  where each visible patch j contributes solid angle
 *    dΩ_j = A_j · cosα_j / r_j²
 *  with A_j the cell area on the ellipsoid, α_j the emission angle at j
 *  (angle between patch j's outward normal and the direction toward P), and
 *  θ_P the angle of incidence at the scatter point P.
 *
 *  Each patch's radiation is reflected into k_out by the facet whose normal
 *  bisects -k̂_j and k̂_out, with the same Rv and Rh for every facet.  The
 *  mean surface normal n_surface decides which patches are visible and gives
 *  cosθ_P.
 *
 *  The directions are given in ECEF.  The Stokes vectors (sources, J and the
 *  result) are in the meridional bases of the local vertical at P, the
 *  ellipsoid normal there, as the radiation arriving at P is in ARTS.
 *
 *  @param coords     Visible (lat, lon) pairs [degrees] from visible_coordinates()
 *  @param sources    Source radiance (stokvec) at each visible coordinate
 *  @param J          Thermal emission Stokes vector at the scatter point
 *  @param Rv         Complex Fresnel amplitude for vertical polarisation
 *  @param Rh         Complex Fresnel amplitude for horizontal polarisation
 *  @param pos        Scatter-point (lat, lon) [degrees]
 *  @param h_pos      Scatter-point height [m]
 *  @param n_surface  Outward unit normal at the scatter point (ECEF)
 *  @param k_out      Unit propagation direction toward the sensor (ECEF)
 *  @param ellipsoid  Ellipsoid semi-axes (a, b) [m]
 *  @param hfield     Height field — provides grid spacing and patch heights
 *  @return           Outgoing Stokes vector toward the sensor
 */
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
                                          const GeodeticField2&     hfield);
}  // namespace rtepack
