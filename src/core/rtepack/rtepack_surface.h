#pragma once

#include <geodetic.h>

#include <span>

#include "rtepack_mueller_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** Reflect and emit a Stokes vector
 * 
 * The reflection is scalar, i.e. the same for all Stokes components.
 * To maintain R + E = 1, (1 - R) * B is added as emission.  Named B
 * because this is often the Planck emission
 * 
 * Circular polarization is mirrored.
 * 
 * @param I Incoming radiation
 * @param R Scalar reflection
 * @param B Vector emission
 * @return  Outgoing radiation - R I + (1 - R) B
 */
stokvec flat_scalar_reflection(stokvec I, const Numeric R, const stokvec B);

/** First order partial derivative wrt reflection of flat_scalar_reflection
 * 
 * @param I Incoming radiation
 * @param B Vector emission
 * @return  Outgoing first order radiation
 */
stokvec dflat_scalar_reflection_dr(stokvec I, const stokvec B);

/** Fresnel style reflectance matrix
 *
 *  @param Rv Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh Reflection AMPLITUDE coefficient for horizontal polarisation.
 */
muelmat fresnel_reflectance(Complex Rv, Complex Rh);

/** Fresnel reflection Mueller matrix for specular reflection off a
 *  surface with arbitrary tilt.
 *
 * UNTESTED
 *
 *  @param Rv        Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh        Reflection AMPLITUDE coefficient for horizontal polarisation.
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param n_surface Outward unit surface normal
 *  @return          4x4 Mueller matrix in the incident beam's polarization frame
 *
 *  For a flat horizontal surface (n_surface = (0,0,1)) this reduces to
 *  fresnel_reflectance(Rv, Rh) with U and V sign flips.
 */
muelmat fresnel_reflectance_specular(Complex Rv,
                                     Complex Rh,
                                     const Vector3& k_inc,
                                     const Vector3& n_surface);

/** Fresnel reflection Mueller matrix for non-specular (e.g. Lambertian,
 *  BRDF, or tilted-facet) reflection where incident and outgoing directions
 *  are independent.
 *
 * UNTESTED
 *
 *  @param Rv        Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param Rh        Reflection AMPLITUDE coefficient for horizontal polarisation.
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param k_out     Unit propagation vector of outgoing beam (toward observer)
 *  @param n_surface Outward unit surface normal
 *  @return          4x4 Mueller matrix mapping incident Stokes to outgoing Stokes
 */
muelmat fresnel_reflectance_nonspecular(Complex Rv,
                                        Complex Rh,
                                        const Vector3& k_inc,
                                        const Vector3& k_out,
                                        const Vector3& n_surface);

/** Reflect and emit a Stokes vector
 * 
 * The reflection is scalar, i.e. the same for all Stokes components.
 * To maintain R + E = 1, (1 - R) * B is added as emission.  Named B
 * because this is often the Planck emission
 * 
 * Circular polarization is mirrored.
 * 
 * @param I Incoming radiation
 * @param R Matrix reflection
 * @param B Vector emission
 * @return  Outgoing radiation - R I + (1 - R) B
 */
stokvec reflection(stokvec I, const muelmat R, const stokvec B);
stokvec dreflection(stokvec I, const muelmat dR, const stokvec B);

/** Return the outgoing propagation direction for specular reflection.
     *
     *  Both `k_inc` and `n_surface` are assumed to be unit vectors with
     *  `k_inc` pointing toward the surface and `n_surface` pointing outwards.
     */
Vector3 specular_reflected_direction(const Vector3& k_inc,
                                     const Vector3& n_surface);

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
stokvec specular_radiance(const stokvec& I_in,
                          const stokvec& J,
                          Complex Rv,
                          Complex Rh,
                          const Vector3& k_inc,
                          const Vector3& n_surface);

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
                             Complex Rv,
                             Complex Rh,
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
 *  @param coords     Visible (lat, lon) pairs [degrees] from visible_coordinates()
 *  @param sources    Source radiance (stokvec) at each visible coordinate
 *  @param J          Thermal emission Stokes vector at the scatter point
 *  @param Rv         Complex Fresnel amplitude for vertical polarisation
 *  @param Rh         Complex Fresnel amplitude for horizontal polarisation
 *  @param pos        Scatter-point (lat, lon) [degrees]
 *  @param h_pos      Scatter-point height [m]
 *  @param n_surface  Outward unit normal at the scatter point
 *  @param k_out      Unit propagation direction toward the sensor
 *  @param ellipsoid  Ellipsoid semi-axes (a, b) [m]
 *  @param hfield     Height field — provides grid spacing and patch heights
 *  @return           Outgoing Stokes vector toward the sensor
 */
stokvec nonspecular_radiance_from_patches(std::span<const Vector2> coords,
                                          stokvec_vector_const_view sources,
                                          const stokvec& J,
                                          Complex Rv,
                                          Complex Rh,
                                          Vector2 pos,
                                          Numeric h_pos,
                                          const Vector3& n_surface,
                                          const Vector3& k_out,
                                          Vector2 ellipsoid,
                                          const GeodeticField2& hfield);
}  // namespace rtepack
