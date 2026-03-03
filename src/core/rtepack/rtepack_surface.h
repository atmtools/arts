#pragma once

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
}  // namespace rtepack
