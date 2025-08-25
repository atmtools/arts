#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
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
 *  @param[out]  Rv    Reflection AMPLITUDE coefficient for vertical polarisation.
 *  @param[out]  Rh    Reflection AMPLITUDE coefficient for vertical polarisation.
 */
muelmat fresnel_reflectance(Complex Rv, Complex Rh);

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
stokvec dreflection_dn2(stokvec I, const muelmat dR, const stokvec B);
}  // namespace rtepack
