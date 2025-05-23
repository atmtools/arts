#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** Reflect and emit a Stokes vector
 * 
 * The reflection is scalar, i.e. the same for all Stokes components.
 * To maintain R + E = 1, (1 - R) * B is added to the I() component.
 * 
 * Circular polarization is mirrored.
 * 
 * @param I Incoming radiation
 * @param R Scalar reflection
 * @param T Surface temperature
 * @param f Frequency of incoming radiation
 * @return  Outgoing radiation
 */
stokvec flat_scalar_reflection(stokvec I, const Numeric R, const Numeric B);

/** See flat_scalar_reflection, this is just the reflection part. */
stokvec flat_scalar_reflection_pure_reflect(stokvec I, const Numeric R);

/** See flat_scalar_reflection, this is the dRdx part. */
stokvec dflat_scalar_reflection_dr(stokvec I, const Numeric dRdx, const Numeric B);

/** See flat_scalar_reflection, this is the dBdx part. */
stokvec dflat_scalar_reflection_db(const Numeric R, const Numeric dBdx);
}  // namespace rtepack
