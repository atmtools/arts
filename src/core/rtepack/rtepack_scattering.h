#pragma once

#include <matpack.h>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** Polarized monostatic-radar return.
 *
 * The order is explicit because polarized outgoing and return transmissions
 * need neither be equal nor commute.
 */
constexpr stokvec radar_return(const muelmat &return_transmission,
                               const muelmat &backscatter,
                               const muelmat &outgoing_transmission,
                               const stokvec &transmitted) {
  return return_transmission * (backscatter * (outgoing_transmission * transmitted));
}

/** Directional derivative of radar_return.
 *
 * Implements the complete product rule for
 * ``Tr * Z * Tf * transmitted``.  Derivatives of the transmitted Stokes
 * vector are deliberately outside this atmospheric-optics helper.
 */
constexpr stokvec radar_return_derivative(const muelmat &return_transmission,
                                          const muelmat &return_transmission_derivative,
                                          const muelmat &backscatter,
                                          const muelmat &backscatter_derivative,
                                          const muelmat &outgoing_transmission,
                                          const muelmat &outgoing_transmission_derivative,
                                          const stokvec &transmitted) {
  const stokvec outgoing             = outgoing_transmission * transmitted;
  const stokvec outgoing_derivative  = outgoing_transmission_derivative * transmitted;
  const stokvec scattered            = backscatter * outgoing;
  const stokvec scattered_derivative = backscatter_derivative * outgoing + backscatter * outgoing_derivative;
  return return_transmission_derivative * scattered + return_transmission * scattered_derivative;
}

Array<muelmat_vector> bulk_backscatter(const ConstTensor5View &Pe, const ConstMatrixView &pnd);

Array<muelmat_matrix> bulk_backscatter_derivative(const ConstTensor5View &Pe, const ArrayOfMatrix &dpnd_dx);

void bulk_backscatter_commutative_transmission_rte(Array<stokvec_vector>        &I,
                                                   Array<Array<stokvec_matrix>> &dI,
                                                   const stokvec_vector         &I_incoming,
                                                   const Array<muelmat_vector>  &T,
                                                   const Array<muelmat_vector>  &PiTf,
                                                   const Array<muelmat_vector>  &PiTr,
                                                   const Array<muelmat_vector>  &Z,
                                                   const Array<muelmat_matrix>  &dT1,
                                                   const Array<muelmat_matrix>  &dT2,
                                                   const Array<muelmat_matrix>  &dZ);

muelmat rayleigh_scattering(const Vector2 &los_in, const Vector2 &los_out, const Numeric depolarization_factor);
}  // namespace rtepack
