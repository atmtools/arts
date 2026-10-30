#pragma once

#include <enums.h>
#include <matpack.h>

#include <span>

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** A quadrature point of a radar range gate inside one layer.
 *
 * The layer spans level layer - 1, closest to the radar, to level layer.
 */
struct radar_sample {
  //! The far level of the layer (>= 1)
  Size layer;

  //! Fractional position from level layer - 1 towards level layer, in [0, 1]
  Numeric x;

  //! Quadrature weight
  Numeric weight;

  //! The output index of the range gate
  Size gate;
};

/** Single-scattering monostatic radar return of a stack of layers.
 *
 * Level 0 is closest to the radar.  The return at a point is Tr * Z * Tf * s,
 * where Tf is the ordered transmission from the radar to the point, Tr the
 * ordered transmission back, Z the backscatter matrix, and s the transmitted
 * Stokes vector.  The outgoing and return propagation matrices need neither
 * be equal nor commute.
 *
 * Inside a layer, the propagation and backscatter matrices are linear in the
 * position.  The transmission to a point inside a layer uses the same
 * exponential as the full layer, so the Magnus operator also applies to
 * partial layers.  There is no source term, so the linear-source options use
 * the transmission of their base option.
 *
 * The derivatives are with regards to the level quantities only.
 *
 * @param[out] I The return of every gate [ngate, nf]
 * @param[out] dI The derivative of every gate [ngate, np, nq, nf]
 * @param[in] transmitted The transmitted Stokes vector
 * @param[in] K_out The propagation matrix away from the radar [np][nf]
 * @param[in] dK_out The derivative of K_out [np][nq, nf]
 * @param[in] K_ret The propagation matrix towards the radar [np][nf]
 * @param[in] dK_ret The derivative of K_ret [np][nq, nf]
 * @param[in] Z The backscatter matrix [np][nf]
 * @param[in] dZ The derivative of Z [np][nq, nf]
 * @param[in] r The layer thicknesses [np - 1]
 * @param[in] samples The gate samples, sorted by layer
 * @param[in] option The transmittance option
 */
void radar_single_scattering(stokvec_matrix_view                    I,
                             stokvec_tensor4_view                   dI,
                             const stokvec                         &transmitted,
                             const std::span<const propmat_vector> &K_out,
                             const std::span<const propmat_matrix> &dK_out,
                             const std::span<const propmat_vector> &K_ret,
                             const std::span<const propmat_matrix> &dK_ret,
                             const std::span<const muelmat_vector> &Z,
                             const std::span<const muelmat_matrix> &dZ,
                             const ConstVectorView                 &r,
                             const std::span<const radar_sample>   &samples,
                             const TransmittanceOption              option);
}  // namespace rtepack
