#pragma once

#include <lin_alg.h>
#include <matpack.h>

#include <array>
#include <cstddef>

namespace rt4 {
/** The work arrays of rt4::radtrano and the routines it calls, for
 * n = nstokes * nummu streams and num_layers layers.
 *
 * Keep one (per thread) for repeated calls with the same streams, as over
 * frequency: they then allocate nothing.  rt4::radtrano sizes it; its
 * values between calls are unspecified, as every routine writes an array
 * before it reads it.  As in RADTRANO, an n x n matrix of two hemispheres
 * is [2, n, n] (the column-major Fortran matrix, so the matpack matrix
 * holds its transpose), and a source or a radiance is [2, n].
 */
struct rt4_workdata {
  /////////////////////////////////////////////////////////////////////////////
  // The layers and the ground: made by RADTRANO's layer loop, combined by its
  // level loop.  REFLECT(KRT), TRANS(KRT) and SOURCE(KS) of layer L are
  // reflect[L-1], trans[L-1] and source[L-1]; L = NUM_LAYERS+1 is the ground.
  /////////////////////////////////////////////////////////////////////////////

  Tensor4 reflect;  //!< [num_layers + 1, 2, n, n]
  Tensor4 trans;    //!< [num_layers + 1, 2, n, n]
  Tensor3 source;   //!< [num_layers + 1, 2, n]

  /////////////////////////////////////////////////////////////////////////////
  // RADTRANO on the streams: the quadrature weights, the sky radiance, the
  // initial sublayer of a scattering layer and its linear source (reflect1,
  // trans1, lin_source), and the atmosphere above (up) and below (down) a
  // level, of which reflect1, trans1 and source1 also hold the copy that
  // COMBINE_LAYERS combines with the next layer.
  /////////////////////////////////////////////////////////////////////////////

  Vector  quad_weights;                      //!< [nummu]
  Matrix  sky_radiance;                      //!< [2, n]
  Matrix  lin_source;                        //!< [2, n]
  Tensor3 reflect1, upreflect, downreflect;  //!< [2, n, n]
  Tensor3 trans1, uptrans, downtrans;        //!< [2, n, n]
  Matrix  source1, upsource, downsource;     //!< [2, n]

  /////////////////////////////////////////////////////////////////////////////
  // DOUBLING_INTEGRATION's linear source: T_LIN, CONST and T_CONST.
  /////////////////////////////////////////////////////////////////////////////

  Matrix t_lin, cnst, t_const;  //!< [2, n]

  /////////////////////////////////////////////////////////////////////////////
  // The scratch of DOUBLING_INTEGRATION, COMBINE_LAYERS and INTERNAL_RADIANCE,
  // one at a time: X, Y and GAMMA (COMMON /SCRATCH1/ and /SCRATCH2/), two
  // vectors (INTERNAL_RADIANCE's S and V), and LAPACK's workspace for the
  // inverse of GAMMA.
  /////////////////////////////////////////////////////////////////////////////

  Matrix       x, y, gamma;  //!< [n, n]
  Vector       xv, yv;       //!< [n]
  inv_workdata inv;          //!< n

  rt4_workdata() = default;
  rt4_workdata(Index nstokes, Index nummu, Index num_layers) { resize(nstokes, nummu, num_layers); }

  //! Sizes the arrays; this allocates only where an array grows
  void resize(Index nstokes, Index nummu, Index num_layers) {
    const Index n = nstokes * nummu;

    reflect.resize(num_layers + 1, 2, n, n);
    trans.resize(num_layers + 1, 2, n, n);
    source.resize(num_layers + 1, 2, n);

    quad_weights.resize(nummu);
    sky_radiance.resize(2, n);
    lin_source.resize(2, n);
    for (Tensor3* t : {&reflect1, &upreflect, &downreflect, &trans1, &uptrans, &downtrans}) t->resize(2, n, n);
    for (Matrix* m : {&source1, &upsource, &downsource, &t_lin, &cnst, &t_const}) m->resize(2, n);

    for (Matrix* m : {&x, &y, &gamma}) m->resize(n, n);
    xv.resize(n);
    yv.resize(n);
    inv.resize(static_cast<std::size_t>(n));
  }

  //! Whether the scratch of the routines (the last two groups) is sized for n streams
  [[nodiscard]] bool scratch_sized(Index n) const {
    const std::array<Index, 2> nn{n, n}, two_n{2, n};
    return x.shape() == nn and y.shape() == nn and gamma.shape() == nn and xv.extent(0) == n and yv.extent(0) == n and
           inv.N == static_cast<std::size_t>(n) and t_lin.shape() == two_n and cnst.shape() == two_n and
           t_const.shape() == two_n;
  }
};
}  // namespace rt4
