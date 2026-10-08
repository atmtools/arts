#pragma once

#include <lin_alg.h>
#include <matpack.h>

#include <array>
#include <cstddef>

#include "rt3_fft.h"

namespace polradtran::rt3 {
/** The work arrays of rt3::radtran and the routines it calls, for
 * n = nstokes * nummu streams, aziorder + 1 azimuth modes, num_layers
 * layers and nsl scattering sets.
 *
 * Keep one (per thread) for repeated calls, as over frequency: with the
 * same sizes they then allocate nothing.  rt3::radtran sizes it (resize),
 * and the routines it calls take it last.  DOUBLING_INTEGRATION,
 * COMBINE_LAYERS and INTERNAL_RADIANCE throw if it is not sized for their
 * streams; SCATTERING, DIRECT_SCATTERING and FOURIER_MATRIX size their own
 * scratch, as only they know how many azimuths they sample.  Its values
 * between calls are unspecified, as every routine writes an array before
 * it reads it, except fft, FFT1DR's phase table, which is kept and
 * extended as needed.  As in RADTRAN, an n x n matrix of two hemispheres
 * is [2, n, n] (the column-major Fortran matrix, so the matpack matrix
 * holds its transpose), and a source or a radiance is [2, n].
 */
struct rt3_workdata {
  /////////////////////////////////////////////////////////////////////////////
  // The scattering sets and the layers' optics: made once per call, before
  // the azimuth modes.  scatbuf[s] and directbuf[s] are the parts of SCATBUF
  // and DIRECTBUF of set s + 1, as SCATTERING and DIRECT_SCATTERING lay
  // them out; legendre_coef is GET_SCAT_SET's COEF, one set at a time.
  /////////////////////////////////////////////////////////////////////////////

  Vector       quad_weights;              //!< [nummu]
  Matrix       legendre_coef;             //!< [max(2 nummu, max nlegen + 1), 6]
  Vector       set_extinct, set_scatter;  //!< [nsl]
  Tensor7      scatbuf;                   //!< [nsl, aziorder + 1, 2, nummu, nummu, nstokes, nstokes]
  Tensor5      directbuf;                 //!< [nsl, aziorder + 1, 2, nummu, nstokes]
  ArrayOfIndex scat_nums;                 //!< [num_layers]
  Vector       extinctions, albedos;      //!< [num_layers]
  Vector       direct_level_flux;         //!< [num_layers + 1]

  /////////////////////////////////////////////////////////////////////////////
  // The scratch of SCATTERING and DIRECT_SCATTERING, one at a time: the
  // phase matrices at the azimuths they sample (SCAT_MATRIX, Fortran
  // layout) and their Fourier modes (BASIS_MATRIX); that of SUM_LEGENDRE
  // (the Legendre polynomials at the scattering angle) and of
  // FOURIER_MATRIX (REAL_VECTOR, BASIS_VECTOR); and FFT1DR's phase table,
  // kept between calls.  The routines size them.
  /////////////////////////////////////////////////////////////////////////////

  Tensor3      scat_matrix;   //!< [numpts (+ 1), 4, 4]
  Tensor3      basis_matrix;  //!< [2 aziorder + 1, 4, 4]
  Vector       legendre_p;    //!< [numlegendre + 1]
  Vector       real_vector;   //!< [numpts]
  Vector       basis_vector;  //!< [2 aziorder + 1]
  fft_workdata fft;

  /////////////////////////////////////////////////////////////////////////////
  // The layers and the ground of a mode: made by RADTRAN's layer loop,
  // combined by its level loop.  REFLECT(KRT), TRANS(KRT) and SOURCE(KS) of
  // layer L are reflect[L-1], trans[L-1] and source[L-1]; L = NUM_LAYERS+1
  // is the ground.
  /////////////////////////////////////////////////////////////////////////////

  Tensor4 reflect;  //!< [num_layers + 1, 2, n, n]
  Tensor4 trans;    //!< [num_layers + 1, 2, n, n]
  Tensor3 source;   //!< [num_layers + 1, 2, n]

  /////////////////////////////////////////////////////////////////////////////
  // RADTRAN on the streams in a mode: the set's scattering matrix
  // (SCATTER_MATRIX) and direct vector, the thermal vector, the initial
  // sublayer of a scattering layer and its sources (reflect1, trans1,
  // source1, exp_source, lin_source), the atmosphere above (up) and below
  // (down) a level, of which reflect1, trans1 and source1 also hold the copy
  // that COMBINE_LAYERS combines with the next layer, and the radiances
  // incident from the sky and the ground.
  /////////////////////////////////////////////////////////////////////////////

  Tensor5 scatter_matrix;                    //!< [4, nummu, nstokes, nummu, nstokes]
  Matrix  direct_vector, thermal_vector;     //!< [2, n]
  Matrix  exp_source, lin_source;            //!< [2, n]
  Tensor3 reflect1, upreflect, downreflect;  //!< [2, n, n]
  Tensor3 trans1, uptrans, downtrans;        //!< [2, n, n]
  Matrix  source1, upsource, downsource;     //!< [2, n]
  Matrix  sky_radiance;                      //!< [2, n]
  Matrix  ground_radiance, direct_radiance;  //!< [nummu, nstokes]

  /////////////////////////////////////////////////////////////////////////////
  // DOUBLING_INTEGRATION's sources: T_EXP, T_LIN, CONST and T_CONST.
  /////////////////////////////////////////////////////////////////////////////

  Matrix t_exp, t_lin, cnst, t_const;  //!< [2, n]

  /////////////////////////////////////////////////////////////////////////////
  // The scratch of DOUBLING_INTEGRATION, COMBINE_LAYERS and INTERNAL_RADIANCE,
  // one at a time: X, Y and GAMMA (COMMON /RT3_SCRATCH1/ and /RT3_SCRATCH2/),
  // two vectors (INTERNAL_RADIANCE's S and V), and LAPACK's workspace for
  // the inverse of GAMMA.
  /////////////////////////////////////////////////////////////////////////////

  Matrix       x, y, gamma;  //!< [n, n]
  Vector       xv, yv;       //!< [n]
  inv_workdata inv;          //!< n

  rt3_workdata() = default;
  rt3_workdata(Index nstokes, Index nummu, Index aziorder, Index num_layers, Index nsl, Index legendre_rows) {
    resize(nstokes, nummu, aziorder, num_layers, nsl, legendre_rows);
  }

  //! Sizes the arrays but the scratch of the scattering routines; this
  //! allocates only where an array grows
  void resize(Index nstokes, Index nummu, Index aziorder, Index num_layers, Index nsl, Index legendre_rows) {
    const Index n = nstokes * nummu;

    quad_weights.resize(nummu);
    legendre_coef.resize(legendre_rows, 6);
    set_extinct.resize(nsl);
    set_scatter.resize(nsl);
    scatbuf.resize(nsl, aziorder + 1, 2, nummu, nummu, nstokes, nstokes);
    directbuf.resize(nsl, aziorder + 1, 2, nummu, nstokes);
    scat_nums.resize(num_layers);
    extinctions.resize(num_layers);
    albedos.resize(num_layers);
    direct_level_flux.resize(num_layers + 1);

    reflect.resize(num_layers + 1, 2, n, n);
    trans.resize(num_layers + 1, 2, n, n);
    source.resize(num_layers + 1, 2, n);

    scatter_matrix.resize(4, nummu, nstokes, nummu, nstokes);
    for (Tensor3* t : {&reflect1, &upreflect, &downreflect, &trans1, &uptrans, &downtrans}) t->resize(2, n, n);
    for (Matrix* m : {&direct_vector,
                      &thermal_vector,
                      &exp_source,
                      &lin_source,
                      &source1,
                      &upsource,
                      &downsource,
                      &sky_radiance,
                      &t_exp,
                      &t_lin,
                      &cnst,
                      &t_const})
      m->resize(2, n);
    ground_radiance.resize(nummu, nstokes);
    direct_radiance.resize(nummu, nstokes);

    for (Matrix* m : {&x, &y, &gamma}) m->resize(n, n);
    xv.resize(n);
    yv.resize(n);
    inv.resize(static_cast<std::size_t>(n));
  }

  //! Whether the scratch of DOUBLING_INTEGRATION, COMBINE_LAYERS and
  //! INTERNAL_RADIANCE (the last two groups) is sized for n streams
  [[nodiscard]] bool scratch_sized(Index n) const {
    const std::array<Index, 2> nn{n, n}, two_n{2, n};
    return x.shape() == nn and y.shape() == nn and gamma.shape() == nn and xv.extent(0) == n and yv.extent(0) == n and
           inv.N == static_cast<std::size_t>(n) and t_exp.shape() == two_n and t_lin.shape() == two_n and
           cnst.shape() == two_n and t_const.shape() == two_n;
  }
};
}  // namespace polradtran::rt3
