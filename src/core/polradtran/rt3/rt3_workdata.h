#pragma once

#include <matpack.h>
#include <polradtran_workdata.h>

#include <array>

#include "rt3_fft.h"

namespace polradtran::rt3 {
/** The work arrays of rt3::radtran and the routines it calls, for
 * n = nstokes * nummu streams, aziorder + 1 azimuth modes, num_layers
 * layers and nsl scattering sets: those it shares with RT4
 * (polradtran::workdata) and its own.
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
struct rt3_workdata : workdata {
  /////////////////////////////////////////////////////////////////////////////
  // The scattering sets and the layers' optics: made once per call, before
  // the azimuth modes.  scatbuf[s] and directbuf[s] are the parts of SCATBUF
  // and DIRECTBUF of set s + 1, as SCATTERING and DIRECT_SCATTERING lay
  // them out; legendre_coef is GET_SCAT_SET's COEF, one set at a time.
  /////////////////////////////////////////////////////////////////////////////

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
  // RADTRAN on the streams in a mode, beyond what it shares with RT4: the
  // set's scattering matrix (SCATTER_MATRIX) and direct vector, the thermal
  // vector, the solar source of the initial sublayer of a scattering layer
  // (exp_source), and the radiance incident from the ground and the part of
  // it reflected from the direct beam.
  /////////////////////////////////////////////////////////////////////////////

  Tensor5 scatter_matrix;                    //!< [4, nummu, nstokes, nummu, nstokes]
  Matrix  direct_vector, thermal_vector;     //!< [2, n]
  Matrix  exp_source;                        //!< [2, n]
  Matrix  ground_radiance, direct_radiance;  //!< [nummu, nstokes]

  /////////////////////////////////////////////////////////////////////////////
  // DOUBLING_INTEGRATION's solar source, T_EXP; the scratch of its linear
  // source is shared.
  /////////////////////////////////////////////////////////////////////////////

  Matrix t_exp;  //!< [2, n]

  rt3_workdata() = default;
  rt3_workdata(Index nstokes, Index nummu, Index aziorder, Index num_layers, Index nsl, Index legendre_rows) {
    resize(nstokes, nummu, aziorder, num_layers, nsl, legendre_rows);
  }

  //! Sizes the arrays but the scratch of the scattering routines; this
  //! allocates only where an array grows
  void resize(Index nstokes, Index nummu, Index aziorder, Index num_layers, Index nsl, Index legendre_rows) {
    const Index n = nstokes * nummu;
    workdata::resize(nstokes, nummu, num_layers);

    legendre_coef.resize(legendre_rows, 6);
    set_extinct.resize(nsl);
    set_scatter.resize(nsl);
    scatbuf.resize(nsl, aziorder + 1, 2, nummu, nummu, nstokes, nstokes);
    directbuf.resize(nsl, aziorder + 1, 2, nummu, nstokes);
    scat_nums.resize(num_layers);
    extinctions.resize(num_layers);
    albedos.resize(num_layers);
    direct_level_flux.resize(num_layers + 1);

    scatter_matrix.resize(4, nummu, nstokes, nummu, nstokes);
    for (Matrix* m : {&direct_vector, &thermal_vector, &exp_source, &t_exp}) m->resize(2, n);
    ground_radiance.resize(nummu, nstokes);
    direct_radiance.resize(nummu, nstokes);
  }

  //! Whether the scratch of DOUBLING_INTEGRATION, COMBINE_LAYERS and
  //! INTERNAL_RADIANCE (the shared one and T_EXP) is sized for n streams
  [[nodiscard]] bool scratch_sized(Index n) const {
    return workdata::scratch_sized(n) and t_exp.shape() == (std::array<Index, 2>{2, n});
  }
};
}  // namespace polradtran::rt3
