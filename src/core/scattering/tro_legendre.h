#pragma once

#include <matpack.h>

/** Conversions between gridded and Legendre (spectral) scattering matrices of totally randomly oriented particles.
 *
 * A TRO scattering matrix depends on the scattering angle Theta only.  Its
 * spectral form holds, for every element F of [F11, F12, F22, F33, F34, F44],
 * the coefficients a_l on the orthonormal Y_l0 = sqrt((2 l + 1) / 4 pi) P_l(cos(Theta)),
 *
 *   a_l = 2 pi int F(Theta) Y_l0(Theta) dcos(Theta),   F = sum_l a_l Y_l0,
 *
 * so that 2 pi int F11 dcos(Theta) = sqrt(4 pi) a_0, and the Legendre
 * coefficients b_l of F = sum_l b_l P_l are b_l = a_l sqrt((2 l + 1) / 4 pi).
 *
 * Gridded data define F between their nodes as ARTS's gridded paths read it:
 * linear in the scattering angle between nodes, and constant beyond the first
 * and the last node.  project() gives the coefficients of exactly that
 * function, and evaluate() the series at any angle, so neither step
 * interpolates or samples beyond what the data themselves define.
 */
namespace scattering::tro_legendre {
/** The coefficients on Y_l0, l = 0..degree, of the gridded function described above
 *
 * @param angles The scattering angles [deg], ascending in [0, 180]
 * @param values [angles.size(), n], n functions of the scattering angle
 * @param degree The highest degree, >= 0
 * @return [degree + 1, n]
 */
Matrix project(const ConstVectorView& angles, const ConstMatrixView& values, Index degree);

/** The series sum_l a_l Y_l0 at the scattering angles [deg]
 *
 * @param coefficients [degree + 1, n]
 * @param angles The scattering angles [deg], in [0, 180]
 * @return [angles.size(), n]
 */
Matrix evaluate(const ConstMatrixView& coefficients, const ConstVectorView& angles);

/** How well a truncated series represents gridded data, for one set of six TRO elements
 *
 * Relative errors are relative to the largest |F11| at the nodes.
 */
struct report {
  //! max over the nodes of |series - data| per element, relative
  Vector reconstruction_error = Vector(6, 0.0);
  //! |a_degree| per element relative to |a_0| of F11: how far the series is from converged
  Vector tail = Vector(6, 0.0);
  //! The smallest F11 of the series on a fine grid, relative; negative values are truncation ringing
  Numeric min_f11{0.0};
  //! The asymmetry parameter <cos(Theta)> of the series, a_1 / (sqrt(3) a_0) of F11
  Numeric asymmetry{0.0};
};

/** The report on coefficients [degree + 1, 6] of gridded data [angles.size(), 6] at angles [deg] */
report assess(const ConstMatrixView& coefficients, const ConstVectorView& angles, const ConstMatrixView& values);
}  // namespace scattering::tro_legendre

namespace scattering {
/** How well a Legendre series represents gridded TRO data, per temperature and frequency (see tro_legendre::report)
 *
 * Errors are relative to the largest |F11| at the nodes of each temperature
 * and frequency.
 */
struct LegendreReport {
  //! [t, f, 6], max over the nodes of |series - data| per element
  Tensor3 reconstruction_error;
  //! [t, f, 6], |a_degree| per element relative to |a_0| of F11
  Tensor3 tail;
  //! [t, f], the smallest F11 of the series on a fine grid; negative values are truncation ringing
  Matrix min_f11;
  //! [t, f], the asymmetry parameter of the series
  Matrix asymmetry;
  //! [t, f], (2 pi int F11 dcos(Theta) - (K11 - a1)) / K11 of the series, NaN where the optics are not known
  Matrix normalisation_error;

  LegendreReport() = default;
  LegendreReport(Index n_temps, Index n_freqs);

  void set(Index i_t, Index i_f, const tro_legendre::report& r);
};
}  // namespace scattering
