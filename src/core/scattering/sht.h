
/* Copyright (C) 2023
   Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*****************************************************************************
 ***  File description
 *****************************************************************************/

/*!
   \file   sht.h
   \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
   \date 2023-02-03

   An interface to the SHTns library
*/
#ifndef ARTS_CORE_SCATTERING_SHT_H_
#define ARTS_CORE_SCATTERING_SHT_H_

#include <array.h>
#include <arts_conversions.h>
#include <matpack.h>
#include <scattering/integration.h>

#include <complex>
#include <map>
#include <memory>
#include <numbers>
#include <utility>

typedef struct shtns_info *shtns_cfg;

namespace scattering {
namespace sht {

using std::numbers::pi_v;

/** FFTW coefficient array
 *
 * SHTns library requires memory to be aligned according to fftw3
 * requirements which is ensured by using fftw_malloc and fftw_free
 * functions to allocate and free memory, respectively.
 */
template <typename Numeric>
class FFTWArray {
 public:
  FFTWArray() = default;
  FFTWArray(Index n);
  [[nodiscard]] operator Numeric *() const { return ptr_.get(); }

 private:
  std::shared_ptr<Numeric> ptr_ = nullptr;
};

class ShtnsHandle {
 public:
  [[nodiscard]] static shtns_cfg get(Index l_max,
                                     Index m_max,
                                     Index n_aa,
                                     Index n_za);

 private:
  static std::array<Index, 4> current_config_;
  static shtns_cfg shtns_;
};

////////////////////////////////////////////////////////////////////////////////
// SHT
////////////////////////////////////////////////////////////////////////////////
/** Spherical harmonics transformation
 *
 * Represents a spherical harmonics transformation (SHT) for fixed spatial- and
 * spectral grid sizes. A SHT object acts as a wrapper around the
 * SHTns library. Each object holds the arrays required to store
 * spatial and spectral coefficients.
 *
 * A specific SHT configuration is defined by its configuration which is
 * an array containing the numbers (l_max, m_max, n_aa, n_za), which
 * mean the following:
 *
 * - l_max: The degree of the spherical harmonics transform.
 * - m_max: The order of the spherical harmonics transform
 * - n_aa: The number of points in the azimuth-angle grid. Must satisfy
 *     n_aa > 2 * m_max + 1
 * - n_za: The number of points in the zenith-angle grid. Must satisfy
 *     n_za > 2 * l_max + 1
 *
 * Note that the underlying shtns library only support one active configuration
 * at a time. Performing transforms with to distinct SHT objects therefore requires
 * reinitialization of shtns. Efficient use of the interface should therefore try
 * to minimize switching between different SHT configurations.
 *
 * NOTE: This interface is not thread safe.
 */
class SHT {
 public:
  /** Return azimuth-angle grid used by the SH transform.
   * @param The number of points in the zenith-angle grid.
   * @param radians If true the zenith-angle grid is returned in radians.
   * @return A vector containing the azimuth angles in degree (or radians).
   */
  [[nodiscard]] static Vector get_azimuth_angle_grid(Index n_aa,
                                                     bool radians = false);

  /// Pointer to the azimuth angle grid in degrees.
  [[nodiscard]] std::shared_ptr<const Vector> get_aa_grid_ptr() const {
    return aa_grid_;
  }

  /** Return zenith-angle grid used by the SH transform.
   * @param The number of points in the zenith-angle grid.
   * @param radians If true the zenith-angle grid is returned in radians.
   * @return A vector containing the zenith-angle grid in degree (or radians).
   */
  [[nodiscard]] static FejerGrid get_zenith_angle_grid(Index n_pts,
                                                       bool radians = false);

  [[nodiscard]] FejerGrid get_zenith_angle_grid(bool radians = false) const;

  /// Pointer to the zenith angle grid in degrees.
  [[nodiscard]] std::shared_ptr<const ZenithAngleGrid> get_za_grid_ptr() const;

  /** Calculates the number of spherical harmonics coefficients for a real
   * transform.
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @return The number of spherical harmonics coefficients.
   */
  [[nodiscard]] static Index calc_n_spectral_coeffs(Index l_max, Index m_max);

  /** Calculates the number of spherical harmonics coefficients for a complex
   * transform.
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @return The number of spherical harmonics coefficients for a complex
   * transform.
   */
  [[nodiscard]] static Index calc_n_spectral_coeffs_cmplx(Index l_max,
                                                          Index m_max);

  /** Calc l_max for m_max == l_max.
   *
   * Calculates the value of l_max that yields the given number of spectral
   * coefficients under the assumption that m_max is equal to l_max.
   *
   * @param n_spectral_coeffs The number of spectral coefficients.
   * @return l_max value yielding the given number of spectral coeffs.
   */
  [[nodiscard]] static Index calc_l_max(Index n_spectral_coeffs);

  /** SHT parameters for a spatial field of given size.
   * @param n_aa: The size of the azimuth-angle grid.
   * @param n_za: The size of the zenith-angle  grid.
   * @return A 4-element array containing parameters l_max,
   * m_max, n_aa and n_za that are valid inputs to initialize
   * the SHTns library.
   */
  [[nodiscard]] static std::array<Index, 4> get_config_lonlat(Index n_aa,
                                                              Index n_za);

  /** SHT parameters for given SHT maximum degree and order of SHT
   * @param l_max: The maximum degree l of the SHT.
   * @param m_lat: The maximum degree m of the SHT.
   * @return A 4-element array containing parameters l_max,
   * m_max, n_aa and n_za that are valid inputs to initialize
   * the SHTns library.
   */
  [[nodiscard]] static std::array<Index, 4> get_config_lm(Index l_max,
                                                          Index m_max);

  /**
   * Create a spherical harmonics transformation object.
   *
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @param n_aa The number of azimuth-angle grid points.
   * @param n_za The number of zenith-angle grid points.
   */
  SHT(Index l_max, Index m_max, Index n_aa, Index n_za);

  /**
   * Create a spherical harmonics transformation object.
   *
   * The values for n_aa and n_za are set to 2 * l_max + 2 and
   * 2 * m_max + 2, respectively.
   *
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   */
  SHT(Index l_max, Index m_max);

  /**
   * Create a spherical harmonics transformation object.
   *
   * Create spherical harmonics transformation object with l_max == m_max
   * and values for n_aa and n_za set to 2 * l_max + 2 and
   * 2 * m_max + 2, respectively.
   * @param l_max The maximum degree of the SHT.
   */
  SHT(Index l_max);

  /// Serialize SHT.
  std::ostream &serialize(std::ostream &output) const;

  /// Deserialize SHT.
  [[nodiscard]] static SHT deserialize(std::istream &input);

  /** Return the cosine of the zenith-angle grid used by SHTns.
   * @return A vector containing the cosine of the zenith-angle grid.
   */
  [[nodiscard]] Vector get_cos_za_grid();

  /** Return azimuth-angle grid used by SHTns.
   * @return A vector containing the azimuth grid in radians.
   */
  [[nodiscard]] Vector get_azimuth_angle_grid(bool radians = false);

  /** L-indices of the SHT modes.
   *
   * @return A vector of indices containing the l-value corresponding to each
   * element in a spectral coefficient vector.
   */
  ArrayOfIndex get_l_indices();

  /** M-indices of the SHT modes.
   *
   * @return A vector of indices containing the m-value corresponding to each
   * element in a spectral coefficient vector.
   */
  ArrayOfIndex get_m_indices();

  /** Return index of coefficient in coefficient vector.
   *
   * @param l The l parameter of SH corresponding to the coefficient.
   * @param m The m parameter of the SH corresponding to the coefficient.
   * @return The index of the coefficient in the spectral coefficient
   * vector.
   */
  Index get_coeff_index(Index l, Index m);

  /**
   * Copy spatial field into the array that holds spatial data for
   * spherical harmonics computations.
   *
   * @param m Matrix or comparable providing read-only access
   * to the input data. Row indices should correspond to azimuth angles
   * and columns to zenith angle.
   */
  template <typename T>
  void set_spatial_coeffs(
      const matpack::strided_view_t<const T, 2> &view) const {
    assert(view.nrows() == n_aa_);
    assert(view.ncols() == n_za_);
    Index index = 0;
    for (int i = 0; i < view.nrows(); ++i) {
      for (int j = 0; j < view.ncols(); ++j) {
        if constexpr (matpack::complex_type<T>) {
          spatial_coeffs_cmplx_[index] = view[i, j];
        } else {
          spatial_coeffs_[index] = view[i, j];
        }
        ++index;
      }
    }
  }

  /**
   * Copy complex spatial field into the array that holds spatial data for
   * real spherical harmonics computations.
   *
   * @param m Matrix or comparable providing read-only access
   * to the input data. Row indices should correspond to azimuth angles
   * and columns to zenith angles.
   */
  void set_spectral_coeffs(
      const matpack::strided_view_t<const Complex, 1> &view) const;

  /**
   * Copy spherical harmonics coefficients into the array that holds spectral
   * data for complex spherical harmonics computations.
   *
   * @param m SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   */
  void set_spectral_coeffs_cmplx(
      const matpack::strided_view_t<const Complex, 1> &view) const;

  /**
   * Return content of the array that holds spatial data for
   * spherical harmonics computations.
   *
   * @return  matrix containing the spatial field. Row indices should
   * correspond to azimuth angles  and columns to zenith angles.
   * angle).
   */
  StridedConstMatrixView get_spatial_coeffs() const;

  /**
   * Return content of the array that holds complex spatial data for
   * spherical harmonics computations.
   *
   * @return  matrix containing the complex spatial field. Row indices
   * should correspond to azimuth angles and columns to zenith angles.
   *
   */
  StridedConstComplexMatrixView get_spatial_coeffs_cmplx() const;

  /**
   * @return The size of the zenith-angle grid.
   */
  Index get_n_zenith_angles() const;
  /**
   * @return The size of the azimuth angle grid.
   */
  Index get_n_azimuth_angles() const;
  /**
   * @return The number of spherical harmonics coefficients.
   */
  Index get_n_spectral_coeffs() const;
  /**
   * @return The number of spherical harmonics coefficients.
   */
  Index get_n_spectral_coeffs_cmplx() const;

  /**
   * @return The maximum degree l of the SHT transformation.
   */
  Index get_l_max() const;

  /**
   * @return The maximum order m of the SHT transformation.
   */
  Index get_m_max() const;

  /**
   * Return content of the array that holds spectral data for
   * spherical harmonics computations.
   *
   * @return m  vector containing the spherical harmonics coefficients
   * representing the data.
   */
  StridedConstComplexVectorView get_spectral_coeffs() const;

  /**
   * Return content of the array that holds spectral data for
   * spherical harmonics computations of complex fields.
   *
   * @return m Vector containing the spherical harmonics coefficients
   * representing the data.
   */
  StridedConstComplexVectorView get_spectral_coeffs_cmplx() const;

  /** Apply forward SHT Transform *
   * Transforms discrete spherical data into spherical harmonics representation.
   * @param view GridCoeffs containing the data. Row indices should correspond to
   * azimuth angles and columns to zenith angles.
   * @return Coefficient vector containing the spherical harmonics coefficients.
   */
  ComplexVector transform(const StridedConstMatrixView &view);

  /** Apply forward SHT Transform
   *
   * Transforms discrete spherical data into spherical harmonics representation.
   * @param view GridCoeffs containing the data. Row indices should correspond to
   * azimuth angles and columns to zenith angles.
   * @return Coefficient vector containing the spherical harmonics coefficients.
   */
  ComplexVector transform_cmplx(const StridedConstComplexMatrixView &view);

  /** Apply inverse SHT Transform
   *
   * Transforms discrete spherical data given in spherical harmonics
   * representation back to spatial domain.
   *
   * @param view SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   * @return GridCoeffs containing the spatial data.
   */
  Matrix synthesize(const StridedConstComplexVectorView &view);

  /** Apply inverse SHT Transform for complex data.
   *
   * Transforms discrete spherical data given in spherical harmonics
   * representation back to spatial domain.
   *
   * @param view SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   * @return GridCoeffs containing the spatial data.
   */
  ComplexMatrix synthesize_cmplx(const StridedConstComplexVectorView &view);

  /** Evaluate spectral representation at given point.
   *
   * @param view Spectral coefficient vector containing the SH coefficients.
   * @param phi The azimuth angles in radians.
   * @return theta The zenith angle in radians.
   */
  Numeric evaluate(const StridedConstComplexVectorView &view,
                   Numeric phi,
                   Numeric theta);

  /** Evaluate spectral representation at given point.
   *
   * @param view Spectral coefficient vector containing the SH coefficients.
   * @param points 2-row matrix containing the points (azimuth angle, zenith ange)
   * at which to evaluate the function.
   * @return A vector containing the values corresponding to the points
   * in points.
   */
  Vector evaluate(const StridedComplexVectorView &view, const StridedMatrixView &points);

  /** Evaluate 1D spectral representation at given point.
   *
   * This method covers the special case of 1D data that varies
   * only along zenith angles. In this case the SH transform degenerates
   * to a Legendre transform.
   *
   * @param view Spectral coefficient vector containing the SH coefficients.
   * @param Vector containing the zenith angles within [0, PI] to evaluate the
   * function.
   * @return A vector containing the values corresponding to the points
   * in points.
   */
  Vector evaluate(const StridedConstComplexVectorView &view, const Vector &thetas);

  template <typename Vec1, typename Vec2>
  friend matpack::data_t<typename Vec1::value_type, 1> add_coeffs(
      const SHT &sht_v, const Vec1 &v, const SHT &sht_w, const Vec2 &w);
  template <typename T>
  friend matpack::data_t<T, 2> add_coeffs(
      const SHT &sht_inc_v,
      const SHT &sht_scat_v,
      const matpack::strided_view_t<const T, 2> &v,
      const SHT &sht_inc_w,
      const SHT &sht_scat_w,
      const matpack::strided_view_t<const T, 2> &w);

 private:
  bool is_trivial_;
  Index l_max_, m_max_, n_aa_, n_za_, n_spectral_coeffs_,
      n_spectral_coeffs_cmplx_;

  sht::FFTWArray<std::complex<double>> spectral_coeffs_, spectral_coeffs_cmplx_,
      spatial_coeffs_cmplx_;
  sht::FFTWArray<double> spatial_coeffs_;
  std::shared_ptr<Vector> aa_grid_;
  std::shared_ptr<ZenithAngleGrid> za_grid_;
};

/** SHT instance provider.
 *
 * Simple cache for SHT instances.
 */
class SHTProvider {
 public:
  using SHTParams = std::array<Index, 4>;

  SHTProvider() {}

  /** Get SHT instance for given SHT parameters.
   * @arg params Length-4 array containing the parameters required to initialize
   * the SHT transform: l_max, m_max, n_aa, n_za. See documention of SHT class
   * for explanation of their significance.
   * @return shared pointer to SHT instance.
   */
  std::shared_ptr<SHT> get_instance(SHTParams params [[maybe_unused]]);

  std::shared_ptr<SHT> get_instance(Index n_aa, Index n_za);

  std::shared_ptr<SHT> get_instance_lm(Index l_max, Index m_max);

 protected:
  std::map<SHTParams, std::shared_ptr<SHT>> sht_instances_;
};

extern SHTProvider provider;

/** Addition of spectral coefficients
 *
 * Add spectral coefficients  from potentially different
 * SHT transforms. The coefficients of the sum are computed with
 * respect to the SHT of the first summand.
 *
 * @param sht_v: The SHT object used to derive the spectral coefficients
 * in v
 * @param v: The first summand
 * @param sht_w: The SHT object used to derive the spectral coefficients
 * in w
 * @param w: The second summand
 * @return A spectral coefficient field containing spectral coefficients
 * of the sum of v and w for the transform object sht_v.
 */
template <typename Vec1, typename Vec2>
matpack::data_t<typename Vec1::value_type, 1> add_coeffs(const SHT &sht_v,
                                                               const Vec1 &v,
                                                               const SHT &sht_w,
                                                               const Vec2 &w) {
  auto result = matpack::data_t<typename Vec1::value_type, 1>(v);

  if (sht_w.is_trivial_) {
    result[0] += w[0];
    return result;
  }

  Index m_max_min = std::min(sht_v.m_max_, sht_w.m_max_);
  Index l_max_min = std::min(sht_v.l_max_, sht_w.l_max_);
  for (Index m = 0; m <= m_max_min; ++m) {
    Index index_r = m * (sht_w.l_max_ + 1) - (m * (m - 1)) / 2;
    Index index_l = m * (sht_v.l_max_ + 1) - (m * (m - 1)) / 2;
    for (Index l = m; l <= l_max_min; ++l) {
      result[index_l] += w[index_r];
      ++index_r;
      ++index_l;
    }
  }
  return result;
}

/** In-place addition of doubly-transformed spectral coefficients.
  *
  * Add spectral coefficient matrices derived from potentially different
  * SHT transforms. The coefficients of the sum are computed with
  * respect to the SHT of the first summand.
  *
  * @param sht_inc_v: The SHT object used to transform the incoming
  * angles of v.
  * @param sht_scat_v: The SHT object used to transform the scattering
  * angles of v.
  * @param v: The first summand
  * @param sht_inc_w: The SHT object used to transform the incoming
  * angles of w.
  * @param sht_scat_w: The SHT object used to transform the scattering
  * angles of v.
  * @param w: The second summand
  * @return A spectral coefficient field containing spectral coefficients
  * of the sum of v and w for the transform objects sht_inc_v and
  * sht_scat_v.
  */
template <typename T>
matpack::data_t<T, 2> add_coeffs(
    const SHT &sht_inc_v,
    const SHT &sht_scat_v,
    const matpack::strided_view_t<const T, 2> &v,
    const SHT &sht_inc_w,
    const SHT &sht_scat_w,
    const matpack::strided_view_t<const T, 2> &w) {
  Index nlm_inc  = sht_inc_v.get_n_spectral_coeffs_cmplx();
  Index nlm_scat = sht_scat_v.get_n_spectral_coeffs();
  auto result    = ComplexMatrix(nlm_inc, nlm_scat);

  Index index_l = 0;
  for (int l = 0; l <= static_cast<int>(sht_inc_v.l_max_); ++l) {
    int m_max = static_cast<int>(
        (l <= static_cast<int>(sht_inc_v.m_max_)) ? l : sht_inc_v.m_max_);
    for (int m = -m_max; m <= m_max; ++m) {
      if ((l > sht_inc_w.l_max_) || (std::abs(m) > sht_inc_w.m_max_)) {
        result[index_l, joker] = v.row(index_l, joker);
      } else {
        int h = static_cast<int>(
            std::min<int>(static_cast<int>(sht_inc_w.m_max_), l));
        int index_r = l * (2 * h + 1) - h * h + m;

        auto r =
            add_coeffs(sht_scat_v, v.row(index_l), sht_scat_w, w.row(index_r));
        result[index_l, joker] = r;
      }
      ++index_l;
    }
  }
  return result;
}

}  // namespace sht

}  // namespace scattering

#endif
