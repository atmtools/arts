/* Copyright (C) 2020-2022
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
   \date 2022-10-14

   Provides an interface to the SHTns library.
*/
#ifndef __SCATTERING_SHT__
#define __SCATTERING_SHT__

#include <fftw3.h>
#include <shtns.h>

#include <complex>
#include <iostream>
#include <map>
#include <memory>

#include "Eigen/Dense"

#include <scattering/math.h>
#include <scattering/integration.h>

using GridCoeffs = scattering::math::Matrix<double>;
using GridCoeffsRef = scattering::math::ConstMatrixRef<double>;
using CmplxGridCoeffs = scattering::math::Matrix<std::complex<double>>;
using CmplxGridCoeffsRef = scattering::math::ConstMatrixRef<std::complex<double>>;
using SpectralCoeffs = scattering::math::Vector<std::complex<double>>;
using SpectralCoeffsRef = scattering::math::ConstVectorRef<std::complex<double>>;
using SpectralCoeffMatrix = scattering::math::Matrix<std::complex<double>>;
using SpectralCoeffMatrixRef =
    scattering::math::ConstMatrixRef<std::complex<double>>;

namespace scattering {
namespace sht {

using scattering::math::Index;


/** FFTW coefficient array
 *
 * SHTns library requires memory to be aligned according to fftw3
 * requirements which is ensured by using fftw_malloc and fftw_free
 * functions to allocate and free memory, respectively.
 */
template <typename Numeric>
class FFTWArray {
 public:
  FFTWArray() {}
  FFTWArray(Index n);
  operator Numeric *() const { return ptr_.get(); }

 private:
  std::shared_ptr<Numeric> ptr_ = nullptr;
};

class ShtnsHandle {
 public:
 static shtns_cfg get(Index l_max, Index m_max, Index n_lon, Index n_lat);

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
 * spectral grid sizes. A SHT object essentially acts as a wrapper around the
 * SHTns library. Each object holds a SHTns configuration as well as the arrays
 * required to store spatial and spectral coefficients.
 */
// pxx :: export
class SHT {
 public:
  using Vector = scattering::math::Vector<double>;
  using LatGrid = QuadratureLatitudeGrid<FejerQuadrature<double>, double>;
  using IndexVector = scattering::math::Vector<Index>;
  using ConstVectorMap = scattering::math::ConstVectorMap<double>;

  /** Addition of spectral coefficients
   *
   * Add spectral coefficients derived from potentially different
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
  static SpectralCoeffs add_coeffs(const SHT &sht_v,
                                   SpectralCoeffsRef v,
                                   const SHT &sht_w,
                                   SpectralCoeffsRef w);

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
  static SpectralCoeffMatrix add_coeffs(const SHT &sht_inc_v,
                                        const SHT &sht_scat_v,
                                        SpectralCoeffMatrixRef v,
                                        const SHT &sht_inc_w,
                                        const SHT &sht_scat_w,
                                        SpectralCoeffMatrixRef w);

  /// The latitude grid (zenith angles) of the spatial grid.
  static LatGrid get_latitude_grid(Index n_lat);
  /// The longitude grid (azimuth angles) of the spatial grid.
  static Vector get_longitude_grid(Index n_lon);

  /** Calculates the number of spherical harmonics coefficients for a real
   * transform.
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @return The number of spherical harmonics coefficients.
   */
  static Index calc_n_spectral_coeffs(Index l_max, Index m_max) {
    return (l_max + 1) * (m_max + 1) - (m_max * (m_max + 1)) / 2;
  }

  /** Calculates the number of spherical harmonics coefficients for a complex
   * transform.
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @return The number of spherical harmonics coefficients for a complex
   * transform.
   */
  static Index calc_n_spectral_coeffs_cmplx(Index l_max, Index m_max) {
    return (2 * m_max + 1) * (l_max + 1) - m_max * (m_max + 1);
  }

  /** Calc l_max for m_max == l_max.
   *
   * Calculates the value of l_max that yields the given number of spectral
   * coefficients under the assumption that m_max is equal to l_max.
   *
   * @param n_spectral_coeffs The number of spectral coefficients.
   * @return l_max value yielding the given number of spectral coeffs.
   */
  static Index calc_l_max(Index n_spectral_coeffs) {
    return static_cast<Index>(sqrt(2.0 * static_cast<double>(n_spectral_coeffs) + 0.25) - 1.5);
  }

  /** SHT parameters for a spatial field of given size.
   * @param n_lon: The size of the longitude (azimuth) grid.
   * @param n_lat: The size of the latitude (zenith) grid.
   * @return A 4-element array containing parameters l_max,
   * m_max, n_lon and n_lat that are valid inputs to initialize
   * the SHTns library.
   */
  static std::array<Index, 4> get_params(Index n_lon, Index n_lat);

  /**
   * Create a spherical harmonics transformation object.
   *
   * @param l_max The maximum degree of the SHT.
   * @param m_max The maximum order of the SHT.
   * @param n_lon The number of longitude grid points.
   * @param n_lat The number of co-latitude grid points.
   */
  SHT(Index l_max, Index m_max, Index n_lon, Index n_lat);

  /**
   * Create a spherical harmonics transformation object.
   *
   * The values for n_lon and n_lat are set to 2 * l_max + 2 and
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
   * and values for n_lon and n_lat set to 2 * l_max + 2 and
   * 2 * m_max + 2, respectively.
   * @param l_max The maximum degree of the SHT.
   */
  SHT(Index l_max);

  /** Return latitude grid used by SHTns.
   * @return Eigen vector containing the latitude grid in radians.
   */
  LatGrid get_latitude_grid();

  /** Return the cosine of the latitude grid used by SHTns.
   * @return Eigen vector containing the co-latitude grid.
   */
  Vector get_colatitude_grid();

  /** Return longitude grid used by SHTns.
   * @return Eigen vector containing the longitude grid in radians.
   */
  Vector get_longitude_grid();

  /** L-indices of the SHT modes.
   *
   * @return A vector of indices containing the l-value corresponding to each
   * element in a spectral coefficient vector.
   */
  IndexVector get_l_indices();


  /** M-indices of the SHT modes.
   *
   * @return A vector of indices containing the m-value corresponding to each
   * element in a spectral coefficient vector.
   */
  IndexVector get_m_indices();

  /** Return index of coefficient in coefficient vector.
   *
   * @param l The l parameter of SH corresponding to the coefficient.
   * @param m The m parameter of the SH corresponding to the coefficient.
   * @return The index of the coefficient in the spectral coefficient
   * vector.
   */
  Index get_coeff_index(int l, int m);

  /**
   * Copy spatial field into the array that holds spatial data for
   * spherical harmonics computations.
   *
   * @param m Eigen::Matrix or comparable providing read-only access
   * to the input data. Row indices should correspond to longitudes
   * (azimuth angle) and columns to latitudes (zenith angle).
   */
  void set_spatial_coeffs(const GridCoeffsRef &m) const;

  /**
   * Copy complex spatial field into the array that holds spatial data for
   * spherical harmonics computations.
   *
   * @param m Eigen::Matrix or comparable providing read-only access
   * to the input data. Row indices should correspond to longitudes
   * (azimuth angle) and columns to latitudes (zenith angle).
   */
  void set_spatial_coeffs(const CmplxGridCoeffsRef &m) const;

  /**
   * Copy spherical harmonics coefficients into the array that holds spectral
   * data for spherical harmonics computations.
   *
   * @param m SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   */
  void set_spectral_coeffs(const SpectralCoeffsRef &m) const;

  /**
   * Copy spherical harmonics coefficients into the array that holds spectral
   * data for spherical harmonics computations involving complex spatial fields.
   *
   * @param m Eigen vector containing the spherical harmonics coefficients
   * representing the data.
   */
  void set_spectral_coeffs_cmplx(const SpectralCoeffsRef &m) const;

  /**
   * Return content of the array that holds spatial data for
   * spherical harmonics computations.
   *
   * @return Eigen matrix containing the spatial field. Row indices should
   * correspond to longitudes (azimuth angle) and columns to latitudes (zenith
   * angle).
   */
  GridCoeffs get_spatial_coeffs() const;

  /**
   * Return content of the array that holds complex spatial data for
   * spherical harmonics computations.
   *
   * @return Eigen matrix containing the complex spatial field. Row indices
   * should correspond to longitudes (azimuth angle) and columns to latitudes
   * (zenith angle).
   */
  CmplxGridCoeffs get_cmplx_spatial_coeffs() const;

  /**
   * @return The size of the co-latitude grid.
   */
  Index get_n_latitudes() const { return n_lat_; }
  /**
   * @return The size of the longitude grid.
   */
  Index get_n_longitudes() const { return n_lon_; }
  /**
   * @return The number of spherical harmonics coefficients.
   */
  Index get_n_spectral_coeffs() const { return n_spectral_coeffs_; }
  /**
   * @return The number of spherical harmonics coefficients.
   */
  Index get_n_spectral_coeffs_cmplx() const { return n_spectral_coeffs_cmplx_; }

  /**
   * @return The maximum degree l of the SHT transformation.
   */
  Index get_l_max() { return l_max_; }

  /**
   * @return The maximum order m of the SHT transformation.
   */
  Index get_m_max() { return m_max_; }

  /**
   * Return content of the array that holds spectral data for
   * spherical harmonics computations.
   *
   * @return m Eigen vector containing the spherical harmonics coefficients
   * representing the data.
   */
  SpectralCoeffs get_spectral_coeffs() const;

  /**
   * Return content of the array that holds spectral data for
   * spherical harmonics computations of complex fields.
   *
   * @return m Eigen vector containing the spherical harmonics coefficients
   * representing the data.
   */
  SpectralCoeffs get_spectral_coeffs_cmplx() const;

  /** Apply forward SHT Transform
   *
   * Transforms discrete spherical data into spherical harmonics representation.
   * @param m GridCoeffs containing the data. Row indices should correspond to
   * longitudes (azimuth angle) and columns to latitudes (zenith angle).
   * @return Coefficient vector containing the spherical harmonics coefficients.
   */
  SpectralCoeffs transform(const GridCoeffsRef &m);

  /** Apply forward SHT Transform
   *
   * Transforms discrete spherical data into spherical harmonics representation.
   * @param m GridCoeffs containing the data. Row indices should correspond to
   * longitudes (azimuth angle) and columns to latitudes (zenith angle).
   * @return Coefficient vector containing the spherical harmonics coefficients.
   */
  SpectralCoeffs transform_cmplx(const CmplxGridCoeffsRef &m);

  /** Apply inverse SHT Transform
   *
   * Transforms discrete spherical data given in spherical harmonics
   * representation back to spatial domain.
   *
   * @param m SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   * @return GridCoeffs containing the spatial data.
   */
  GridCoeffs synthesize(const SpectralCoeffsRef &m);

  /** Apply inverse SHT Transform for complex data.
   *
   * Transforms discrete spherical data given in spherical harmonics
   * representation back to spatial domain.
   *
   * @param m SpectralCoeffs The spherical harmonics coefficients
   * representing the data.
   * @return GridCoeffs containing the spatial data.
   */
  CmplxGridCoeffs synthesize_cmplx(const SpectralCoeffsRef &m);

  /** Evaluate spectral representation at given point.
   *
   * @param m Spectral coefficient vector containing the SH coefficients.
   * @param phi The azimuth angles in radians.
   * @return theta The zenith angle in radians.
   */
  double evaluate(
      const SpectralCoeffsRef &m,
      double phi,
      double theta);

  /** Evaluate spectral representation at given point.
   *
   * @param m Spectral coefficient vector containing the SH coefficients.
   * @param points 2-row matrix containing the points (lon, lat) at which
   * to evaluate the function.
   * @return A vector containing the values corresponding to the points
   * in points.
   */
  Vector evaluate(
      const SpectralCoeffsRef &m,
      const math::MatrixFixedRows<double, 2> &points);

  /** Evaluate 1D spectral representation at given point.
   *
   * This method covers the special case of 1D data that varies
   * only along latitudes. In this case the SH transform degenerates
   * to a Legendre transform.
   *
   * @param m Spectral coefficient vector containing the SH coefficients.
   * @param Vector containing the latitudes within [0, PI] to evaluate the
   * function.
   * @return A vector containing the values corresponding to the points
   * in points.
   */
  Vector evaluate(const SpectralCoeffsRef &m,
                                 const Vector &thetas);

 private:
  bool is_trivial_;
  Index l_max_, m_max_, n_lon_, n_lat_, n_spectral_coeffs_,
      n_spectral_coeffs_cmplx_;

  sht::FFTWArray<std::complex<double>> spectral_coeffs_, spectral_coeffs_cmplx_,
      cmplx_spatial_coeffs_;
  sht::FFTWArray<double> spatial_coeffs_;
};

/** SHT instance provider.
 *
 * Simple cache that caches created SHT instances.
 */
class SHTProvider {
 public:
  using SHTParams = std::array<Index, 4>;

  SHTProvider();

  /** Get SHT instance for given SHT parameters.
   * @arg params Length-4 array containing the parameters required to initialize
   * the SHT transform: l_max, m_max, n_lat, n_lon. See documention of SHT class
   * for explanation of their significance.
   * @return Reference to SHT instance.
   */
  SHT &get_sht_instance(SHTParams params);

 protected:
  std::map<SHTParams, std::unique_ptr<SHT>> sht_instances_;
};

}  // namespace sht
}  // namespace scattering

#endif
