#ifndef SCATTERING_MIE_H_
#define SCATTERING_MIE_H_

/** \file mie.h
 *
 * This file provides the MieSphere class for the calculation of scattering
 * properties of Mie spheres.
 *
 * The implementation is based on:
 *
 * [1] Bohren, C. F., & Huffman, D. R. (2008). Absorption and scattering of
 * light by small particles. John Wiley & Sons.
 *
 * [2] Wiscombe, W. J. (1979). Mie scattering calculations: Advances in
 * technique and fast, vector-speed computer codes (Vol. 10). National Technical
 * Information Service, US Department of Commerce.
 *
 * Calculation of refractive indices is taken from:
 *
 * [3] W. J. Ellison , "Permittivity of Pure Water, at Standard Atmospheric
 * Pressure, over the Frequency Range 0–25THz and the Temperature Range
 * 0–100°C", Journal of Physical and Chemical Reference Data 36, 1-18
 * (2007) https://doi.org/10.1063/1.2360986
 *
 * [4] Mätzler, Christian, et al. "Microwave dielectric properties of
 *  ice." Thermal microwave radiation: applications for remote sensing
 *  52 (2006): 455-462.
 *
 * Author: Simon Pfreundschuh, 2022
 */
#pragma once

#include <matpack.h>
#include <rtepack.h>

namespace scattering {

/** Calculates the refractive index of water following [3].
 *
 * The implementation is essentially copied from refraction.cc.
 *
 * @param frequency The frequency in Hz
 * @param temperature The temperature in K.
 */
Complex refr_index_water_ellison07(Numeric frequency, Numeric temperature);

/** Calculates the refractive index of ice following [4].
 *
 * The implementation is copied from refraction.cc.
 *
 * @param frequency The frequency in Hz
 * @param temperature The temperature in K.
 */
Complex refr_index_ice_matzler06(Numeric frequency, Numeric temperature);

/** Calculation of the log-derivative using downward recursion.
 *
 * The logarithmic derivative is defined by Eq. (4.89) in [1]. The suggested
 * number of steps in the recurrence is:
 *
 *  n_steps = max(x + 4 * x^(1/3) + 2, abs(x * n)) + 15
 *
 * @param rho: The product of relative refractive index n and size parameter x.
 * @param n_steps: The number of steps in the downward recursion.
 * @return A vector containing the 'n_steps' calculated complex values
 * of the logarithmic derivative
 */
ComplexVector log_derivative(Complex rho, Size n_steps);

/** A Mie sphere
 *
 * The MieSphere class calculates the Mie parameters of a given sphere and
 * provides access to to its scattering properties.
 */
class MieSphere {
 public:
  /** Create a Mie sphere.
     *
     * @param lambda The wavelength of the radiation in meters.
     * @param radius The radius of the particle in meters.
     * @param n The relative refractive index of the sphere
     * @param theta Vector of scattering angles in degree for which to compute
     * the scattering matrix.
     */
  MieSphere(Numeric lambda, Numeric radius, Complex n, StridedConstVectorView theta);

  /** Calculate scattering properties of a spherical water particle.
     *
     * @param frequency The frequency in Hz.
     * @param temperature The ambient temperature in K.
     * @param radius The radius of the sphere.
     * @param theta The scattering angles for which to calculate the
     *     the scattering properties.
     * @return A MieSphere object representing the scattering properties
     * of the corresponding spherical water droplet.
     */
  static MieSphere Liquid(Numeric frequency, Numeric temperature, Numeric radius, StridedConstVectorView theta);

  /** Calculate scattering properties of a spherical ice particle.
     *
     * @param frequency The frequency in Hz.
     * @param temperature The ambient temperature in K.
     * @param radius The radius of the sphere.
     * @param theta The scattering angles for which to calculate the
     *     the scattering properties.
     * @return A MieSphere object representing the scattering properties
     * of the corresponding spherical ice particle.
     */
  static MieSphere Ice(Numeric frequency, Numeric temperature, Numeric radius, StridedConstVectorView theta);

  /// The scattering coefficient.
  [[nodiscard]] Numeric get_scattering_coeff() const;

  /// The scattering efficiency.
  [[nodiscard]] Numeric get_scattering_eff() const;

  /// The back-scattering coefficient.
  [[nodiscard]] Numeric get_backscattering_coeff() const;

  /// The back-scattering efficiency.
  [[nodiscard]] Numeric get_backscattering_eff() const;

  /// The extinction coefficient.
  [[nodiscard]] Numeric get_extinction_coeff() const;

  /// The extinction efficiency.
  [[nodiscard]] Numeric get_extinction_eff() const;

  /// The absorption coefficient.
  [[nodiscard]] Numeric get_absorption_coeff() const;

  /// The absorption efficiency.
  [[nodiscard]] Numeric get_absorption_eff() const;

  /// The phase function
  [[nodiscard]] Vector get_phase_function() const;

  /** Returns compact representation of scattering matrix
     *
     * The resulting matrix contains the elements F11, F12, F22, F33, F34, F44 of
     * the scattering matrix, as differential cross sections, for each angle for
     * which the scattering coefficients have been calculated.
     *
     * The elements are those of Eq. 4.77 from [1] in ARTS's Stokes convention
     * (concept.polarization): F is the Mueller matrix of the Jones matrix
     * diag(S_2, S_1) of the scattering-plane bases, so F34 = -Im(S_2 S_1*), the
     * negative of S_34 in [1], whose U has the other sign.
     *
     * @return The scattering matrix at every requested scattering angle
     */
  [[nodiscard]] CompactPlanarMuelmatVector get_scattering_matrix_compact() const;

 private:
  /** Calculation of the principal scattering properties.
     *
     * The implementation is essentially a port of the Fortran version
     * from [1].
     */
  void calculate_mie_parameters();

  Numeric       lambda_;  // The wavelength.
  Numeric       r_;       // The particle radius.
  Numeric       x_;       // The size parameter.
  Complex       n_;       // The relative refractive index.
  Numeric       q_sca_{0.0}, g_sca_{0.0}, q_ext_{0.0}, q_back_{0.0};
  Vector        theta_;
  ComplexVector s_1_, s_2_;
};
}  // namespace scattering

#endif  // SCATTERING_MIE_H_
