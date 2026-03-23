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

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <matpack.h>

#include <complex>
#include <numbers>

namespace scattering {

using Math::pow2;
using Math::pow3;

/** Calculates the complex refractive index of liquid water following Liebe 1993.
 *
 * This is the widely-used double-Debye model from C. Mätzler's epswater93,
 * following the paper version in AGARD CP-May93.  Used in radar applications
 * for computing the dielectric factor |K|² at a reference temperature.
 *
 * @param frequency The frequency in Hz
 * @param temperature The temperature in K
 * @return The complex refractive index sqrt(epsilon)
 */
template <typename Scalar>
std::complex<Scalar> refr_index_water_liebe93(Scalar frequency,
                                              Scalar temperature) {
  const Scalar theta = 1.0 - 300.0 / temperature;
  const Scalar e0    = 77.66 - 103.3 * theta;
  const Scalar e1    = 0.0671 * e0;
  const Scalar f1    = 20.2 + 146.0 * theta + 316.0 * theta * theta;
  const Scalar e2    = 3.52;
  const Scalar f2    = 39.8 * f1;

  const std::complex<Scalar> ifGHz(0.0, frequency / 1e9);

  return std::sqrt(e2 + (e1 - e2) / (1.0 - ifGHz / f2) +
                   (e0 - e1) / (1.0 - ifGHz / f1));
}

/** Calculates the complex dielectric factor K = (n²-1)/(n²+2).
 *
 * @param n Complex refractive index
 * @return The complex K factor
 */
template <typename Scalar>
std::complex<Scalar> dielectric_factor(std::complex<Scalar> n) {
  auto n2 = n * n;
  return (n2 - 1.0) / (n2 + 2.0);
}

/** Calculates the refractive index of water following [3].
 *
 * The implementation is essentially copied from refraction.cc but included
 * here as a template function to simplify testing.
 *
 * @param frequency The frequency in GHz
 * @param temperature The temperature in K.
 */
template <typename Scalar>
std::complex<Scalar> refr_index_water_ellison07(Scalar frequency,
                                                Scalar temperature) {
  Numeric pi     = std::numbers::pi_v<double>;
  Numeric two_pi = pi * 2.0;

  // ELL07 model parameters - table 2 in Ellison (2007)
  constexpr Numeric a1  = 79.23882;
  constexpr Numeric a2  = 3.815866;
  constexpr Numeric a3  = 1.634967;
  constexpr Numeric tc  = 133.1383;
  constexpr Numeric b1  = 0.004300598;
  constexpr Numeric b2  = 0.01117295;
  constexpr Numeric b3  = 0.006841548;
  constexpr Numeric c1  = 1.382264e-13;
  constexpr Numeric c2  = 3.510354e-16;
  constexpr Numeric c3  = 6.30035e-15;
  constexpr Numeric d1  = 652.7648;
  constexpr Numeric d2  = 1249.533;
  constexpr Numeric d3  = 405.5169;
  constexpr Numeric p0  = 0.8379692;
  constexpr Numeric p1  = -0.006118594;
  constexpr Numeric p2  = -0.000012936798;
  constexpr Numeric p3  = 4235901000000.0;
  constexpr Numeric p4  = -14260880000.0;
  constexpr Numeric p5  = 273815700.0;
  constexpr Numeric p6  = -1246943.0;
  constexpr Numeric p7  = 9.618642e-14;
  constexpr Numeric p8  = 1.795786e-16;
  constexpr Numeric p9  = -9.310017E-18;
  constexpr Numeric p10 = 1.655473e-19;
  constexpr Numeric p11 = 0.6165532;
  constexpr Numeric p12 = 0.007238532;
  constexpr Numeric p13 = -0.00009523366;
  constexpr Numeric p14 = 15983170000000.0;
  constexpr Numeric p15 = -74413570000.0;
  constexpr Numeric p16 = 497448000.0;
  constexpr Numeric p17 = 2.882476e-14;
  constexpr Numeric p18 = -3.142118e-16;
  constexpr Numeric p19 = 3.528051e-18;

  // Temperature in celsius
  const Numeric t_cels = temperature - 273.15;
  // static permittivity
  const Numeric epsilon_s = 87.9144 - 0.404399 * t_cels -
                            9.58726e-4 * pow2(t_cels) -
                            1.32802e-6 * pow3(t_cels);
  // Model parameters
  const Numeric delta1 = a1 * exp(-b1 * t_cels);
  const Numeric delta2 = a2 * exp(-b2 * t_cels);
  const Numeric delta3 = a3 * exp(-b3 * t_cels);
  const Numeric tau1   = c1 * exp(d1 / (t_cels + tc));
  const Numeric tau2   = c2 * exp(d2 / (t_cels + tc));
  const Numeric tau3   = c3 * exp(d3 / (t_cels + tc));
  const Numeric delta4 = p0 + p1 * t_cels + p2 * pow2(t_cels);
  const Numeric f0 = p3 + p4 * t_cels + p5 * pow2(t_cels) + p6 * pow3(t_cels);
  const Numeric tau4 =
      p7 + p8 * t_cels + p9 * pow2(t_cels) + p10 * pow3(t_cels);
  const Numeric delta5 = p11 + p12 * t_cels + p13 * pow2(t_cels);
  const Numeric f1     = p14 + p15 * t_cels + p16 * pow2(t_cels);
  const Numeric tau5   = p17 + p18 * t_cels + p19 * pow2(t_cels);

  const Numeric epsilon_real =
      epsilon_s -
      pow2((two_pi * frequency)) *
          (pow2(tau1) * delta1 / (1. + pow2(two_pi * frequency * tau1)) +
           pow2(tau2) * delta2 / (1. + pow2(two_pi * frequency * tau2)) +
           pow2(tau3) * delta3 / (1. + pow2(two_pi * frequency * tau3))) -
      pow2(two_pi * tau4) * delta4 / 2. *
          (frequency * (f0 + frequency) /
               (1. + pow2(two_pi * tau4 * (f0 + frequency))) -
           frequency * (f0 - frequency) /
               (1. + pow2(two_pi * tau4 * (f0 - frequency)))) -
      pow2(two_pi * tau5) * delta5 / 2. *
          (frequency * (f1 + frequency) /
               (1. + pow2(two_pi * tau5 * (f1 + frequency))) -
           frequency * (f1 - frequency) /
               (1. + pow2(two_pi * tau5 * (f1 - frequency))));
  // imaginary part of the complex permittivity of water (triple-debye + 2 resonances)
  const Numeric epsilon_imag =
      two_pi * frequency *
          (tau1 * delta1 / (1. + pow2(two_pi * frequency * tau1)) +
           tau2 * delta2 / (1. + pow2(two_pi * frequency * tau2)) +
           tau3 * delta3 / (1. + pow2(two_pi * frequency * tau3))) +
      pi * frequency * tau4 * delta4 *
          (1. / (1. + pow2(two_pi * tau4 * (f0 + frequency))) +
           1. / (1. + pow2(two_pi * tau4 * (f0 - frequency)))) +
      pi * frequency * tau5 * delta5 *
          (1. / (1. + pow2(two_pi * tau5 * (f1 + frequency))) +
           1. / (1. + pow2(two_pi * tau5 * (f1 - frequency))));

  return sqrt(std::complex<Scalar>{epsilon_real, epsilon_imag});
}

/** Calculates the refractive index of water following [4].
 *
 * The implementation is copied from refraction.cc but included
 * here as a template function to simplify testing.
 *
 * @param frequency The frequency in Hz
 * @param temperature The temperature in K.
 */
template <std::floating_point Scalar>
std::complex<Scalar> refr_index_ice_matzler06(Scalar frequency,
                                              Scalar temperature) {
  // some parametrization constants
  const Scalar B1 = 0.0207;
  const Scalar B2 = 1.16e-11;
  const Scalar b  = 335.;

  const Scalar deltabeta = exp(-9.963 + 0.0372 * (temperature - 273));
  const Scalar ebdt      = exp(b / temperature);
  const Scalar betam = (B1 / temperature) * ebdt / ((ebdt - 1.) * (ebdt - 1.));

  const Scalar theta = 300. / temperature - 1;
  const Scalar alfa  = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);
  const Scalar reps  = 3.1884 + 9.1e-4 * (temperature - 273);

  Scalar f    = frequency / 1e9;
  Scalar beta = betam + B2 * f * f + deltabeta;
  Scalar ieps = alfa / f + beta * f;

  std::complex<Scalar> n = sqrt(std::complex<Scalar>{reps, ieps});
  return n;
}

template <typename Scalar>
using VectorGen = matpack::data_t<Scalar, 1>;
template <typename Scalar>
using MatrixGen = matpack::data_t<Scalar, 2>;

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
template <std::floating_point Scalar>
VectorGen<std::complex<Scalar>> log_derivative(std::complex<Scalar> rho,
                                               size_t n_steps) {
  VectorGen<std::complex<Scalar>> result(n_steps);

  result[n_steps - 1] = {0.0, 0.0};
  for (size_t step = n_steps - 1; step > 0; --step) {
    Scalar n_f       = static_cast<Scalar>(step + 1);
    result[step - 1] = n_f / rho - static_cast<std::complex<Scalar>>(1.0) /
                                       (result[step] + n_f / rho);
  }
  return result;
}

/** A Mie sphere
 *
 * The MieSphere class calculates the Mie parameters of a given sphere and
 * provides access to to its scattering properties.
 */
template <std::floating_point Scalar>
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
  MieSphere(Scalar lambda,
            Scalar radius,
            std::complex<Scalar> n,
            StridedConstVectorView theta)
      : lambda_(lambda),
        r_(radius),
        x_(2.0 * std::numbers::pi_v<Scalar> * radius / lambda_),
        n_(n),
        theta_(theta),
        s_1_(theta.size()),
        s_2_(theta.size()) {
    s_1_ = std::complex<Scalar>(0.0, 0.0);
    s_2_ = std::complex<Scalar>(0.0, 0.0);
    calculate_mie_parameters();
  }

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
  static MieSphere Liquid(Scalar frequency,
                          Scalar temperature,
                          Scalar radius,
                          StridedConstVectorView theta) {
    Scalar c       = 2.99792458e8;
    Scalar lambda  = c / frequency;
    std::complex n = refr_index_water_ellison07(frequency, temperature);
    return MieSphere(lambda, radius, n, theta);
  }

  /** Calculate scattering properties of a spherical ice particle.
     *
     * @param temperature The ambient temperature in K.
     * @param frequency The frequency in Hz.
     * @param radius The radius of the sphere.
     * @param theta The scattering angles for which to calculate the
     *     the scattering properties.
     * @return A MieSphere object representing the scattering properties
     * of the corresponding spherical ice droplet.
     */
  static MieSphere Ice(Scalar frequency,
                       Scalar temperature,
                       Scalar radius,
                       VectorGen<Scalar> theta) {
    Scalar c       = 2.99792458e8;
    Scalar lambda  = c / frequency;
    std::complex n = refr_index_ice_matzler06(frequency, temperature);
    return MieSphere(lambda, radius, n, theta);
  }

  /// The scattering coefficient.
  Scalar get_scattering_coeff() {
    return q_sca_ * std::numbers::pi_v<Scalar> * pow(r_, 2);
  }

  /// The scattering efficiency.
  Scalar get_scattering_eff() { return q_sca_; }

  /// The back-scattering coefficient.
  Scalar get_backscattering_coeff() {
    return q_back_ * std::numbers::pi_v<Scalar> * pow(r_, 2);
  }

  /// The back-scattering efficiency.
  Scalar get_backscattering_eff() { return q_back_; }

  /// The extinction coefficient.
  Scalar get_extinction_coeff() {
    return q_ext_ * std::numbers::pi_v<Scalar> * pow(r_, 2);
  }

  /// The extinction efficiency.
  Scalar get_extinction_eff() { return q_ext_; }

  /// The extinction coefficient.
  Scalar get_absorption_coeff() {
    return (q_ext_ - q_sca_) * std::numbers::pi_v<Scalar> * pow(r_, 2);
  }

  /// The extinction efficiency.
  Scalar get_absorption_eff() { return q_ext_ - q_sca_; }

  /// The phase function
  VectorGen<Scalar> get_phase_function() {
    Scalar k = std::numbers::pi_v<Scalar> * 2.0 / lambda_;
    VectorGen<Scalar> results(s_1_.size());
    for (Index i = 0; i < results.size(); ++i) {
      results[i] =
          0.5 *
          (std::pow(std::abs(s_1_[i]), 2) + std::pow(std::abs(s_2_[i]), 2)) /
          (k * k);
    }
    return results;
  }

  /** Returns compact representation of scattering matrix
     *
     * The resulting matrix contains the four elements S_11, S_12, S_33, S_34 of
     * the scattering matrix for each angle for which the scattering coefficients have
     * been calculated.
     *
     * This method implements Eq. 4.77 from [1].
     *
     * @return A matrix holding the four scattering matrix elements along its columns
     * for all requested scattering angles.
     */
  MatrixGen<Scalar> get_scattering_matrix_compact() {
    Scalar k       = std::numbers::pi_v<Scalar> * 2.0 / lambda_;
    Index n_angles = s_1_.size();
    VectorGen<Scalar> s11(n_angles);
    VectorGen<Scalar> s12(n_angles);
    VectorGen<Scalar> s33(n_angles);
    VectorGen<Scalar> s34(n_angles);

    for (Index i = 0; i < n_angles; ++i) {
      s11[i] = 0.5 * (std::pow(std::abs(s_1_[i]), 2) +
                      std::pow(std::abs(s_2_[i]), 2));
      s12[i] = 0.5 * (std::pow(std::abs(s_2_[i]), 2) -
                      std::pow(std::abs(s_1_[i]), 2));
      s33[i] = 0.5 * std::real(s_1_[i] * std::conj(s_2_[i]) +
                               s_2_[i] * std::conj(s_1_[i]));
      s34[i] = -0.5 * std::imag(s_1_[i] * std::conj(s_2_[i]) -
                                s_2_[i] * std::conj(s_1_[i]));
    }

    MatrixGen<Scalar> result{theta_.size(), 6};
    result[joker, 0]  = s11;
    result[joker, 1]  = s12;
    result[joker, 2]  = s11;
    result[joker, 3]  = s33;
    result[joker, 4]  = s34;
    result[joker, 5]  = s33;
    result           /= (k * k);
    return result;
  }

 private:
  /** Calculation of the principal scattering properties.
     *
     * The implementation is essentially a port of the Fortran version
     * from [1].
     */
  void calculate_mie_parameters() {
    std::complex<Numeric> rho = x_ * n_;
    size_t n_steps_x   = static_cast<size_t>(x_ + 4.0 * pow(x_, 1.0 / 3.0) + 2);
    size_t n_steps_rho = static_cast<size_t>(abs(rho));
    size_t n_steps_D   = std::max(n_steps_x, n_steps_rho) + 15;

    auto D = log_derivative<Numeric>(rho, n_steps_D);

    Numeric psi_2 = cos(x_);  // psi_{n - 2}
    Numeric psi_1 = sin(x_);  // psi_{n - 1}
    Numeric chi_2 = -sin(x_);
    Numeric chi_1 = cos(x_);
    std::complex xi_1{psi_1, -chi_1};
    std::complex<Numeric> a_n{0.0}, b_n{0.0}, a_n_1{0.0}, b_n_1{0.0},
        q_back_acc{0.0};

    VectorGen<Numeric> mu = theta_;
    std::transform(mu.begin(), mu.end(), mu.begin(), [](Numeric x) {
      return cos(Conversion::deg2rad(x));
    });
    Index n_angs = theta_.size();

    VectorGen<Numeric> pi_1(n_angs);
    pi_1                   = 0.0;
    VectorGen<Numeric> pi  = pi_1;
    pi                    += 1.0;
    VectorGen<Numeric> pi_2;

    // First step of iteration corresponds to n = 1
    for (size_t step = 0; step < n_steps_x; ++step) {
      Numeric step_f = static_cast<Numeric>(step + 1);

      // Recurrence relation for Bessel functions.
      Numeric psi = (step_f * 2.0 - 1.0) * psi_1 / x_ - psi_2;
      Numeric chi = (step_f * 2.0 - 1.0) * chi_1 / x_ - chi_2;
      std::complex xi{psi, -chi};

      a_n_1 = a_n;
      b_n_1 = b_n;

      // Mie parameters
      a_n  = (D[step] / n_ + step_f / x_) * psi - psi_1;
      a_n /= ((D[step] / n_ + step_f / x_) * xi - xi_1);
      b_n  = (D[step] * n_ + step_f / x_) * psi - psi_1;
      b_n /= ((D[step] * n_ + step_f / x_) * xi - xi_1);

      // Scattering parameters
      Numeric n2p1       = 2.0 * step_f + 1.0;
      Numeric n2p1_nnp1  = n2p1 / (step_f * (step_f + 1.0));
      q_sca_            += n2p1 * (pow(abs(a_n), 2) + pow(abs(b_n), 2));
      q_ext_            += n2p1 * (a_n.real() + b_n.real());
      q_back_acc        += n2p1 * std::pow(-1, step + 1) * (a_n - b_n);
      g_sca_ +=
          (n2p1_nnp1 * (a_n.real() * b_n.real() + a_n.imag() * b_n.imag()));
      if (step > 0) {
        g_sca_ += (((step_f - 1.0) * (step_f + 1.0) / step_f) *
                   (a_n_1.real() * a_n.real() + a_n_1.imag() * a_n.imag() +
                    b_n_1.real() * b_n.real() + b_n_1.imag() * b_n.imag()));
      }

      VectorGen<Numeric> tau(n_angs);
      for (Index i = 0; i < n_angs; ++i) {
        tau[i]   = step_f * mu[i] * pi[i] - (step_f + 1.0) * pi_1[i];
        s_1_[i] += n2p1_nnp1 * (a_n * pi[i] + b_n * tau[i]);
        s_2_[i] += n2p1_nnp1 * (a_n * tau[i] + b_n * pi[i]);
      }

      // Loop update.
      psi_2 = psi_1;
      psi_1 = psi;
      chi_2 = chi_1;
      chi_1 = chi;
      xi_1  = xi;

      pi_2 = pi_1;
      pi_1 = pi;

      for (Index i = 0; i < n_angs; ++i) {
        pi[i] = (n2p1 * mu[i] * pi_1[i] - (step_f + 1.0) * pi_2[i]) / step_f;
      }
    }

    g_sca_  *= 2.0 / q_sca_;
    q_sca_  *= 2.0 / (x_ * x_);
    q_back_  = pow(abs(q_back_acc), 2) / (x_ * x_);
    q_ext_  *= 2.0 / (x_ * x_);
  }

  Scalar lambda_;           // The wavelength.
  Scalar r_;                // The particle radius.
  Scalar x_;                // The size parameter.
  std::complex<Scalar> n_;  // The relative refractive index.
  Scalar q_sca_{0.0}, g_sca_{0.0}, q_ext_{0.0}, q_back_{0.0};
  VectorGen<Scalar> theta_;
  VectorGen<std::complex<Scalar>> s_1_, s_2_;
};

/** Rayleigh sphere — analytic scattering in the small-particle limit.
 *
 * Valid when the size parameter x = π D / λ ≪ 1.
 * Provides the same interface as MieSphere for the quantities needed by
 * the radar forward model (extinction, backscattering, phase matrix).
 *
 * All formulae follow Bohren & Huffman (1983), Eqs. 5.7–5.16.
 */
template <std::floating_point Scalar>
class RayleighSphere {
 public:
  /** Create a Rayleigh sphere.
   *
   * @param frequency Frequency in Hz
   * @param temperature Temperature in K
   * @param diameter Particle diameter in m
   * @param n Complex refractive index of the particle material
   */
  RayleighSphere(Scalar frequency,
                 [[maybe_unused]] Scalar temperature,
                 Scalar diameter,
                 std::complex<Scalar> n)
      : frequency_(frequency), diameter_(diameter), K_(dielectric_factor(n)) {
    constexpr Scalar pi = std::numbers::pi_v<Scalar>;
    const Scalar la     = Constant::speed_of_light / frequency;
    const Scalar r      = diameter / 2;
    const Scalar x      = pi * diameter / la;

    // Absorption cross-section (Bohren & Huffman Eq. 5.11)
    c_abs_ = pi * r * r * 4.0 * x * K_.imag();

    // Scattering cross-section (Bohren & Huffman Eq. 5.8)
    c_sca_ = pi * r * r * (8.0 / 3.0) * x * x * x * x * std::norm(K_);

    // Backscatter cross-section per steradian:
    //   dσ_back/dΩ = π⁴ D⁶ |K|² / (4 λ⁴)
    // This is the (1,1) element of the 4×4 phase matrix at 180°.
    z11_ = Math::pow4(pi) * std::pow(diameter, Scalar{6}) * std::norm(K_) /
           (4.0 * Math::pow4(la));
  }

  /** Create a Rayleigh sphere for liquid water (Liebe 1993 model).
   *
   * @param frequency Frequency in Hz
   * @param temperature Temperature in K
   * @param diameter Particle diameter in m
   */
  static RayleighSphere Liquid(Scalar frequency,
                               Scalar temperature,
                               Scalar diameter) {
    auto n = refr_index_water_liebe93(frequency, temperature);
    return RayleighSphere(frequency, temperature, diameter, n);
  }

  /// Extinction cross-section [m²]
  Scalar get_extinction_coeff() const { return c_abs_ + c_sca_; }

  /// Absorption cross-section [m²]
  Scalar get_absorption_coeff() const { return c_abs_; }

  /// Scattering cross-section [m²]
  Scalar get_scattering_coeff() const { return c_sca_; }

  /// Backscatter cross-section per steradian (Z11 at 180°) [m² sr⁻¹]
  Scalar get_backscatter_z11() const { return z11_; }

  /// Dielectric factor K = (n²-1)/(n²+2)
  std::complex<Scalar> get_K() const { return K_; }

  /// |K|²
  Scalar get_K2() const { return std::norm(K_); }

  /** 4×4 backscatter phase matrix at exactly 180°.
   *
   * For a Rayleigh sphere, the backscatter matrix is:
   *   diag(z11, z11, -z11, -z11)
   * (sign flip in S33, S44 due to exact backscatter geometry).
   */
  std::array<Scalar, 16> get_backscatter_matrix() const {
    // Row-major 4×4
    return {z11_, 0, 0, 0, 0, z11_, 0, 0, 0, 0, -z11_, 0, 0, 0, 0, -z11_};
  }

 private:
  Scalar frequency_;
  Scalar diameter_;
  std::complex<Scalar> K_;
  Scalar c_abs_{0};
  Scalar c_sca_{0};
  Scalar z11_{0};
};

}  // namespace scattering

#endif  // SCATTERING_MIE_H_
