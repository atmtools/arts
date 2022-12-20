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

#include "Eigen/Core"
#include <complex>
#include <iostream>
#include <numbers>
#include <vector>

#include "scattering/maths.h"
#include "arts_constexpr_math.h"

namespace scattering {

    using Math::pow2;
    using Math::pow3;

/** Calculates the refractive index of water following [3].
 *
 * The implementation is essentially copied from refraction.cc but included
 * here as a template function to simplify testing.
 *
 * @param frequency The frequency in GHz
 * @param temperature The temperature in K.
 */
template<typename Scalar>
std::complex<Scalar> refr_index_water_ellison07(
    Scalar frequency,
    Scalar temperature
  ) {

  Numeric pi = std::numbers::pi_v<double>;
  Numeric two_pi = pi * 2.0;

  // ELL07 model parameters - table 2 in Ellison (2007)
  constexpr Numeric a1 = 79.23882;
  constexpr Numeric a2 = 3.815866;
  constexpr Numeric a3 = 1.634967;
  constexpr Numeric tc = 133.1383;
  constexpr Numeric b1 = 0.004300598;
  constexpr Numeric b2 = 0.01117295;
  constexpr Numeric b3 = 0.006841548;
  constexpr Numeric c1 = 1.382264e-13;
  constexpr Numeric c2 = 3.510354e-16;
  constexpr Numeric c3 = 6.30035e-15;
  constexpr Numeric d1 = 652.7648;
  constexpr Numeric d2 = 1249.533;
  constexpr Numeric d3 = 405.5169;
  constexpr Numeric p0 = 0.8379692;
  constexpr Numeric p1 = -0.006118594;
  constexpr Numeric p2 = -0.000012936798;
  constexpr Numeric p3 = 4235901000000.0;
  constexpr Numeric p4 = -14260880000.0;
  constexpr Numeric p5 = 273815700.0;
  constexpr Numeric p6 = -1246943.0;
  constexpr Numeric p7 = 9.618642e-14;
  constexpr Numeric p8 = 1.795786e-16;
  constexpr Numeric p9 = -9.310017E-18;
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
  const Numeric tau1 = c1 * exp(d1 / (t_cels + tc));
  const Numeric tau2 = c2 * exp(d2 / (t_cels + tc));
  const Numeric tau3 = c3 * exp(d3 / (t_cels + tc));
  const Numeric delta4 = p0 + p1 * t_cels + p2 * pow2(t_cels);
  const Numeric f0 = p3 + p4 * t_cels + p5 * pow2(t_cels) + p6 * pow3(t_cels);
  const Numeric tau4 =
      p7 + p8 * t_cels + p9 * pow2(t_cels) + p10 * pow3(t_cels);
  const Numeric delta5 = p11 + p12 * t_cels + p13 * pow2(t_cels);
  const Numeric f1 = p14 + p15 * t_cels + p16 * pow2(t_cels);
  const Numeric tau5 = p17 + p18 * t_cels + p19 * pow2(t_cels);

  const Numeric epsilon_real =
        epsilon_s -
        pow2((two_pi * frequency)) * (pow2(tau1) * delta1 / (1. + pow2(two_pi * frequency * tau1)) +
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
std::complex<Scalar> refr_index_ice_matzler06(
    Scalar frequency,
    Scalar temperature
    ) {
  // some parametrization constants
  const Scalar B1 = 0.0207;
  const Scalar B2 = 1.16e-11;
  const Scalar b = 335.;

  const Scalar deltabeta = exp(-9.963 + 0.0372 * (temperature - 273));
  const Scalar ebdt = exp(b / temperature);
  const Scalar betam = (B1 / temperature) * ebdt / ((ebdt - 1.) * (ebdt - 1.));

  const Scalar theta = 300. / temperature - 1;
  const Scalar alfa = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);
  const Scalar reps = 3.1884 + 9.1e-4 * (temperature - 273);

  Scalar f = frequency / 1e9;
  Scalar beta = betam + B2 * f * f + deltabeta;
  Scalar ieps = alfa / f + beta * f;

  std::complex<Scalar> n = sqrt(std::complex<Scalar>{reps, ieps});
  return n;
}

template <typename Scalar>
using Array = Eigen::ArrayX<Scalar>;

/** Calculation of the log-derivative using downward recursion.
 *
 * The logarithmic derivative is defined by Eq. (4.89) in [1]. The suggested
 * number of steps in the recurrence is:
 *
 *  n_steps = max(x + 4 * x^(1/3) + 2, abs(x * n)) + 15
 *
 * @param rho: The product of relative refractive index n and size parameter x.
 * @param n_steps: The number of steps in the downward recursion.
 * @return An Eigen vector containing the 'n_steps' calculated complex values
 * of the logarithmic derivative
 */
template <std::floating_point Scalar>
Array<std::complex<Scalar>> log_derivative(std::complex<Scalar> rho, size_t n_steps) {
    Array<std::complex<Scalar>> result{n_steps};

    result[n_steps - 1] = {0.0, 0.0};
    for (size_t step = n_steps - 1; step > 0; --step) {
        Scalar n_f = static_cast<Scalar>(step + 1);
        result[step - 1] = n_f / rho - static_cast<std::complex<Scalar>>(1.0) / (result[step] + n_f / rho);
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
     * @param theta Array of scattering angles for which to compute the scattering matrix.
     */
  MieSphere(
      Scalar lambda,
      Scalar radius,
      std::complex<Scalar> n,
      Array<Scalar> theta)
      : lambda_(lambda), r_(radius), x_(2.0 * std::numbers::pi_v<Scalar> * radius / lambda_), n_(n), theta_(theta), s_1_{Array<std::complex<Scalar>>::Constant(theta.size(), {0.0, 0.0})},
        s_2_{Array<std::complex<Scalar>>::Constant(theta.size(), {0.0, 0.0})} {
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
    static MieSphere Liquid(
        Scalar frequency,
        Scalar temperature,
        Scalar radius,
        Array<Scalar> theta
        ) {
        Scalar c = 2.99792458e8;
        Scalar lambda = c / frequency;
        std::complex n = refr_index_water_ellison07(frequency, temperature);
        return MieSphere(
            lambda,
            radius,
            n,
            theta);
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
    static MieSphere Ice(
        Scalar frequency,
        Scalar  temperature,
        Scalar radius,
        Array<Scalar> theta
        ) {
        Scalar c = 2.99792458e8;
        Scalar lambda = c / frequency;
        std::complex n = refr_index_ice_matzler06(frequency, temperature);
        return MieSphere(
            lambda,
            radius,
            n,
            theta);
    }

    /// The scattering coefficient.
    Scalar get_scattering_coeff() {
        return q_sca_ * std::numbers::pi_v<Scalar> * pow(r_, 2);
    }

    /// The scattering efficiency.
    Scalar get_scattering_eff() {
        return q_sca_;
    }

    /// The back-scattering coefficient.
    Scalar get_backscattering_coeff() {
        return q_back_ * std::numbers::pi_v<Scalar> * pow(r_, 2);
    }

    /// The back-scattering efficiency.
    Scalar get_backscattering_eff() {
        return q_back_;
    }

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
    Array<Scalar> get_phase_function() {
        Scalar k = std::numbers::pi_v<Scalar> * 2.0 / lambda_;
        return 0.5 * (s_1_.abs().pow(2) + s_2_.abs().pow(2)) / (k * k);
    }

    /** Returns compact representation of scattering matrix
     *
     * The resulting matrix contains the four elements S_11, S_12, S_33, S_34 of
     * the scattering matrix for each angle for which the scattering coefficients have
     * been calculated.
     *
     * This method implements Eq. 4.77 from [1].
     *
     * @return An Eigen::Matrix holding the four scattering matrix elements along its columns
     * for all requested scattering angles.
     */
    math::Matrix<Scalar> get_scattering_matrix_compact() {
        Scalar k = std::numbers::pi_v<Scalar> * 2.0 / lambda_;
        Array<Scalar> s_1_abs = s_1_.abs();
        Array<Scalar> s_2_abs = s_2_.abs();
        Array<Scalar> s11 = 0.5 * (s_1_abs.pow(2) + s_2_abs.pow(2));
        Array<Scalar> s12 = 0.5 * (s_2_abs.pow(2) - s_1_abs.pow(2));
        Array<Scalar> s33 = 0.5 * (s_1_ * s_2_.conjugate() + s_2_ * s_1_.conjugate()).real();
        Array<Scalar> s34 =  -0.5 * (s_1_ * s_2_.conjugate() - s_2_ * s_1_.conjugate()).imag();

        math::Matrix<Scalar> result{theta_.size(), 6};
        result(Eigen::all, 0) = s11;
        result(Eigen::all, 1) = s12;
        result(Eigen::all, 2) = s11;
        result(Eigen::all, 3) = s33;
        result(Eigen::all, 4) = s34;
        result(Eigen::all, 5) = s33;
        return result / (k * k);
    }

  private:

    /** Calculation of the principal scattering properties.
     *
     * The implementation is essentially a port of the Fortran version
     * from [1].
     */
    void calculate_mie_parameters() {

        std::complex<Scalar> rho = x_ * n_;
        size_t n_steps_x = static_cast<size_t>(x_ + 4.0 * pow(x_, 1.0 / 3.0) + 2);
        size_t n_steps_rho = static_cast<size_t>(abs(rho));
        size_t n_steps_D = std::max(n_steps_x, n_steps_rho) + 15;


        auto D = log_derivative<Scalar>(rho, n_steps_D);

        Scalar psi_2 = cos(x_); // psi_{n - 2}
        Scalar psi_1 = sin(x_); // psi_{n - 1}
        Scalar chi_2 = -sin(x_);
        Scalar chi_1 = cos(x_);
        std::complex xi_1{psi_1, -chi_1};
        std::complex<Scalar> a_n{0.0}, b_n{0.0}, a_n_1{0.0}, b_n_1{0.0}, q_back_acc{0.0};

        Array<Scalar> mu = theta_.array().cos();
        Array<Scalar> pi_1 = Array<Scalar>::Zero(mu.size());
        Array<Scalar> pi = pi_1 + 1.0;
        Array<Scalar> pi_2;

        // First step of iteration corresponds to n = 1
        for (size_t step = 0; step < n_steps_x; ++step) {

            Scalar step_f = static_cast<Scalar>(step + 1);

            // Recurrence relation for Bessel functions.
            Scalar psi = (step_f * 2.0 - 1.0) * psi_1 / x_ - psi_2;
            Scalar chi = (step_f * 2.0 - 1.0) * chi_1 / x_ - chi_2;
            std::complex xi{psi, -chi};

            a_n_1 = a_n;
            b_n_1 = b_n;

            // Mie parameters
            a_n = (D[step] / n_ + step_f / x_) * psi - psi_1;
            a_n /= ((D[step] / n_ + step_f / x_) * xi - xi_1);
            b_n = (D[step] * n_ + step_f / x_) * psi - psi_1;
            b_n /= ((D[step] * n_ + step_f / x_) * xi - xi_1);

            // Scattering parameters
            Scalar n2p1 = 2.0 * step_f + 1.0;
            Scalar n2p1_nnp1 = n2p1 / (step_f * (step_f + 1.0));
            q_sca_ +=
                n2p1 * (pow(abs(a_n), 2) + pow(abs(b_n), 2));
            q_ext_ += n2p1 * (a_n.real() + b_n.real());
            q_back_acc += n2p1 * pow(-1, step + 1) * (a_n - b_n);
            g_sca_ += (
                n2p1_nnp1 *
                (a_n.real() * b_n.real() + a_n.imag() * b_n.imag())
                );
            if (step > 0) {
              g_sca_ += (((step_f - 1.0) * (step_f + 1.0) / step_f) *
                        (a_n_1.real() * a_n.real() + a_n_1.imag() * a_n.imag() +
                         b_n_1.real() * b_n.real() + b_n_1.imag() * b_n.imag()));
            }


            Array<Scalar> tau = step_f * mu * pi - (step_f + 1.0) * pi_1;
            s_1_ += n2p1_nnp1 * (a_n * pi + b_n * tau);
            s_2_ += n2p1_nnp1 * (a_n * tau + b_n * pi);

            // Loop update.
            psi_2 = psi_1;
            psi_1 = psi;
            chi_2 = chi_1;
            chi_1 = chi;
            xi_1 = xi;

            pi_2 = pi_1;
            pi_1 = pi;
            pi = (n2p1 * mu * pi_1 - (step_f + 1.0) * pi_2) / step_f;
        }

        g_sca_ *= 2.0 / q_sca_;
        q_sca_ *= 2.0 / (x_ * x_);
        q_back_ = pow(abs(q_back_acc), 2) / (x_ * x_);
        q_ext_ *= 2.0 / (x_ * x_);
    }

    Scalar lambda_; // The wavelength.
    Scalar r_; // The particle radius.
    Scalar x_; // The size parameter.
    std::complex<Scalar> n_; // The relative refractive index.
    Scalar q_sca_{0.0}, g_sca_{0.0}, q_ext_{0.0}, q_back_{0.0};
    Array<Scalar> theta_;
    Array<std::complex<Scalar>> s_1_, s_2_;
};

}
