#include "mie.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>

#include <algorithm>
#include <cmath>

namespace scattering {

using Math::pow2;
using Math::pow3;

Complex refr_index_water_ellison07(Numeric frequency, Numeric temperature) {
  using Constant::pi;
  using Constant::two_pi;

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
  const Numeric t_cels = Conversion::kelvin2celsius(temperature);
  // static permittivity
  const Numeric epsilon_s = 87.9144 - 0.404399 * t_cels - 9.58726e-4 * pow2(t_cels) - 1.32802e-6 * pow3(t_cels);
  // Model parameters
  const Numeric delta1 = a1 * std::exp(-b1 * t_cels);
  const Numeric delta2 = a2 * std::exp(-b2 * t_cels);
  const Numeric delta3 = a3 * std::exp(-b3 * t_cels);
  const Numeric tau1   = c1 * std::exp(d1 / (t_cels + tc));
  const Numeric tau2   = c2 * std::exp(d2 / (t_cels + tc));
  const Numeric tau3   = c3 * std::exp(d3 / (t_cels + tc));
  const Numeric delta4 = p0 + p1 * t_cels + p2 * pow2(t_cels);
  const Numeric f0     = p3 + p4 * t_cels + p5 * pow2(t_cels) + p6 * pow3(t_cels);
  const Numeric tau4   = p7 + p8 * t_cels + p9 * pow2(t_cels) + p10 * pow3(t_cels);
  const Numeric delta5 = p11 + p12 * t_cels + p13 * pow2(t_cels);
  const Numeric f1     = p14 + p15 * t_cels + p16 * pow2(t_cels);
  const Numeric tau5   = p17 + p18 * t_cels + p19 * pow2(t_cels);

  const Numeric epsilon_real =
      epsilon_s -
      pow2((two_pi * frequency)) * (pow2(tau1) * delta1 / (1. + pow2(two_pi * frequency * tau1)) +
                                    pow2(tau2) * delta2 / (1. + pow2(two_pi * frequency * tau2)) +
                                    pow2(tau3) * delta3 / (1. + pow2(two_pi * frequency * tau3))) -
      pow2(two_pi * tau4) * delta4 / 2. *
          (frequency * (f0 + frequency) / (1. + pow2(two_pi * tau4 * (f0 + frequency))) -
           frequency * (f0 - frequency) / (1. + pow2(two_pi * tau4 * (f0 - frequency)))) -
      pow2(two_pi * tau5) * delta5 / 2. *
          (frequency * (f1 + frequency) / (1. + pow2(two_pi * tau5 * (f1 + frequency))) -
           frequency * (f1 - frequency) / (1. + pow2(two_pi * tau5 * (f1 - frequency))));
  // imaginary part of the complex permittivity of water (triple-debye + 2 resonances)
  const Numeric epsilon_imag =
      two_pi * frequency *
          (tau1 * delta1 / (1. + pow2(two_pi * frequency * tau1)) +
           tau2 * delta2 / (1. + pow2(two_pi * frequency * tau2)) +
           tau3 * delta3 / (1. + pow2(two_pi * frequency * tau3))) +
      pi * frequency * tau4 * delta4 *
          (1. / (1. + pow2(two_pi * tau4 * (f0 + frequency))) + 1. / (1. + pow2(two_pi * tau4 * (f0 - frequency)))) +
      pi * frequency * tau5 * delta5 *
          (1. / (1. + pow2(two_pi * tau5 * (f1 + frequency))) + 1. / (1. + pow2(two_pi * tau5 * (f1 - frequency))));

  return std::sqrt(Complex{epsilon_real, epsilon_imag});
}

Complex refr_index_ice_matzler06(Numeric frequency, Numeric temperature) {
  // some parametrization constants
  const Numeric B1 = 0.0207;
  const Numeric B2 = 1.16e-11;
  const Numeric b  = 335.;

  const Numeric deltabeta = std::exp(-9.963 + 0.0372 * (temperature - 273));
  const Numeric ebdt      = std::exp(b / temperature);
  const Numeric betam     = (B1 / temperature) * ebdt / ((ebdt - 1.) * (ebdt - 1.));

  const Numeric theta = 300. / temperature - 1;
  const Numeric alfa  = (0.00504 + 0.0062 * theta) * std::exp(-22.1 * theta);
  const Numeric reps  = 3.1884 + 9.1e-4 * (temperature - 273);

  Numeric f    = frequency / 1e9;
  Numeric beta = betam + B2 * f * f + deltabeta;
  Numeric ieps = alfa / f + beta * f;

  return std::sqrt(Complex{reps, ieps});
}

ComplexVector log_derivative(Complex rho, Size n_steps) {
  ComplexVector result(n_steps);

  result[n_steps - 1] = {0.0, 0.0};
  for (Size step = n_steps - 1; step > 0; --step) {
    Numeric n_f      = static_cast<Numeric>(step + 1);
    result[step - 1] = n_f / rho - Complex{1.0} / (result[step] + n_f / rho);
  }
  return result;
}

MieSphere::MieSphere(Numeric lambda, Numeric radius, Complex n, StridedConstVectorView theta)
    : lambda_(lambda),
      r_(radius),
      x_(Constant::two_pi * radius / lambda_),
      n_(n),
      theta_(theta),
      s_1_(theta.size(), Complex{0.0, 0.0}),
      s_2_(theta.size(), Complex{0.0, 0.0}) {
  calculate_mie_parameters();
}

MieSphere MieSphere::Liquid(Numeric frequency, Numeric temperature, Numeric radius, StridedConstVectorView theta) {
  return MieSphere(Constant::c / frequency, radius, refr_index_water_ellison07(frequency, temperature), theta);
}

MieSphere MieSphere::Ice(Numeric frequency, Numeric temperature, Numeric radius, StridedConstVectorView theta) {
  return MieSphere(Constant::c / frequency, radius, refr_index_ice_matzler06(frequency, temperature), theta);
}

Numeric MieSphere::get_scattering_coeff() const { return q_sca_ * Constant::pi * pow2(r_); }

Numeric MieSphere::get_scattering_eff() const { return q_sca_; }

Numeric MieSphere::get_backscattering_coeff() const { return q_back_ * Constant::pi * pow2(r_); }

Numeric MieSphere::get_backscattering_eff() const { return q_back_; }

Numeric MieSphere::get_extinction_coeff() const { return q_ext_ * Constant::pi * pow2(r_); }

Numeric MieSphere::get_extinction_eff() const { return q_ext_; }

Numeric MieSphere::get_absorption_coeff() const { return (q_ext_ - q_sca_) * Constant::pi * pow2(r_); }

Numeric MieSphere::get_absorption_eff() const { return q_ext_ - q_sca_; }

Vector MieSphere::get_phase_function() const {
  const Numeric k = Constant::two_pi / lambda_;
  Vector        results(s_1_.size());
  for (Size i = 0; i < results.size(); ++i) { results[i] = 0.5 * (std::norm(s_1_[i]) + std::norm(s_2_[i])) / (k * k); }
  return results;
}

CompactPlanarMuelmatVector MieSphere::get_scattering_matrix_compact() const {
  const Numeric              k = Constant::two_pi / lambda_;
  CompactPlanarMuelmatVector result(s_1_.size());
  for (Size i = 0; i < s_1_.size(); ++i) {
    const Numeric s11 = 0.5 * (std::norm(s_1_[i]) + std::norm(s_2_[i]));
    const Numeric s12 = 0.5 * (std::norm(s_2_[i]) - std::norm(s_1_[i]));
    const Numeric s33 = 0.5 * std::real(s_1_[i] * std::conj(s_2_[i]) + s_2_[i] * std::conj(s_1_[i]));
    const Numeric s34 = 0.5 * std::imag(s_1_[i] * std::conj(s_2_[i]) - s_2_[i] * std::conj(s_1_[i]));
    result[i]         = CompactPlanarMuelmat{s11, s12, s11, s33, s34, s33} / (k * k);
  }
  return result;
}

void MieSphere::calculate_mie_parameters() {
  const Complex rho         = x_ * n_;
  const Size    n_steps_x   = static_cast<Size>(x_ + 4.0 * std::cbrt(x_) + 2);
  const Size    n_steps_rho = static_cast<Size>(std::abs(rho));
  const Size    n_steps_D   = std::max(n_steps_x, n_steps_rho) + 15;

  const ComplexVector D = log_derivative(rho, n_steps_D);

  Numeric psi_2 = std::cos(x_);  // psi_{n - 2}
  Numeric psi_1 = std::sin(x_);  // psi_{n - 1}
  Numeric chi_2 = -std::sin(x_);
  Numeric chi_1 = std::cos(x_);
  Complex xi_1{psi_1, -chi_1};
  Complex a_n{0.0}, b_n{0.0}, a_n_1{0.0}, b_n_1{0.0}, q_back_acc{0.0};

  Vector mu = theta_;
  std::transform(mu.begin(), mu.end(), mu.begin(), [](Numeric x) { return std::cos(Conversion::deg2rad(x)); });
  const Size n_angs = theta_.size();

  Vector pi_1(n_angs, 0.0);
  Vector pi(n_angs, 1.0);
  Vector pi_2;

  // First step of iteration corresponds to n = 1
  for (Size step = 0; step < n_steps_x; ++step) {
    Numeric step_f = static_cast<Numeric>(step + 1);

    // Recurrence relation for Bessel functions.
    Numeric psi = (step_f * 2.0 - 1.0) * psi_1 / x_ - psi_2;
    Numeric chi = (step_f * 2.0 - 1.0) * chi_1 / x_ - chi_2;
    Complex xi{psi, -chi};

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
    q_sca_            += n2p1 * (std::norm(a_n) + std::norm(b_n));
    q_ext_            += n2p1 * (a_n.real() + b_n.real());
    q_back_acc        += n2p1 * (step % 2 == 0 ? -1.0 : 1.0) * (a_n - b_n);
    g_sca_            += (n2p1_nnp1 * (a_n.real() * b_n.real() + a_n.imag() * b_n.imag()));
    if (step > 0) {
      g_sca_ += (((step_f - 1.0) * (step_f + 1.0) / step_f) * (a_n_1.real() * a_n.real() + a_n_1.imag() * a_n.imag() +
                                                               b_n_1.real() * b_n.real() + b_n_1.imag() * b_n.imag()));
    }

    Vector tau(n_angs);
    for (Size i = 0; i < n_angs; ++i) {
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

    for (Size i = 0; i < n_angs; ++i) { pi[i] = (n2p1 * mu[i] * pi_1[i] - (step_f + 1.0) * pi_2[i]) / step_f; }
  }

  g_sca_  *= 2.0 / q_sca_;
  q_sca_  *= 2.0 / (x_ * x_);
  q_back_  = std::norm(q_back_acc) / (x_ * x_);
  q_ext_  *= 2.0 / (x_ * x_);
}
}  // namespace scattering
