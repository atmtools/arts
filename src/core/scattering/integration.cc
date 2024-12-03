#include <integration.h>

#ifndef ARTS_NO_SHTNS
#include <fftw3.h>
#endif

namespace scattering {
void GaussLegendreQuadrature::calculate_nodes_and_weights() {
  const Index n            = degree_;
  const Index n_half_nodes = (n + 1) / 2;
  const Index n_max_iter   = 10;
  Numeric x, x_old, p_l, p_l_1, p_l_2, dp_dx, n_f;
  n_f = static_cast<Numeric>(n);

  for (int i = 1; i <= n_half_nodes; ++i) {
    p_l   = pi_v<Numeric>;
    p_l_1 = 2.0 * n_f;
    //
    // Initial guess.
    //
    x = -(1.0 - (n_f - 1.0) / (p_l_1 * p_l_1 * p_l_1)) *
        cos((p_l * (4.0 * i - 1.0)) / (4.0 * n_f + 2.0));

    //
    // Evaluate Legendre Polynomial and its derivative at node.
    //
    for (Index j = 0; j < n_max_iter; ++j) {
      p_l   = x;
      p_l_1 = 1.0;
      for (int l = 2; l <= n; ++l) {
        // Legendre recurrence relation
        p_l_2 = p_l_1;
        p_l_1 = p_l;
        p_l   = ((2.0 * l - 1.0) * x * p_l_1 - (l - 1.0) * p_l_2) / l;
      }
      dp_dx = ((1.0 - x) * (1.0 + x)) / (n_f * (p_l_1 - x * p_l));
      x_old = x;

      //
      // Perform Newton step.
      //
      x       -= p_l * dp_dx;
      auto dx  = x - x_old;
      if (detail::small(std::abs(dx * (x + x_old)), 1e-10)) {
        break;
      }
    }
    nodes_[i - 1]   = x;
    weights_[i - 1] = 2.0 * dp_dx * dp_dx / ((1.0 - x) * (1.0 + x));
    nodes_[n - i]   = -x;
    weights_[n - i] = weights_[i - 1];
  }
}

DoubleGaussQuadrature::DoubleGaussQuadrature(Index degree)
    : degree_(degree), nodes_(degree), weights_(degree) {
  ARTS_ASSERT(degree % 2 == 0);
  auto gq      = GaussLegendreQuadrature(degree / 2);
  auto nodes   = gq.get_nodes();
  auto weights = gq.get_weights();

  for (Index i = 0; i < degree / 2; ++i) {
    nodes_[i]                = -0.5 + nodes[i] / 2.0;
    nodes_[degree / 2 + i]   = 0.5 + nodes[i] / 2.0;
    weights_[i]              = 0.5 * weights[i];
    weights_[degree / 2 + i] = 0.5 * weights[i];
  }
}

void LobattoQuadrature::calculate_nodes_and_weights() {
  const Index n            = degree_;
  const Index n_half_nodes = (n + 1) / 2;
  const Index n_max_iter   = 10;
  Numeric x, x_old, p_l, p_l_1, p_l_2, dp_dx, d2p_dx;

  const Index left  = (n + 1) / 2 - 1;
  const Index right = n / 2;
  Numeric n_f       = static_cast<Numeric>(n);

  for (int i = 0; i < n_half_nodes - 1; ++i) {
    //
    // Initial guess.
    //
    Numeric d_i = ((n % 2) == 0) ? 0.5 : 0.0;
    x           = sin(pi_v<Numeric> * (i + d_i) / (n_f - 0.5));

    //
    // Evaluate Legendre Polynomial and its derivative at node.
    //
    for (Index j = 0; j < n_max_iter; ++j) {
      p_l   = x;
      p_l_1 = 1.0;
      for (int l = 2; l < n; ++l) {
        // Legendre recurrence relation
        p_l_2 = p_l_1;
        p_l_1 = p_l;
        p_l   = ((2.0 * l - 1.0) * x * p_l_1 - (l - 1.0) * p_l_2) / l;
      }
      dp_dx  = (n_f - 1.0) * (x * p_l - p_l_1) / (x * x - 1.0);
      d2p_dx = (2.0 * x * dp_dx - (n_f - 1.0) * n_f * p_l) / (1.0 - x * x);

      x_old = x;

      //
      // Perform Newton step.
      //
      x       -= dp_dx / d2p_dx;
      auto dx  = x - x_old;
      if (detail::small(std::abs(dx * (x + x_old)), 1e-10)) {
        break;
      }
    }
    nodes_[right + i]   = x;
    weights_[right + i] = 2.0 / (n_f * (n_f - 1.0) * p_l * p_l);
    nodes_[left - i]    = -x;
    weights_[left - i]  = weights_[right + i];
  }
  nodes_[0]       = -1.0;
  weights_[0]     = 2.0 / (n_f * (n_f - 1.0));
  nodes_[n - 1]   = -nodes_[0];
  weights_[n - 1] = weights_[0];
}

IrregularZenithAngleGrid::IrregularZenithAngleGrid(const Vector& zenith_angles)
    : ZenithAngleGrid(zenith_angles),
      weights_(zenith_angles.size()),
      cos_theta_(zenith_angles),
      type_(QuadratureType::Trapezoidal) {
  std::transform(
      cos_theta_.begin(),
      cos_theta_.end(),
      cos_theta_.begin(),
      [](Numeric lat) { return -1.0 * cos(Conversion::deg2rad(lat)); });
  weights_ = 0.0;
  Index n  = static_cast<Index>(Vector::size());
  for (Index i = 0; i < n - 1; ++i) {
    auto dx          = 0.5 * (cos_theta_[i + 1] - cos_theta_[i]);
    weights_[i]     += dx;
    weights_[i + 1] += dx;
  }
  weights_[0]     += cos_theta_[0] + 1.0;
  weights_[n - 1] += 1.0 - cos_theta_[n - 1];
}

void ClenshawCurtisQuadrature::calculate_nodes_and_weights() {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with FFTW support.");
#else
  Index n = degree_ - 1;
  fftw_plan ifft;
  double* weights =
      reinterpret_cast<double*>(fftw_malloc(2 * (n / 2 + 1) * sizeof(double)));
  std::complex<double>* coeffs =
      reinterpret_cast<std::complex<double>*>(weights);
  Numeric n_f = static_cast<Numeric>(n);

  ifft = fftw_plan_dft_c2r_1d(static_cast<int>(n),
                              reinterpret_cast<double(*)[2]>(coeffs),
                              weights,
                              FFTW_ESTIMATE);
  // Calculate DFT input.
  for (int i = 0; i < n / 2 + 1; ++i) {
    coeffs[i] = 2.0 / (1.0 - 4.0 * i * i);
  }
  fftw_execute_dft_c2r(ifft, reinterpret_cast<double(*)[2]>(coeffs), weights);

  weights[0] *= 0.5;
  for (Index i = 0; i < n; ++i) {
    weights_[i] = weights[i] / n_f;
  }
  weights_[n] = weights[0];
  fftw_destroy_plan(ifft);
  fftw_free(weights);

  // Calculate nodes.
  for (Index i = 0; i <= n; i++) {
    nodes_[i] = -cos(pi_v<Numeric> * static_cast<Numeric>(i) / n_f);
  }
#endif
}

void FejerQuadrature::calculate_nodes_and_weights() {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with FFTW support.");
#else
  const Index n = degree_;
  fftw_plan ifft;
  double* weights =
      reinterpret_cast<double*>(fftw_malloc(2 * (n + 1) * sizeof(double)));
  std::complex<double>* coeffs =
      reinterpret_cast<std::complex<double>*>(weights);

  ifft = fftw_plan_dft_c2r_1d(static_cast<int>(n),
                              reinterpret_cast<double(*)[2]>(coeffs),
                              weights,
                              FFTW_ESTIMATE);
  // Calculate DFT input.
  for (int i = 0; i < n / 2 + 1; ++i) {
    Numeric x  = (pi_v<Numeric> * i) / static_cast<Numeric>(n);
    coeffs[i]  = std::complex<double>(cos(x), sin(x));
    coeffs[i] *= 2.0 / (1.0 - 4.0 * i * i);
  }
  fftw_execute_dft_c2r(ifft, reinterpret_cast<double(*)[2]>(coeffs), weights);
  for (Index i = 0; i < n; ++i) {
    weights_[n - i - 1] = weights[i] / static_cast<Numeric>(n);
  }

  fftw_destroy_plan(ifft);
  fftw_free(weights);

  // Calculate nodes.
  for (Index i = 0; i < n; i++) {
    nodes_[i] = -cos(pi_v<Numeric> * (static_cast<double>(i) + 0.5) /
                     static_cast<double>(n));
  }
#endif
}
}  // namespace scattering
