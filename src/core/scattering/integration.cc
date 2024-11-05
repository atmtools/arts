#include <integration.h>

#ifndef ARTS_NO_SHTNS
#include <fftw3.h>
#endif

namespace scattering {
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
