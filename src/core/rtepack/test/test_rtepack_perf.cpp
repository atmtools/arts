#include <omp.h>
#include <rng.h>
#include <rtepack.h>
#include <time_report.h>

#include <iostream>
#include <vector>

namespace {
Numeric test_expm(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_expm");

  Numeric sum = 0.0;

  for (auto& k : K) {
    sum += rtepack::tran(k, k, 1.0)()[0, 0];
  }

  return sum;
}

Numeric test_expm_deriv(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_expm_deriv");

  Numeric sum = 0.0;

  for (auto& k : K) {
    const auto t  = rtepack::tran(k, k, 1.0);
    sum          += t.deriv(t(), k, k, k, 1.0, 0.0)[0, 0];
  }

  return sum;
}

Numeric test_linsrc(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_linsrc");

  Numeric sum = 0.0;

  for (auto& k : K) {
    sum += rtepack::tran(k, k, 1.0).linsrc()[0, 0];
  }

  return sum;
}

Numeric test_linsrc_deriv(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_linsrc_deriv");

  Numeric sum = 0.0;

  for (auto& k : K) {
    sum += rtepack::tran(k, k, 1.0).linsrc_deriv(k, 1.0, 0.0)[0, 0];
  }

  return sum;
}

Numeric test_linprop(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_linprop");

  Numeric sum = 0.0;

  for (auto& k : K) {
    const auto t  = rtepack::tran(k, k, 1.0);
    sum          += t.linsrc_linprop(t(), k, k, 1.0)[0, 0];
  }

  return sum;
}

Numeric test_linprop_deriv_k1(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_linprop_deriv_k1");

  Numeric sum = 0.0;

  for (auto& k : K) {
    const auto t     = rtepack::tran(k, k, 1.0);
    const auto exp_t = t();
    const auto la    = t.linsrc_linprop(exp_t, k, k, 1.0)[0, 0];
    sum +=
        t.linsrc_linprop_deriv(la, exp_t, k, k, k, exp_t, 1.0, 1.0, true)[0, 0];
  }

  return sum;
}

Numeric test_linprop_deriv_k2(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_linprop_deriv_k2");

  Numeric sum = 0.0;

  for (auto& k : K) {
    const auto t      = rtepack::tran(k, k, 1.0);
    const auto exp_t  = t();
    const auto la     = t.linsrc_linprop(exp_t, k, k, 1.0)[0, 0];
    sum              += t.linsrc_linprop_deriv(
        la, exp_t, k, k, k, exp_t, 1.0, 1.0, false)[0, 0];
  }

  return sum;
}

Numeric test_logk(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_logk");

  Numeric sum = 0.0;

  for (auto& k : K) {
    sum += rtepack::logK(rtepack::tran(k, k, 1.0)()).A();
  }

  return sum;
}

Numeric test_sqrt(const PropmatVector& K) {
  ARTS_NAMED_TIME_REPORT("test_sqrt");

  Numeric sum = 0.0;

  for (auto& k : K) {
    sum += rtepack::sqrt(k)[0, 0].real();
  }

  return sum;
}

Numeric test_transmittance_matrix_init_constant(
    const std::vector<PropmatVector>& K,
    const std::vector<PropmatMatrix>& dK,
    const Vector& r,
    const Tensor3& dr) {
  ARTS_NAMED_TIME_REPORT(
      std::format("test_transmittance_matrix_init_constant; threads {}",
                  omp_get_max_threads()));

  Numeric sum{};
  rtepack::TransmittanceMatrix tm;
  tm.init(K, dK, r, dr, TransmittanceOption::constant);
  sum += tm.T[0, 0][0, 0];

  return sum;
}

Numeric test_transmittance_matrix_init_linsrc(
    const std::vector<PropmatVector>& K,
    const std::vector<PropmatMatrix>& dK,
    const Vector& r,
    const Tensor3& dr) {
  ARTS_NAMED_TIME_REPORT(
      std::format("test_transmittance_matrix_init_linsrc; threads {}",
                  omp_get_max_threads()));

  Numeric sum{};
  rtepack::TransmittanceMatrix tm;
  tm.init(K, dK, r, dr, TransmittanceOption::linsrc);
  sum += tm.T[0, 0][0, 0];

  return sum;
}

Numeric test_transmittance_matrix_init_linprop(
    const std::vector<PropmatVector>& K,
    const std::vector<PropmatMatrix>& dK,
    const Vector& r,
    const Tensor3& dr) {
  ARTS_NAMED_TIME_REPORT(
      std::format("test_transmittance_matrix_init_linprop; threads {}",
                  omp_get_max_threads()));

  Numeric sum{};
  rtepack::TransmittanceMatrix tm;
  tm.init(K, dK, r, dr, TransmittanceOption::linprop);
  sum += tm.T[0, 0][0, 0];

  return sum;
}
}  // namespace

int main() {
  Numeric buf = 0.0;

  {
    constexpr Index M = 10'000'000;
    PropmatVector K(M);
    MatrixView Kv{MatrixView::base{reinterpret_cast<Numeric*>(K.data_handle()),
                                   std::array<Index, 2>{M, 7}}};
    random_numbers(Kv, 0.0, 1.0);
    const Vector r = random_numbers<1>(10);

    buf += test_expm(K);
    buf += test_logk(K);
    buf += test_sqrt(K);
    buf += test_linsrc(K);
    buf += test_linprop(K);
    buf += test_expm_deriv(K);
    buf += test_linsrc_deriv(K);
    buf += test_linprop_deriv_k1(K);
    buf += test_linprop_deriv_k2(K);
  }

  {
    constexpr Index M = 10'000;
    constexpr Index N = 100;
    constexpr Index P = 10;

    std::vector<PropmatVector> K(M);
    for (auto& kv : K) {
      kv = PropmatVector(N);
      MatrixView Kv{
          MatrixView::base{reinterpret_cast<Numeric*>(kv.data_handle()),
                           std::array<Index, 2>{N, 7}}};
      random_numbers(Kv, 0.0, 1.0);
    }
    const Vector r = random_numbers<1>({M - 1});

    std::vector<PropmatMatrix> dK(M);
    for (auto& kv : dK) {
      kv = PropmatMatrix(P, N);
      Tensor3View dKv{
          Tensor3View::base{reinterpret_cast<Numeric*>(kv.data_handle()),
                            std::array<Index, 3>{N, P, 7}}};
      random_numbers(dKv, 0.0, 1.0);
    }
    const Tensor3 dr = random_numbers<3>({2, M - 1, P});

    buf += test_transmittance_matrix_init_constant(K, dK, r, dr);
    buf += test_transmittance_matrix_init_linsrc(K, dK, r, dr);
    buf += test_transmittance_matrix_init_linprop(K, dK, r, dr);

    const int x = omp_get_max_threads();
    omp_set_num_threads(1);

    buf += test_transmittance_matrix_init_constant(K, dK, r, dr);
    buf += test_transmittance_matrix_init_linsrc(K, dK, r, dr);
    buf += test_transmittance_matrix_init_linprop(K, dK, r, dr);

    omp_set_num_threads(x);
  }

  std::println(std::cerr, "Prevent optimizing away: {}", buf);
  arts::print_report();
}