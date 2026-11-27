#include <arts_constants.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ranges>
#include <stdexcept>

#include "debug.h"
#include "legendre.h"

void sumlegtest() {
  const Vector c = {1.,          3.,          5.,          7.,         9.,          11.,         13.,
                    15.,         17.,         19.,         21.,        23.,         25.,         27.,
                    29.,         31.,         33.,         32.8125,    29.91796875, 25.96508789, 21.84068298,
                    17.98672199, 14.59263057, 11.71004179, 9.31995848, 7.37111478,  5.80115833,  4.54775559,
                    3.55391083,  2.77005711,  2.15442342,  1.67254542};

  const Vector nu{0.01985507, 0.10166676, 0.2372338, 0.40828268, 0.59171732, 0.7627662, 0.89833324, 0.98014493};

  // From numpy.polynomial.legendre.Legendre
  const Vector res{0.048454, 0.38029738, -0.52891716, 0.44314246, -1.33650515, 1.33803551, 8.08787746, -29.07110553};

  Vector sum(8);
  for (Index i = 0; i < 8; i++) { sum[i] = Legendre::legendre_sum(c, nu[i]); }

  const Numeric sum_squared_diff = std::transform_reduce(
      res.begin(), res.end(), sum.begin(), 0., std::plus<>(), [](auto a, auto b) { return (a - b) * (a - b); });

  if (sum_squared_diff > 1e-10) { throw std::runtime_error("Legendre sum test failed"); }
}

//! legendre_polynomials against boost's single-order legendre, and the closed form P_l(1) = 1, P_l(-1) = (-1)^l
void polynomialstest() {
  Vector p(40);
  for (Numeric x : {-1.0, -0.7, -0.1, 0.0, 0.3, 0.95, 1.0}) {
    Legendre::legendre_polynomials(p, x);
    for (Index l = 0; l < static_cast<Index>(p.size()); l++) {
      const Numeric ref = Legendre::legendre(l, x);
      if (std::abs(p[l] - ref) > 1e-13) throw std::runtime_error("Legendre polynomials test failed");
      if (std::abs(x) == 1.0 and p[l] != (x < 0 and l % 2 == 1 ? -1.0 : 1.0))
        throw std::runtime_error("Legendre polynomials test failed at x = +-1");
    }
  }
}

//! Wigner's d^l_mn(theta) by its explicit sum, accurate for the small degrees it is used for here
Numeric wigner_d_sum(Index l, Index m, Index n, Numeric theta) {
  if (std::abs(m) > l or std::abs(n) > l) return 0.0;
  const auto    f = [](Index k) { return Legendre::factorial(k); };
  const Numeric c = std::cos(theta / 2), s = std::sin(theta / 2);
  Numeric       sum = 0.0;
  for (Index k = std::max<Index>(0, n - m); k <= std::min(l + n, l - m); k++)
    sum += ((m - n + k) % 2 == 0 ? 1.0 : -1.0) * std::pow(c, static_cast<Numeric>(2 * l + n - m - 2 * k)) *
           std::pow(s, static_cast<Numeric>(m - n + 2 * k)) / (f(l + n - k) * f(k) * f(m - n + k) * f(l - m - k));
  return std::sqrt(f(l + m) * f(l - m) * f(l + n) * f(l - n)) * sum;
}

/* wigner_d against
   - Wigner's explicit sum, every order |m|, |n| <= 5 to degree 12, relative to the largest value of the order pair
     and, near the poles, where the functions vanish as powers of the angle, relative to each value;
   - the Legendre polynomials, d^l_00 = P_l;
   - the Schmidt functions to degree 150, d^l_0m = S^m_l / sqrt(2) for m > 0. */
void wignerdtest() {
  constexpr Index L = 13;
  Vector          d(L);
  for (Numeric theta : {0.0, 1e-6, 0.3, 1.1, Constant::pi / 2, 2.5, Constant::pi - 1e-6, Constant::pi}) {
    const bool pole = std::sin(theta) < 1e-3;
    for (Index m = -5; m <= 5; m++) {
      for (Index n = -5; n <= 5; n++) {
        Legendre::wigner_d(d, m, n, theta);
        Numeric scale = 0.0;
        for (Index l = 0; l < L; l++) scale = std::max(scale, std::abs(wigner_d_sum(l, m, n, theta)));
        for (Index l = 0; l < L; l++) {
          const Numeric ref   = wigner_d_sum(l, m, n, theta);
          const Numeric error = std::abs(d[l] - ref);
          if (error > 1e-13 * scale or (pole and error > 1e-12 * std::abs(ref)))
            throw std::runtime_error(
                arts::format("wigner_d({}, {}, {}, {}) = {}, but Wigner's sum is {}", l, m, n, theta, d[l], ref));
        }
      }
    }
  }

  Vector p(40), d00(40);
  for (Numeric x : {-1.0, -0.7, 0.0, 0.3, 1.0}) {
    Legendre::legendre_polynomials(p, x);
    Legendre::wigner_d(d00, 0, 0, std::acos(x));
    for (Index l = 0; l < 40; l++)
      if (std::abs(d00[l] - p[l]) > 1e-13) throw std::runtime_error("wigner_d(l, 0, 0) is not P_l");
  }

  constexpr Index N = 151;
  Matrix          S(N, N), dS(N, N);
  Vector          order(N);
  for (Numeric theta : {0.01, 0.7, 1.9, 3.1}) {
    Legendre::schmidt(S, dS, theta);
    for (Index m = 0; m < N; m++) {
      Legendre::wigner_d(order, 0, m, theta);
      const Numeric scale = m == 0 ? 1.0 : Constant::inv_sqrt_2;
      for (Index l = m; l < N; l++)
        if (std::abs(order[l] - scale * S[l, m]) > 1e-12)
          throw std::runtime_error(
              arts::format("wigner_d({}, 0, {}, {}) = {} is not the Schmidt function {} over sqrt(2)",
                           l,
                           m,
                           theta,
                           order[l],
                           S[l, m]));
    }
  }
}

int main() {
  sumlegtest();
  polynomialstest();
  wignerdtest();
}
