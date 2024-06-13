#include "legendre.h"

#include <boost/math/special_functions/legendre.hpp>
#include <cmath>
#include <iostream>

#include "arts_conversions.h"
#include "debug.h"
#include "fastgl.h"

namespace Legendre {
//! Clamps the longitude in the range [-180, 180)
constexpr Numeric longitude_clamp(Numeric lon) {
  while (lon <= -180) lon += 360;
  while (lon > 180) lon -= 360;
  return lon;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

std::pair<Matrix, Matrix> schmidt(const Numeric theta, const Index nmax) {
  ARTS_USER_ERROR_IF(
      theta < 0 or theta > Constant::pi, "Theta=", theta, " must be in [0, pi]")
  ARTS_USER_ERROR_IF(nmax <= 0, "nmax=", nmax, " must be > 0")

  const Index N = 1 + nmax;

  const Numeric ct = std::cos(theta);
  const Numeric st = std::sin(theta);
  Matrix P(N, N, 0);
  Matrix dP(N, N, 0);
  Matrix S(N, N, 0);
  S(0, 0) = 1.0;
  P(0, 0) = 1.0;

  for (Index n = 1; n < N; ++n) {
    for (Index m = 0; m < n + 1; ++m) {
      if (n == m) {
        P(n, n)  = st * P(n - 1, m - 1);
        dP(n, n) = st * dP(n - 1, m - 1) + ct * P(n - 1, n - 1);
      } else {
        if (n == 1) {
          P(n, m)  = ct * P(n - 1, m);
          dP(n, m) = st * dP(n - 1, m) - st * P(n - 1, m);

        } else {
          const Numeric Knm = static_cast<Numeric>((n - 1 + m) * (n - 1 - m)) /
                              static_cast<Numeric>((2 * n - 1) * (2 * n - 3));
          P(n, m)  = ct * P(n - 1, m) - Knm * P(n - 2, m);
          dP(n, m) = ct * dP(n - 1, m) - st * P(n - 1, m) - Knm * dP(n - 2, m);
        }
      }

      // compute Schmidt normalization
      if (m == 0)
        S(n, 0) = S(n - 1, 0) * (2. * n - 1) / n;
      else
        S(n, m) =
            S(n, m - 1) * std::sqrt((n - m + 1) * (int(m == 1) + 1.) / (n + m));
    }
  }

  P  *= S;
  dP *= S;

  return {P, dP};
}

Vector3 schmidt_fieldcalc(const Matrix& g,
                          const Matrix& h,
                          const Numeric r0,
                          const Vector3 pos) {
  const auto [r, lat, lon] = pos;

  ARTS_USER_ERROR_IF(
      lat < -90 and lat > 90, "Latitude is ", lat, " should be in [-90, 90]")

  const Index N = h.nrows();

  // Take care of boundary issues
  const auto colat        = Conversion::deg2rad(90.0 - lat);
  const Numeric sin_theta = std::sin(colat);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(colat, N - 1);

  // Pre-compute the cosine/sine values
  std::vector<Numeric> cosm(N);
  std::vector<Numeric> sinm(N);
  const Numeric clon = longitude_clamp(lon);

  for (Index m = 0; m < N; ++m) {
    cosm[m] = Conversion::cosd(m * clon);
    sinm[m] = Conversion::sind(m * clon);
  }

  const Numeric r_ratio = r0 / r;
  Vector3 B             = {0, 0, 0};
  Numeric ratn          = r_ratio * r_ratio;
  for (Index n = 1; n < N; ++n) {
    ratn *= r_ratio;
    for (Index m = 0; m < n + 1; ++m) {
      B[0] +=
          (g(n, m) * cosm[m] + h(n, m) * sinm[m]) * P(n, m) * (n + 1) * ratn;
      B[1] -= (g(n, m) * cosm[m] + h(n, m) * sinm[m]) * dP(n, m) * ratn;
      B[2] += (g(n, m) * sinm[m] - h(n, m) * cosm[m]) * P(n, m) * m * ratn;
    }
  }

  // Fix the phi component if sin_theta is not zero
  if (std::abs(sin_theta) > 1e-6) {
    B[2] /= sin_theta;
  } else {
    B[2] = 0.0;
  }

  return B;
}
#pragma GCC diagnostic pop

constexpr Numeric next(Numeric p0, Numeric p1, Numeric ni, Numeric x) {
  return ((2.0 * ni + 1.0) * x * p1 - ni * p0) / (ni + 1.0);
}

Numeric legendre_sum(const ExhaustiveConstVectorView& s, const Numeric& x) {
  using boost::math::legendre_next;

  ARTS_USER_ERROR_IF(x < -1 or x > 1, "x=", x, " must be in [-1, 1]")

  const Index n = s.nelem();
  if (n == 0) return 0.0;

  Numeric p0  = 1.0;
  Numeric p1  = x;
  Numeric out = s[0] * p0;

  for (Index i = 1; i < n; i++) {
    out += s[i] * p1;
    std::swap(p0, p1);
    p1 = legendre_next(static_cast<int>(i), x, p0, p1);
  }

  return out;
}

Numeric legendre(Index n, Numeric x) {
  using boost::math::legendre_p;
  return legendre_p(static_cast<int>(n), x);
}

Numeric factorial(Index i) {
  return boost::math::factorial<Numeric>(static_cast<unsigned>(i));
}

//! port of boost tgamma_ratio_imp
Numeric tgamma_ratio(Numeric x, Numeric y) {
  using boost::math::tgamma_ratio;
  return tgamma_ratio(x, y);
}

// Port of Boost.Math legendre_p_imp
Numeric assoc_legendre(Index l, Index m, Numeric x) {
  using boost::math::legendre_p;
  return legendre_p(static_cast<int>(l), static_cast<int>(m), x);
}

void PositiveDoubleGaussLegendre(ExhaustiveVectorView x,
                                 ExhaustiveVectorView w) {
  const Index n = x.size();
  ARTS_ASSERT(n == w.size());  // same size

  for (Index k = 0; k < n; k++) {
    auto p = fastgl::GLPair(n, n - k);
    x[k]   = 0.5 * (1.0 + p.x());
    w[k]   = 0.5 * p.weight;
  }
}

void PositiveGaussLegendre(ExhaustiveVectorView x, ExhaustiveVectorView w) {
  const Index n = x.size();
  ARTS_ASSERT(n == w.size());

  for (Index k = 0; k < n; k++) {
    auto p = fastgl::GLPair(2 * n, 2 * n - k - n);
    x[k]   = p.x();
    w[k]   = p.weight;
  }
}

void GaussLegendre(ExhaustiveVectorView x, ExhaustiveVectorView w) {
  const Index n = x.size();
  ARTS_ASSERT(n == w.size());

  for (Index k = 0; k < n; k++) {
    auto p = fastgl::GLPair(n, n - k);
    x[k]   = p.x();
    w[k]   = p.weight;
  }
}
}  // namespace Legendre
