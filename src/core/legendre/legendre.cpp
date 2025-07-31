#include "legendre.h"

#include <arts_conversions.h>
#include <debug.h>

#include <boost/math/special_functions/legendre.hpp>
#include <cmath>

#include "fastgl.h"

namespace Legendre {
Numeric SphericalField::total() const noexcept { return std::hypot(U, S, E); }

Numeric SphericalField::total_horizontal() const noexcept {
  return std::hypot(S, E);
}

namespace {
//! Clamps the longitude in the range [-180, 180)
constexpr Numeric longitude_clamp(Numeric lon) {
  while (lon <= -180) lon += 360;
  while (lon > 180) lon -= 360;
  return lon;
}
}  // namespace

#ifndef _MSC_VER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

std::pair<Matrix, Matrix> schmidt(const Numeric theta, const Index nmax) {
  ARTS_USER_ERROR_IF(
      theta < 0 or theta > Constant::pi, "Theta={} must be in [0, pi]", theta)
  ARTS_USER_ERROR_IF(nmax <= 0, "nmax={} must be > 0", nmax)

  const Index N = 1 + nmax;

  const Numeric ct = std::cos(theta);
  const Numeric st = std::sin(theta);
  Matrix P(N, N, 0);
  Matrix dP(N, N, 0);
  Matrix S(N, N, 0);
  S[0, 0] = 1.0;
  P[0, 0] = 1.0;

  for (Index n = 1; n < N; ++n) {
    for (Index m = 0; m < n + 1; ++m) {
      if (n == m) {
        P[n, n]  = st * P[n - 1, m - 1];
        dP[n, n] = st * dP[n - 1, m - 1] + ct * P[n - 1, n - 1];
      } else {
        if (n == 1) {
          P[n, m]  = ct * P[n - 1, m];
          dP[n, m] = st * dP[n - 1, m] - st * P[n - 1, m];

        } else {
          const Numeric Knm = static_cast<Numeric>((n - 1 + m) * (n - 1 - m)) /
                              static_cast<Numeric>((2 * n - 1) * (2 * n - 3));
          P[n, m]  = ct * P[n - 1, m] - Knm * P[n - 2, m];
          dP[n, m] = ct * dP[n - 1, m] - st * P[n - 1, m] - Knm * dP[n - 2, m];
        }
      }

      // compute Schmidt normalization
      if (m == 0)
        S[n, 0] = S[n - 1, 0] * (2. * n - 1) / n;
      else
        S[n, m] =
            S[n, m - 1] * std::sqrt((n - m + 1) * (int(m == 1) + 1.) / (n + m));
    }
  }

  P  *= S;
  dP *= S;

  return {P, dP};
}

Tensor4 dschmidt_fieldcalc(const Size N, const Numeric r0, const Vector3 pos) {
  const auto [r, lat, lon] = pos;

  ARTS_USER_ERROR_IF(
      lat < -90 and lat > 90, "Latitude is {} should be in [-90, 90]", lat)

  Tensor4 dB(2, 3, N, N, 0.0);

  // Take care of boundary issues
  const auto colat        = Conversion::deg2rad(90.0 - lat);
  const Numeric sin_theta = std::sin(colat);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(colat, N - 1);

  // Pre-compute the cosine/sine values
  std::vector<Numeric> cosm(N);
  std::vector<Numeric> sinm(N);
  const Numeric clon = longitude_clamp(lon);

  for (Size m = 0; m < N; ++m) {
    cosm[m] = Conversion::cosd(m * clon);
    sinm[m] = Conversion::sind(m * clon);
  }

  const Numeric r_ratio = r0 / r;
  Numeric ratn          = r_ratio * r_ratio;
  for (Size n = 1; n < N; ++n) {
    ratn *= r_ratio;
    for (Size m = 0; m < n + 1; ++m) {
      dB[0, 0, n, m] = cosm[m] * P[n, m] * (n + 1) * ratn;
      dB[1, 0, n, m] = sinm[m] * P[n, m] * (n + 1) * ratn;

      dB[0, 1, n, m] = -cosm[m] * dP[n, m] * ratn;
      dB[1, 1, n, m] = -sinm[m] * dP[n, m] * ratn;

      dB[0, 2, n, m] = sinm[m] * P[n, m] * m * ratn;
      dB[1, 2, n, m] = -cosm[m] * P[n, m] * m * ratn;
    }
  }

  // Fix the phi component if sin_theta is not zero
  if (std::abs(sin_theta) > 1e-6) {
    dB[joker, 2] /= sin_theta;
  } else {
    dB[joker, 2] = 0.0;
  }

  return dB;
}

Vector3 schmidt_fieldcalc(const ConstMatrixView& g,
                          const ConstMatrixView& h,
                          const Numeric r0,
                          const Vector3 pos) {
  const auto [r, lat, lon] = pos;

  ARTS_USER_ERROR_IF(
      lat < -90 and lat > 90, "Latitude is {} should be in [-90, 90]", lat)

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
          (g[n, m] * cosm[m] + h[n, m] * sinm[m]) * P[n, m] * (n + 1) * ratn;
      B[1] -= (g[n, m] * cosm[m] + h[n, m] * sinm[m]) * dP[n, m] * ratn;
      B[2] += (g[n, m] * sinm[m] - h[n, m] * cosm[m]) * P[n, m] * m * ratn;
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
#ifndef _MSC_VER
#pragma GCC diagnostic pop
#endif

Numeric legendre_sum(const ConstVectorView& s, const Numeric& x) {
  using boost::math::legendre_next;

  ARTS_USER_ERROR_IF(x < -1 or x > 1, "x={} not in [-1, 1]", x)

  const Index n = s.size();
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

void PositiveDoubleGaussLegendre(VectorView x, VectorView w) {
  const Size n = x.size();
  assert(n == w.size());  // same size

  for (Size k = 0; k < n; k++) {
    auto p = fastgl::GLPair(n, n - k);
    x[k]   = 0.5 * (1.0 + p.x());
    w[k]   = 0.5 * p.weight;
  }
}

void PositiveGaussLegendre(VectorView x, VectorView w) {
  const Size n = x.size();
  assert(n == w.size());

  for (Size k = 0; k < n; k++) {
    auto p = fastgl::GLPair(2 * n, 2 * n - k - n);
    x[k]   = p.x();
    w[k]   = p.weight;
  }
}

void GaussLegendre(VectorView x, VectorView w) {
  const Size n = x.size();
  assert(n == w.size());

  for (Size k = 0; k < n; k++) {
    auto p = fastgl::GLPair(n, n - k);
    x[k]   = p.x();
    w[k]   = p.weight;
  }
}
}  // namespace Legendre
