#include "legendre.h"

#include <arts_conversions.h>
#include <debug.h>

#include <algorithm>
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>

#include "fastgl.h"

namespace Legendre {
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

void schmidt(MatrixView P, MatrixView dP, const Numeric theta) {
  const Index N = P.nrows();
  if (theta < 0 or theta > Constant::pi) arts::user_error("Theta={} must be in [0, pi]", theta);
  if (N <= 1) arts::user_error("nmax={} must be > 0", N - 1);
  if (P.ncols() != N or dP.shape() != P.shape())
    arts::user_error("The Schmidt matrices must be (nmax+1) x (nmax+1), got {:B,} and {:B,}", P.shape(), dP.shape());

  const Numeric ct = std::cos(theta);
  const Numeric st = std::sin(theta);
  P                = 0;
  dP               = 0;
  P[0, 0]          = 1.0;

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
          const Numeric Knm =
              static_cast<Numeric>((n - 1 + m) * (n - 1 - m)) / static_cast<Numeric>((2 * n - 1) * (2 * n - 3));
          P[n, m]  = ct * P[n - 1, m] - Knm * P[n - 2, m];
          dP[n, m] = ct * dP[n - 1, m] - st * P[n - 1, m] - Knm * dP[n - 2, m];
        }
      }
    }
  }

  // The Schmidt normalization S of degree n, order by order: S[0, 0] = 1
  Numeric S0 = 1.0;
  for (Index n = 1; n < N; ++n) {
    S0         = S0 * (2. * n - 1) / n;
    Numeric S  = S0;
    P[n, 0]   *= S;
    dP[n, 0]  *= S;
    for (Index m = 1; m < n + 1; ++m) {
      S         = S * std::sqrt((n - m + 1) * (m == 1 ? 2.0 : 1.0) / (n + m));
      P[n, m]  *= S;
      dP[n, m] *= S;
    }
  }
}

std::pair<TriangularMatrix, TriangularMatrix> schmidt(const Numeric theta, const Index nmax) {
  if (nmax <= 0) arts::user_error("nmax={} must be > 0", nmax);

  const Index N = nmax + 1;
  Matrix      P(N, N), dP(N, N);
  schmidt(P, dP, theta);

  std::pair out{TriangularMatrix(N), TriangularMatrix(N)};
  for (Index n = 0; n < N; n++) {
    for (Index m = 0; m <= n; m++) {
      out.first[n, m]  = P[n, m];
      out.second[n, m] = dP[n, m];
    }
  }
  return out;
}

std::array<std::pair<TriangularMatrix, TriangularMatrix>, 3> dschmidt_fieldcalc(const Size    N,
                                                                                const Numeric r0,
                                                                                const Vector3 pos) try {
  const auto [r, lat, lon] = pos;

  if (lat < -90 or lat > 90) arts::user_error("Latitude is {} should be in [-90, 90]", lat);

  std::array<std::pair<TriangularMatrix, TriangularMatrix>, 3> dB{
      std::pair{TriangularMatrix(N, 0), TriangularMatrix(N, 0)},
      std::pair{TriangularMatrix(N, 0), TriangularMatrix(N, 0)},
      std::pair{TriangularMatrix(N, 0), TriangularMatrix(N, 0)}};

  // Take care of boundary issues
  const auto    colat     = Conversion::deg2rad(90.0 - lat);
  const Numeric sin_theta = std::sin(colat);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(colat, N - 1);

  // Pre-compute the cosine/sine values
  std::vector<Numeric> cosm(N);
  std::vector<Numeric> sinm(N);
  const Numeric        clon = longitude_clamp(lon);

  for (Size m = 0; m < N; ++m) {
    cosm[m] = Conversion::cosd(m * clon);
    sinm[m] = Conversion::sind(m * clon);
  }

  const Numeric r_ratio = r0 / r;
  Numeric       ratn    = r_ratio * r_ratio;
  for (Size n = 1; n < N; ++n) {
    ratn *= r_ratio;
    for (Size m = 0; m < n + 1; ++m) {
      dB[0].first[n, m]  = cosm[m] * P[n, m] * (n + 1) * ratn;
      dB[0].second[n, m] = sinm[m] * P[n, m] * (n + 1) * ratn;

      dB[1].first[n, m]  = -cosm[m] * dP[n, m] * ratn;
      dB[1].second[n, m] = -sinm[m] * dP[n, m] * ratn;

      dB[2].first[n, m]  = sinm[m] * P[n, m] * m * ratn;
      dB[2].second[n, m] = -cosm[m] * P[n, m] * m * ratn;
    }
  }

  // Fix the phi component if sin_theta is not zero
  if (std::abs(sin_theta) > 1e-6) {
    dB[2].first.flat()  /= sin_theta;
    dB[2].second.flat() /= sin_theta;
  } else {
    dB[2].first.flat()  = 0.0;
    dB[2].second.flat() = 0.0;
  }

  return dB;
}
ARTS_METHOD_ERROR_CATCH

Vector3 schmidt_fieldcalc(const ConstTriangularMatrixView& g,
                          const ConstTriangularMatrixView& h,
                          const Numeric                    r0,
                          const Vector3                    pos) {
  const auto [r, lat, lon] = pos;

  if (lat < -90 or lat > 90) arts::user_error("Latitude is {} should be in [-90, 90]", lat);
  if (g.extent(0) != h.extent(0))
    arts::user_error("The g-coefficients are {0} x {0} but the h-coefficients are {1} x {1}", g.extent(0), h.extent(0));

  const Index N = h.extent(0);

  // Take care of boundary issues
  const auto    colat     = Conversion::deg2rad(90.0 - lat);
  const Numeric sin_theta = std::sin(colat);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(colat, N - 1);

  // Pre-compute the cosine/sine values
  std::vector<Numeric> cosm(N);
  std::vector<Numeric> sinm(N);
  const Numeric        clon = longitude_clamp(lon);

  for (Index m = 0; m < N; ++m) {
    cosm[m] = Conversion::cosd(m * clon);
    sinm[m] = Conversion::sind(m * clon);
  }

  const Numeric r_ratio = r0 / r;
  Vector3       B       = {0, 0, 0};
  Numeric       ratn    = r_ratio * r_ratio;
  for (Index n = 1; n < N; ++n) {
    ratn *= r_ratio;

    // The h with m = 0 do not affect the field, and the m = 0 terms not the longitudinal component
    B[0] += g[n, 0] * P[n, 0] * (n + 1) * ratn;
    B[1] -= g[n, 0] * dP[n, 0] * ratn;

    for (Index m = 1; m < n + 1; ++m) {
      B[0] += (g[n, m] * cosm[m] + h[n, m] * sinm[m]) * P[n, m] * (n + 1) * ratn;
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

  if (x < -1 or x > 1) arts::user_error("x={} not in [-1, 1]", x);

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

void legendre_polynomials(VectorView p, Numeric x) {
  using boost::math::legendre_next;

  if (x < -1 or x > 1) arts::user_error("x={} not in [-1, 1]", x);

  const Index n = p.size();
  if (n > 0) p[0] = 1.0;
  if (n > 1) p[1] = x;
  for (Index l = 1; l + 1 < n; l++) p[l + 1] = legendre_next(static_cast<unsigned>(l), x, p[l], p[l - 1]);
}

void wigner_d(VectorView d, Index m, Index n, Numeric theta) {
  const auto  L  = static_cast<Index>(d.size());
  const Index l0 = std::max(std::abs(m), std::abs(n));
  d              = 0.0;
  if (l0 >= L) return;

  /* d^M_MN, M >= |N|, as sqrt(binomial(2 M, M + N)) cos^(M + N)(theta / 2) sin^(M - N)(theta / 2), the root first
     and then the powers, which only decrease it.  The binomial coefficient is below 4^M, which a double holds to
     M = 511; above, the root is the product of the sqrt((M - N + k) / k), k = 1, ..., M + N, each with a cosine. */
  const Numeric c      = std::cos(0.5 * theta);
  const Numeric s      = std::sin(0.5 * theta);
  const auto    parity = [](Index k) { return k % 2 == 0 ? 1.0 : -1.0; };
  const auto    first  = [&](Index M, Index N) {
    Numeric v = 1.0;
    if (M <= 511) {
      for (Index k = 1; k <= M + N; k++) v *= static_cast<Numeric>(M - N + k) / static_cast<Numeric>(k);
      v = std::sqrt(v);
      for (Index k = 0; k < M + N; k++) v *= c;
    } else {
      for (Index k = 1; k <= M + N; k++) v *= c * std::sqrt(static_cast<Numeric>(M - N + k) / static_cast<Numeric>(k));
    }
    for (Index k = 0; k < M - N; k++) v *= s;
    return parity(M - N) * v;
  };
  if (std::abs(m) >= std::abs(n))
    d[l0] = m >= 0 ? first(m, n) : parity(m - n) * first(-m, -n);
  else
    d[l0] = n > 0 ? parity(m - n) * first(n, m) : first(-n, -m);

  const Numeric x     = std::cos(theta);
  Index         start = l0;
  if (l0 == 0) {
    if (L > 1) d[1] = x;
    start = 1;
  }
  // The root of the d^(l+1) term of degree l is that of the d^(l-1) term of degree l + 1
  const auto dm = static_cast<Numeric>(m), dn = static_cast<Numeric>(n);
  const auto root   = [dm, dn](Numeric k) { return std::sqrt((k * k - dm * dm) * (k * k - dn * dn)); };
  Numeric    root_l = root(static_cast<Numeric>(start));
  for (Index l = start; l + 1 < L; l++) {
    const auto    k         = static_cast<Numeric>(l);
    const Numeric root_next = root(k + 1);
    const Numeric b         = (2 * k + 1) * (k * (k + 1) * x - dm * dn);
    d[l + 1]                = (b * d[l] - (k + 1) * root_l * d[l - 1]) / (k * root_next);
    root_l                  = root_next;
  }
}

Numeric factorial(Index i) { return boost::math::factorial<Numeric>(static_cast<unsigned>(i)); }

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
