#include "legendre.h"

#include "arts_conversions.h"
#include "matpack_data.h"

namespace Legendre {
//! Converts latitude to co-latitude with a small distance from the poles and flags if this was activated
struct ColatitudeConversion {
  static constexpr Numeric dlat_lim = 1e-4;

  const bool high;
  const bool low;
  const Numeric theta;

  constexpr ColatitudeConversion(const Numeric latitude) noexcept
      : high(latitude > (90.0 - dlat_lim)),
        low(latitude < (dlat_lim - 90.0)),
        theta(high ? Conversion::deg2rad(dlat_lim)
                   : (low ? Conversion::deg2rad(180.0 - dlat_lim)
                          : Conversion::deg2rad(90.0 - latitude))) {}
};

//! Clamps the longitude in the range [-180, 180)
constexpr Numeric longitude_clamp(const Numeric lon) {
  if (lon >= -180 and lon < 180) return lon;
  if (lon < -180) return longitude_clamp(lon + 360);
  return longitude_clamp(lon - 360);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
/** Returns the Schmidt normalized Lagrange polynominal and its derivative
 *
 * The derivative is undefined if theta is 0 or pi, so the function cannot
 * be called with these values
 *
 * This function is taken from the implementation by at https://github.com/klaundal/ppigrf
 * (2024-06-13).  It implements the calculations step-by-step, updating both main and derivative
 * based on previous values.
 *
 * @param[in] theta Colatitude in radians
 * @param[in] nmax Max number of n
 * @return The pair of (nmax+1) x (nmax+1) matrices, order: main and derivative
 */

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
        P(n, n) = st * P(n - 1, m - 1);
        dP(n, n) = st * dP(n - 1, m - 1) + ct * P(n - 1, n - 1);
      } else {
        if (n == 1) {
          P(n, m) = ct * P(n - 1, m);
          dP(n, m) = st * dP(n - 1, m) - st * P(n - 1, m);

        } else {
          const Numeric Knm = static_cast<Numeric>((n - 1 + m) * (n - 1 - m)) /
                              static_cast<Numeric>((2 * n - 1) * (2 * n - 3));
          P(n, m) = ct * P(n - 1, m) - Knm * P(n - 2, m);
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

  P *= S;
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
  const auto colat = Conversion::deg2rad(90.0 - lat);
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
  Vector3 B = {0, 0, 0};
  Numeric ratn = r_ratio * r_ratio;
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

//! Computes the altitude, latitude and longitude in relation to the ellopsiod using non-iterative method
std::array<Numeric, 3> to_geodetic(const std::array<Numeric, 3> xyz,
                                   const std::array<Numeric, 2> ell) noexcept {
  using Math::pow2;
  using Math::pow3;
  using Math::pow4;

  const Numeric X = std::get<0>(xyz);
  const Numeric Y = std::get<1>(xyz);
  const Numeric Z = std::get<2>(xyz);
  const Numeric a = std::get<0>(ell);
  const Numeric e = std::get<1>(ell);
  const Numeric b = a * std::sqrt(1 - pow2(e));

  // Output
  std::array<Numeric, 3> pos;

  if (std::abs(X) > 1 / a or std::abs(Y) > 1 / a) {
    const Numeric DZ = std::sqrt(1 - pow2(e)) * Z;
    const Numeric r = std::hypot(X, Y);
    const Numeric e2p = (pow2(a) - pow2(b)) / pow2(b);
    const Numeric F = 54 * pow2(b * Z);
    const Numeric G = pow2(r) + pow2(DZ) - pow2(e) * (pow2(a) - pow2(b));
    const Numeric c = pow4(e) * F * pow2(r) / pow3(G);
    const Numeric s = std::cbrt(1 + c + std::sqrt(pow2(c) + 2 * c));
    const Numeric fP = F / (3 * pow2(G * (s + 1 / s + 1)));
    const Numeric Q = std::sqrt(1 + 2 * pow4(e) * fP);
    const Numeric r0 = (-fP * pow2(e) * r) / (1 + Q) +
                       sqrt(0.5 * pow2(a) * (1 + 1 / Q) -
                            fP * pow2(DZ) / (Q * (1 + Q)) - 0.5 * fP * pow2(r));
    const Numeric U = std::hypot(r - pow2(e) * r0, Z);
    const Numeric V = std::hypot(r - pow2(e) * r0, DZ);
    const Numeric z0 = pow2(b) * Z / (a * V);

    pos[0] = U * (1 - pow2(b) / (a * V));
    pos[1] = Conversion::atan2d(Z + e2p * z0, r);
    pos[2] = Conversion::atan2d(Y, X);
  } else if (std::abs(Z) < 1 / b) {
    pos[0] = -a;
    pos[1] = 0;
    pos[2] = 180;
  } else {
    pos[0] = std::abs(Z) - b;
    pos[1] = Z < 0 ? -90 : 90;
    pos[2] = 0;
  }

  return pos;
}
}  // namespace Legendre
