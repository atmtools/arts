#include "legendre.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "debug.h"

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
constexpr Numeric longitude_clamp(Numeric lon) {
  while (lon <= -180) lon += 360;
  while (lon > 180) lon -= 360;
  return lon;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
/** Returns the Schmidt normalized Lagrange polynominal and its derivative
 *
 * The derivative is undefined if theta is 0 or pi, so the function cannot
 * be called with these values
 *
 * This function is taken from the implementation by Isabela de Oliveira Martins
 * at https://github.com/de-oliveira/IsabelaFunctions/blob/master/IsabelaFunctions/fieldmodel.py
 * (2021-05-06).  It implements the calculations step-by-step, updating both main and derivative
 * based on previous values.
 *
 * @param[in] theta Colatitude in radians
 * @param[in] nmax Max number of n
 * @return The pair of (nmax+1) x (nmax+1) matrices, order: main and derivative
 */
std::pair<Matrix, Matrix> schmidt(const Numeric theta,
                                  const Index nmax) ARTS_NOEXCEPT {
  ARTS_ASSERT(theta > 0 and theta < Constant::pi)
  ARTS_ASSERT(nmax > 0)

  using Math::pow2;
  using std::sqrt;

  Index N = 1 + nmax;

  const Numeric x = std::cos(theta);
  Matrix P(N, N, 0);
  Matrix dP(N, N, 0);

  // Start values
  P(0, 0) = 1.0;
  P(1, 0) = x;
  P(1, 1) = -sqrt(1.0 - pow2(x));
  dP(0, 0) = 0.0;
  dP(1, 0) = 1.0;

  for (Index n = 2; n < N; ++n) {
    P(n, 0) = ((2 * (n - 1) + 1) * x * P(n - 1, 0) - (n - 1) * P(n - 2, 0)) / n;
  }

  dP(nmax, 0) = nmax / (pow2(x) - 1) * (x * P(nmax, 0) - P(nmax - 1, 0));
  for (Index n = 2; n < nmax; ++n) {
    dP(n, 0) = (n + 1) / (pow2(x) - 1) * (P(n + 1, 0) - x * P(n, 0));
  }

  Numeric Cm = sqrt(2.0);
  for (Index m = 1; m < N; ++m) {
    Cm /= sqrt(Numeric(2 * m * (2 * m - 1)));
    P(m, m) = std::pow(1 - pow2(x), 0.5 * m) * Cm;

    for (Index i = 1; i < m; ++i) {
      P(m, m) *= 2 * i + 1;
    }

    dP(m, m) = -P(m, m) * m * x / sqrt(1 - pow2(x));

    if (nmax > m) {
      Numeric twoago = 0.0;
      for (Index n = m + 1; n < N; ++n) {
        P(n, m) = (x * (2 * n - 1) * P(n - 1, m) -
                   sqrt(Numeric((n + m - 1) * (n - m - 1))) * twoago) /
                  sqrt(Numeric(pow2(n) - pow2(m)));
        twoago = P(n - 1, m);
      }
    }
  }

  for (Index n = 2; n < N; ++n) {
    for (Index m = 1; m < n; ++m) {
      dP(n, m) = sqrt(Numeric((n - m) * (n + m + 1))) * P(n, m + 1) -
                 P(n, m) * m * x / sqrt(1 - pow2(x));
    }
  }

  return {P, dP};
}

SphericalField schmidt_fieldcalc(const Matrix& g,
                                 const Matrix& h,
                                 const Numeric r0,
                                 const Numeric r,
                                 const Numeric lat,
                                 const Numeric lon) ARTS_NOEXCEPT {
  ARTS_ASSERT(
      lat >= -90 and lat <= 90, "Latitude is ", lat, " should be in [-90, 90]")

  const Index N = h.nrows();

  // Take care of boundary issues
  const auto [low_lat, high_lat, effective_theta] = ColatitudeConversion(lat);
  const Numeric sin_theta = std::sin(effective_theta);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(effective_theta, N - 1);

  // Pre-compute the cosine/sine values
  std::vector<Numeric> cosm(N);
  std::vector<Numeric> sinm(N);
  const Numeric clon = longitude_clamp(lon);
  if (low_lat or high_lat) {
    for (Index m = 0; m < N; ++m) {
      cosm[m] = 1.0;
      sinm[m] = 0.0;
    }
  } else {
    for (Index m = 0; m < N; ++m) {
      cosm[m] = Conversion::cosd(m * clon);
      sinm[m] = Conversion::sind(m * clon);
    }
  }

  // Compute the field
  SphericalField F;  // Auto init to zeroes
  for (Index n = 1; n < N; ++n) {
    const Numeric ratn = std::pow(r0 / r, n + 2);
    for (Index m = 0; m < n + 1; ++m) {
      F.U += (g(n, m) * cosm[m] + h(n, m) * sinm[m]) * P(n, m) * (n + 1) * ratn;
      F.S += (g(n, m) * cosm[m] + h(n, m) * sinm[m]) * dP(n, m) * ratn;
      F.E += (g(n, m) * sinm[m] + h(n, m) * cosm[m]) * P(n, m) * m * ratn;
    }
  }
  F.S *= sin_theta;
  F.E /= sin_theta;

  // Adjust by theta (nb. can be nan when undefined, but set to zero instead)
  if (high_lat or low_lat) {
    F.S = F.total_horizontal();
    F.E = 0;
  }

  return F;
}

MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g,
                                         const Matrix& h,
                                         const Numeric r0,
                                         const Vector& rv,
                                         const Numeric lat,
                                         const Vector& lonv) ARTS_NOEXCEPT {
  ARTS_ASSERT(
      lat >= -90 and lat <= 90, "Latitude is ", lat, " should be in [-90, 90]")

  const Index N = h.nrows();
  const Index nlon = lonv.nelem();
  const Index nr = rv.nelem();

  // Take care of boundary issues
  const auto [low_lat, high_lat, effective_theta] = ColatitudeConversion(lat);
  const Numeric sin_theta = std::sin(effective_theta);

  // Compute the legendre polynominal with Schmidt renormalization
  const auto [P, dP] = schmidt(effective_theta, N - 1);

  // Pre-compute the radius-ratio-power
  Matrix ratnm(nr, N);
  for (Index ir = 0; ir < nr; ir++) {
    ARTS_ASSERT(rv[ir] > 0, "Radius are [", rv, "] must be in (0, inf]")
    for (Index n = 1; n < N; ++n) {
      ratnm(ir, n) = std::pow(r0 / rv[ir], n + 2);
    }
  }

  MatrixOfSphericalField out(nlon, nr);  // Auto init to zeroes
  for (Index ilon = 0; ilon < nlon; ilon++) {
    // Pre-compute the cosine/sine values
    std::vector<Numeric> cosm(N);
    std::vector<Numeric> sinm(N);
    const Numeric clon = longitude_clamp(lonv[ilon]);
    if (low_lat or high_lat) {
      for (Index m = 0; m < N; ++m) {
        cosm[m] = 1.0;
        sinm[m] = 0.0;
      }
    } else {
      for (Index m = 0; m < N; ++m) {
        cosm[m] = Conversion::cosd(m * clon);
        sinm[m] = Conversion::sind(m * clon);
      }
    }

    // Compute the field
    for (Index ir = 0; ir < nr; ir++) {
      SphericalField& F = out(ilon, ir);  // Auto init'd to zeroes
      for (Index n = 1; n < N; ++n) {
        const Numeric ratn = ratnm(ir, n);
        for (Index m = 0; m < n + 1; ++m) {
          const Numeric gnm{g(n, m)};
          const Numeric hnm{h(n, m)};
          const Numeric cosmm{cosm[m]};
          const Numeric sinmm{sinm[m]};
          const Numeric Pnm{P(n, m)};

          F.U += (gnm * cosmm + hnm * sinmm) * Pnm * (n + 1) * ratn;
          F.S += (gnm * cosmm + hnm * sinmm) * dP(n, m) * ratn;
          F.E += (gnm * sinmm + hnm * cosmm) * Pnm * m * ratn;
        }
      }
      F.S *= sin_theta;
      F.E /= sin_theta;

      // Adjust by theta (nb. can be nan when undefined, but set to zero instead)
      if (high_lat or low_lat) {
        F.S = F.total_horizontal();
        F.E = 0;
      }
    }
  }

  return out;
}
#pragma GCC diagnostic pop

constexpr Numeric next(Numeric p0, Numeric p1, Numeric ni, Numeric x) {
  return ((2.0 * ni + 1.0) * x * p1 - ni * p0) / (ni + 1.0);
}

Numeric legendre_sum(const ExhaustiveConstVectorView& s, const Numeric& x) {
  ARTS_USER_ERROR_IF(x < -1 or x > 1, "x=", x, " must be in [-1, 1]")

  const Index n = s.nelem();
  if (n == 0) return 0.0;

  Numeric p0 = 1.0;
  Numeric p1 = x;
  Numeric out = s[0] * p0;

  for (Index i = 1; i < n; i++) {
    out += s[i] * p1;
    std::swap(p0, p1);
    p1 = next(p1, p0, static_cast<Numeric>(i), x);
  }

  return out;
}

Numeric legendre(Index n, Numeric x) {
  if (n == 0) return 1.0;

  Numeric p0 = 1.0;
  Numeric p1 = x;
  for (Index i = 1; i < n; i++) {
    std::swap(p0, p1);
    p1 = next(p1, p0, static_cast<Numeric>(i), x);
  }

  return p1;
}

static constexpr Index max_factorial = 170;

constexpr Numeric factorial(Index i) {
  Numeric out = 1.0;
  for (Index j = 1; j <= i; j++) {
    out *= static_cast<Numeric>(j);
  }
  return out;
}

constexpr std::array<Numeric, max_factorial> factorial_table() {
  std::array<Numeric, max_factorial> table;
  for (Index i = 0; i < max_factorial; i++) {
    table[i] = factorial(i);
  }
  return table;
}

constexpr static auto factab = factorial_table();

//! Port of boost::math::double_factorial
Numeric double_factorial(Index i) {
  ARTS_USER_ERROR_IF(i < 0, "Must be postive. i=", i)
  if (i & 1) {
    // odd i:
    if (i < max_factorial) {
      Index n = (i - 1) / 2;
      return std::ceil(
          factab[i] / (std::ldexp(1.0, static_cast<int>(n)) * factab[n]) - 0.5);
    }
    Numeric result =
        std::tgamma(static_cast<Numeric>(i) / 2 + 1) * Constant::inv_sqrt_pi;
    if (std::ldexp(std::numeric_limits<Numeric>::max(),
                   -static_cast<int>(i + 1) / 2) > result)
      return std::ceil(result * std::ldexp(1.0, static_cast<int>(i + 1) / 2) -
                       0.5);
  } else {
    // even i:
    Index n = i / 2;
    Numeric result = factorial(n);
    if (std::ldexp(std::numeric_limits<Numeric>::max(), -(int)n) > result)
      return result * std::ldexp(1.0, (int)n);
  }

  ARTS_USER_ERROR("Overflow. i=", i);
}

//! Port of boost evaluate_rational
Numeric evaluate_rational(const std::array<Numeric, 13> num,
                          const std::array<std::uint32_t, 13>& denom,
                          const Numeric& z_) {
  static constexpr std::size_t count = 13;

  Numeric z(z_);
  Numeric s1, s2;
  if (z <= 1) {
    s1 = static_cast<Numeric>(num[count - 1]);
    s2 = static_cast<Numeric>(denom[count - 1]);
    for (int i = (int)count - 2; i >= 0; --i) {
      s1 *= z;
      s2 *= z;
      s1 += num[i];
      s2 += denom[i];
    }
  } else {
    z = 1 / z;
    s1 = static_cast<Numeric>(num[0]);
    s2 = static_cast<Numeric>(denom[0]);
    for (unsigned i = 1; i < count; ++i) {
      s1 *= z;
      s2 *= z;
      s1 += num[i];
      s2 += denom[i];
    }
  }
  return s1 / s2;
}

//! Port of boost lanczos_sum
Numeric lanczos_sum(const Numeric& z) {
  // LCOV_EXCL_START
  static const std::array num{
      static_cast<Numeric>(23531376880.41075968857200767445163675473L),
      static_cast<Numeric>(42919803642.64909876895789904700198885093L),
      static_cast<Numeric>(35711959237.35566804944018545154716670596L),
      static_cast<Numeric>(17921034426.03720969991975575445893111267L),
      static_cast<Numeric>(6039542586.35202800506429164430729792107L),
      static_cast<Numeric>(1439720407.311721673663223072794912393972L),
      static_cast<Numeric>(248874557.8620541565114603864132294232163L),
      static_cast<Numeric>(31426415.58540019438061423162831820536287L),
      static_cast<Numeric>(2876370.628935372441225409051620849613599L),
      static_cast<Numeric>(186056.2653952234950402949897160456992822L),
      static_cast<Numeric>(8071.672002365816210638002902272250613822L),
      static_cast<Numeric>(210.8242777515793458725097339207133627117L),
      static_cast<Numeric>(2.506628274631000270164908177133837338626L)};
  static const std::array denom{static_cast<std::uint32_t>(0u),
                                static_cast<std::uint32_t>(39916800u),
                                static_cast<std::uint32_t>(120543840u),
                                static_cast<std::uint32_t>(150917976u),
                                static_cast<std::uint32_t>(105258076u),
                                static_cast<std::uint32_t>(45995730u),
                                static_cast<std::uint32_t>(13339535u),
                                static_cast<std::uint32_t>(2637558u),
                                static_cast<std::uint32_t>(357423u),
                                static_cast<std::uint32_t>(32670u),
                                static_cast<std::uint32_t>(1925u),
                                static_cast<std::uint32_t>(66u),
                                static_cast<std::uint32_t>(1u)};
  // LCOV_EXCL_STOP
  return evaluate_rational(num, denom, z);
}

//! port of boost tgamma_delta_ratio_imp_lanczos
Numeric tgamma_delta_ratio_imp_lanczos(Numeric z, Numeric delta) {
  if (z < std::numeric_limits<Numeric>::epsilon()) {
    //
    // We get spurious numeric overflow unless we're very careful, this
    // can occur either inside Lanczos::lanczos_sum(z) or in the
    // final combination of terms, to avoid this, split the product up
    // into 2 (or 3) parts:
    //
    // G(z) / G(L) = 1 / (z * G(L)) ; z < eps, L = z + delta = delta
    //    z * G(L) = z * G(lim) * (G(L)/G(lim)) ; lim = largest factorial
    //
    if (max_factorial < delta) {
      Numeric ratio = tgamma_delta_ratio_imp_lanczos(
          delta, static_cast<Numeric>(max_factorial) - delta);
      ratio *= z;
      ratio *= factab[max_factorial - 1];
      return 1 / ratio;
    }
    return 1 / (z * std::tgamma(z + delta));
  }
  Numeric zgh = static_cast<Numeric>(z + 6.024680040776729583740234375) - 0.5;
  Numeric result;
  if (z + delta == z) {
    if (std::fabs(delta / zgh) < std::numeric_limits<Numeric>::epsilon()) {
      // We have:
      // result = exp((constants::half<T>() - z) * boost::math::log1p(delta / zgh, pol));
      // 0.5 - z == -z
      // log1p(delta / zgh) = delta / zgh = delta / z
      // multiplying we get -delta.
      result = std::exp(-delta);
    } else
      // from the pow formula below... but this may actually be wrong, we just can't really calculate it :(
      result = 1;
  } else {
    if (std::fabs(delta) < 10) {
      result = std::exp((0.5 - z) * std::log1p(delta / zgh));
    } else {
      result = pow(static_cast<Numeric>(zgh / (zgh + delta)),
                   static_cast<Numeric>(z - 0.5));
    }
    // Split the calculation up to avoid spurious overflow:
    result *= lanczos_sum(z) / lanczos_sum(z + delta);
  }
  result *=
      std::pow(static_cast<Numeric>(Constant::euler / (zgh + delta)), delta);
  return result;
}

//! port of boost tgamma_delta_ratio_imp
Numeric tgamma_delta_ratio(Numeric z, Numeric delta) {
  if ((z <= 0) || (z + delta <= 0)) {
    // This isn't very sophisticated, or accurate, but it does work:
    return std::tgamma(z) / std::tgamma(z + delta);
  }

  if (std::floor(delta) == delta) {
    if (std::floor(z) == z) {
      //
      // Both z and delta are integers, see if we can just use table lookup
      // of the factorials to get the result:
      //
      if ((z <= max_factorial) && (z + delta <= max_factorial)) {
        return factab[static_cast<Index>(std::trunc(z)) - 1] /
               factab[static_cast<Index>(std::trunc(z + delta)) - 1];
      }
    }
    if (std::fabs(delta) < 20) {
      //
      // delta is a small integer, we can use a finite product:
      //
      if (delta == 0) return 1;
      if (delta < 0) {
        z -= 1;
        Numeric result = z;
        while (0 != (delta += 1)) {
          z -= 1;
          result *= z;
        }
        return result;
      }
      Numeric result = 1 / z;
      while (0 != (delta -= 1)) {
        z += 1;
        result /= z;
      }
      return result;
    }
  }

  return tgamma_delta_ratio_imp_lanczos(z, delta);
}

//! port of boost tgamma_ratio_imp
Numeric tgamma_ratio(Numeric x, Numeric y) {
  ARTS_USER_ERROR_IF((x <= 0) || std::isinf(x), "Bad x=", x)
  ARTS_USER_ERROR_IF((y <= 0) || std::isinf(y), "Bad y=", y)

  if (x <= std::numeric_limits<Numeric>::min()) {
    // Special case for denorms...Ugh.
    Numeric shift = std::ldexp(1.0, std::numeric_limits<Numeric>::digits);
    return shift * tgamma_ratio(x * shift, y);
  }

  if ((x < max_factorial) && (y < max_factorial)) {
    // Rather than subtracting values, lets just call the gamma functions directly:
    return std::tgamma(x) / std::tgamma(y);
  }

  Numeric prefix = 1;
  if (x < 1) {
    if (y < 2 * max_factorial) {
      // We need to sidestep on x as well, otherwise we'll underflow
      // before we get to factor in the prefix term:
      prefix /= x;
      x += 1;
      while (y >= max_factorial) {
        y -= 1;
        prefix /= y;
      }
      return prefix * std::tgamma(x) / std::tgamma(y);
    }
    //
    // result is almost certainly going to underflow to zero, try logs just in case:
    //
    return std::exp(std::lgamma(x) - std::lgamma(y));
  }
  if (y < 1) {
    if (x < 2 * max_factorial) {
      // We need to sidestep on y as well, otherwise we'll overflow
      // before we get to factor in the prefix term:
      prefix *= y;
      y += 1;
      while (x >= max_factorial) {
        x -= 1;
        prefix *= x;
      }
      return prefix * std::tgamma(x) / std::tgamma(y);
    }
    //
    // Result will almost certainly overflow, try logs just in case:
    //
    return std::exp(std::lgamma(x) - std::lgamma(y));
  }
  //
  // Regular case, x and y both large and similar in magnitude:
  //
  return tgamma_delta_ratio(x, y - x);
}

//! port of boost legendre_next
Numeric next(Numeric l, Numeric m, Numeric x, Numeric Pl, Numeric Plm1) {
  return ((2 * l + 1) * x * Pl - (l + m) * Plm1) / (l + 1 - m);
}

// Port of Boost.Math legendre_p_imp
Numeric assoc_legendre(Index l, Index m, Numeric x) {
  if (x < -1 or x > 1) throw std::runtime_error("x must be in [-1, 1]");

  if (l < 0) return assoc_legendre(-l - 1, m, x);
  if ((l == 0) && (m == -1)) return std::sqrt((1 - x) / (1 + x));
  if ((l == 1) && (m == 0)) return x;
  if (-m == l)
    return std::pow((1 - x * x) / 4, static_cast<Numeric>(l) / 2) /
           std::tgamma(l + 1);
  if (m < 0) {
    return ((m & 1) ? -1 : 1) *
           tgamma_ratio(static_cast<Numeric>(l + m + 1),
                        static_cast<Numeric>(l + 1 - m)) *
           assoc_legendre(l, -m, x);
  }

  // Special cases:
  if (m > l) return 0;
  if (m == 0) return legendre(l, x);

  Numeric p0 = double_factorial(2 * m - 1) *
               std::pow(1.0 - x * x, static_cast<Numeric>(std::abs(m)) / 2);
  if (m & 1) p0 *= -1;
  if (m == l) return p0;

  Numeric p1 = x * static_cast<Numeric>(2 * m + 1) * p0;
  Index n = m + 1;
  while (n < l) {
    std::swap(p0, p1);
    p1 = next(static_cast<Numeric>(n), static_cast<Numeric>(m), x, p0, p1);
    ++n;
  }

  return p1;
}
}  // namespace Legendre
