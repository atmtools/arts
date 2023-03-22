#include "arts_conversions.h"
#include "legendre.h"

namespace Legendre {
//! Converts latitude to co-latitude with a small distance from the poles and flags if this was activated
struct ColatitudeConversion {
  static constexpr Numeric dlat_lim = 1e-4;

  const bool high;
  const bool low;
  const Numeric theta;

  constexpr ColatitudeConversion(const Numeric latitude) noexcept :
  high(latitude > (90.0 - dlat_lim)), low(latitude < (dlat_lim - 90.0)),
  theta(high ? Conversion::deg2rad(dlat_lim) : (low ? Conversion::deg2rad(180.0 - dlat_lim) : Conversion::deg2rad(90.0 - latitude))) {}
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
std::pair<Matrix, Matrix> schmidt(const Numeric theta, const Index nmax) ARTS_NOEXCEPT {
  ARTS_ASSERT(theta > 0 and theta < Constant::pi)
  ARTS_ASSERT(nmax > 0)

  using std::sqrt;
  using Math::pow2;

  Index N = 1 + nmax;

  const Numeric x = std::cos(theta);
  Matrix P(N, N, 0);
  Matrix dP(N, N, 0);

  // Start values
  P(0, 0) = 1.0;
  P(1, 0) = x;
  P(1, 1) = - sqrt(1.0 - pow2(x));
  dP(0, 0) = 0.0;
  dP(1, 0) = 1.0;

  for (Index n=2; n<N; ++n) {
    P(n, 0) = ((2 * (n-1) + 1) * x * P(n-1, 0) - (n - 1) * P(n-2, 0)) / n;
  }

  dP(nmax, 0) = nmax / (pow2(x) - 1) * (x * P(nmax, 0) - P(nmax-1, 0));
  for (Index n=2; n<nmax; ++n) {
    dP(n, 0) = (n + 1) / (pow2(x) - 1) * (P(n+1, 0) - x * P(n, 0));
  }

  Numeric Cm = sqrt(2.0);
  for (Index m=1; m<N; ++m) {
    Cm /= sqrt(Numeric(2*m * (2*m - 1)));
    P(m, m) = std::pow(1 - pow2(x), 0.5*m) * Cm;

    for (Index i=1; i<m; ++i) {
      P(m, m) *= 2*i + 1;
    }

    dP(m, m) = -P(m, m) * m * x / sqrt(1 - pow2(x));

    if (nmax > m) {
      Numeric twoago = 0.0;
      for (Index n=m+1; n<N; ++n) {
        P(n, m) = (x * (2*n - 1) * P(n-1, m) - sqrt(Numeric((n+m-1) * (n-m-1))) * twoago) / sqrt(Numeric(pow2(n) - pow2(m)));
        twoago = P(n-1, m);
      }
    }
  }

  for (Index n=2; n<N; ++n) {
    for (Index m=1; m<n; ++m) {
      dP(n, m) = sqrt(Numeric((n-m) * (n+m+1))) * P(n, m+1) - P(n, m) * m * x / sqrt(1 - pow2(x));
    }
  }

  return {P, dP};
}


SphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0,
                                 const Numeric r, const Numeric lat, const Numeric lon) ARTS_NOEXCEPT {
  ARTS_ASSERT(lat >= -90 and lat <= 90, "Latitude is ", lat, " should be in [-90, 90]")

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
    for (Index m=0; m<N; ++m) {
      cosm[m] = 1.0;
      sinm[m] = 0.0;
    }
  } else {
    for (Index m=0; m<N; ++m) {
      cosm[m] = Conversion::cosd(m * clon);
      sinm[m] = Conversion::sind(m * clon);
    }
  }

  // Compute the field
  SphericalField F;  // Auto init to zeroes
  for (Index n=1; n<N; ++n) {
    const Numeric ratn = std::pow(r0 / r, n + 2);
    for (Index m=0; m<n+1; ++m) {
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


MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Vector& rv, const Numeric lat, const Vector& lonv) ARTS_NOEXCEPT {
  ARTS_ASSERT(lat >= -90 and lat <= 90, "Latitude is ", lat, " should be in [-90, 90]")

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
  for (Index ir=0; ir<nr; ir++) {
    ARTS_ASSERT(rv[ir] > 0, "Radius are [", rv, "] must be in (0, inf]")
    for (Index n=1; n<N; ++n) {
      ratnm(ir, n) = std::pow(r0 / rv[ir], n + 2);
    }
  }

  MatrixOfSphericalField out(nlon, nr);  // Auto init to zeroes
  for (Index ilon=0; ilon<nlon; ilon++) {

    // Pre-compute the cosine/sine values
    std::vector<Numeric> cosm(N);
    std::vector<Numeric> sinm(N);
    const Numeric clon = longitude_clamp(lonv[ilon]);
    if (low_lat or high_lat) {
      for (Index m=0; m<N; ++m) {
        cosm[m] = 1.0;
        sinm[m] = 0.0;
      }
    } else {
      for (Index m=0; m<N; ++m) {
        cosm[m] = Conversion::cosd(m * clon);
        sinm[m] = Conversion::sind(m * clon);
      }
    }

    // Compute the field
    for (Index ir=0; ir<nr; ir++) {
      SphericalField& F = out(ilon, ir);  // Auto init'd to zeroes
      for (Index n=1; n<N; ++n) {
        const Numeric ratn = ratnm(ir, n);
        for (Index m=0; m<n+1; ++m) {
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
} // namespace Legendre
