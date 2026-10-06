#include "tro_legendre.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <legendre.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>

namespace scattering::tro_legendre {
namespace {
//! sqrt((2 l + 1) / 4 pi), Y_l0 = this P_l
Numeric ylm0_norm(Index l) { return std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * Constant::pi)); }

//! int P_l(x) dx = (P_{l+1}(x) - P_{l-1}(x)) / (2 l + 1), and x for l = 0, for l = 0..degree
void legendre_antiderivatives(VectorView q, Numeric x, Vector& p) {
  Legendre::legendre_polynomials(p, x);
  q[0] = x;
  for (Size l = 1; l < q.size(); l++) q[l] = (p[l + 1] - p[l - 1]) / static_cast<Numeric>(2 * l + 1);
}

//! The Gauss-Legendre rule of n nodes on [-1, 1]
const std::pair<Vector, Vector>& gauss_legendre(std::map<Index, std::pair<Vector, Vector>>& cache, Index n) {
  auto [it, inserted] = cache.try_emplace(n, Vector(n), Vector(n));
  if (inserted) Legendre::GaussLegendre(it->second.first, it->second.second);
  return it->second;
}
}  // namespace

Matrix project(const ConstVectorView& angles, const ConstMatrixView& values, Index degree) {
  const Index n = static_cast<Index>(angles.size()), k = values.ncols();
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)
  ARTS_USER_ERROR_IF(n < 1, "A Legendre projection needs at least one scattering angle")
  ARTS_USER_ERROR_IF(values.nrows() != n,
                     "The gridded data have {} scattering angles, but the scattering-angle grid has {}",
                     values.nrows(),
                     n)
  ARTS_USER_ERROR_IF(not(angles[0] >= 0.0 and angles[n - 1] <= 180.0),
                     "The scattering angles must be in [0, 180] deg, got [{}, {}]",
                     angles[0],
                     angles[n - 1])
  for (Index i = 0; i + 1 < n; i++)
    ARTS_USER_ERROR_IF(not(angles[i] < angles[i + 1]),
                       "The scattering angles must ascend strictly, but angle {} is {} deg and angle {} is {} deg",
                       i,
                       angles[i],
                       i + 1,
                       angles[i + 1])

  // b[l, j] = int_{-1}^{1} F_j(x) P_l(x) dx
  Matrix b(degree + 1, k, 0.0);
  Vector p(degree + 2), qa(degree + 1), qb(degree + 1);

  // The constant ends, Theta in [0, angles[0]] and [angles[n - 1], 180], i.e. x in [cos(angles[0]), 1] and
  // [-1, cos(angles[n - 1])], integrated in closed form
  const auto constant = [&](Numeric x_low, Numeric x_high, Index node) {
    if (not(x_low < x_high)) return;
    legendre_antiderivatives(qa, x_low, p);
    legendre_antiderivatives(qb, x_high, p);
    for (Index l = 0; l <= degree; l++)
      for (Index j = 0; j < k; j++) b[l, j] += (qb[l] - qa[l]) * values[node, j];
  };
  constant(std::cos(Conversion::deg2rad(angles[0])), 1.0, 0);
  constant(-1.0, std::cos(Conversion::deg2rad(angles[n - 1])), n - 1);

  // The linear segments, int F(Theta) P_l(cos(Theta)) sin(Theta) dTheta.  P_l(cos(Theta)) sin(Theta) is a
  // trigonometric polynomial of degree l + 1; a Gauss rule of (degree + 2) h + 10 nodes on a segment of h
  // radians integrates it, times the linear F, to rounding
  std::map<Index, std::pair<Vector, Vector>> rules;
  for (Index i = 0; i + 1 < n; i++) {
    const Numeric t0 = Conversion::deg2rad(angles[i]), t1 = Conversion::deg2rad(angles[i + 1]);
    const Numeric h = t1 - t0, mid = 0.5 * (t0 + t1);
    const auto& [xi, wi] =
        gauss_legendre(rules, static_cast<Index>(std::ceil(static_cast<Numeric>(degree + 2) * h)) + 10);
    for (Size q = 0; q < xi.size(); q++) {
      const Numeric t = mid + 0.5 * h * xi[q], w = 0.5 * h * wi[q] * std::sin(t), r = (t - t0) / h;
      Legendre::legendre_polynomials(p, std::clamp(std::cos(t), -1.0, 1.0));
      for (Index j = 0; j < k; j++) {
        const Numeric f = w * ((1.0 - r) * values[i, j] + r * values[i + 1, j]);
        for (Index l = 0; l <= degree; l++) b[l, j] += f * p[l];
      }
    }
  }

  // a_l = 2 pi sqrt((2 l + 1) / 4 pi) b_l
  for (Index l = 0; l <= degree; l++) b[l, joker] *= 2.0 * Constant::pi * ylm0_norm(l);
  return b;
}

Matrix evaluate(const ConstMatrixView& coefficients, const ConstVectorView& angles) {
  const Index L = coefficients.nrows() - 1, k = coefficients.ncols();
  ARTS_USER_ERROR_IF(L < 0, "A Legendre series needs at least one coefficient")
  Matrix out(angles.size(), k, 0.0);
  Vector p(L + 1);
  for (Size i = 0; i < angles.size(); i++) {
    ARTS_USER_ERROR_IF(not(angles[i] >= 0.0 and angles[i] <= 180.0),
                       "The scattering angles must be in [0, 180] deg, got {} deg",
                       angles[i])
    Legendre::legendre_polynomials(p, std::clamp(std::cos(Conversion::deg2rad(angles[i])), -1.0, 1.0));
    for (Index l = 0; l <= L; l++) {
      const Numeric y = ylm0_norm(l) * p[l];
      for (Index j = 0; j < k; j++) out[i, j] += coefficients[l, j] * y;
    }
  }
  return out;
}

report assess(const ConstMatrixView& coefficients, const ConstVectorView& angles, const ConstMatrixView& values) {
  ARTS_USER_ERROR_IF(coefficients.ncols() != 6 or values.ncols() != 6,
                     "A TRO report needs the six elements [F11, F12, F22, F33, F34, F44]")
  const Index L = coefficients.nrows() - 1;
  constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

  report out;
  Numeric scale = 0.0;
  for (Index i = 0; i < values.nrows(); i++) scale = std::max(scale, std::abs(values[i, 0]));
  const Numeric inv_scale = scale > 0.0 ? 1.0 / scale : nan;

  const Matrix series = evaluate(coefficients, angles);
  for (Index i = 0; i < values.nrows(); i++)
    for (Index j = 0; j < 6; j++)
      out.reconstruction_error[j] =
          std::max(out.reconstruction_error[j], std::abs(series[i, j] - values[i, j]) * inv_scale);

  const Numeric a0 = coefficients[0, 0];
  for (Index j = 0; j < 6; j++) out.tail[j] = a0 != 0.0 ? std::abs(coefficients[L, j] / a0) : nan;
  out.asymmetry = (L >= 1 and a0 != 0.0) ? coefficients[1, 0] / (std::numbers::sqrt3 * a0) : nan;

  // The series between and beyond the nodes, on a uniform grid fine enough to resolve degree L
  const Index nfine = 8 * (L + 1) + 1;
  Vector      fine(nfine + angles.size());
  for (Index i = 0; i < nfine; i++) fine[i] = 180.0 * static_cast<Numeric>(i) / static_cast<Numeric>(nfine - 1);
  for (Size i = 0; i < angles.size(); i++) fine[nfine + i] = angles[i];
  const Matrix fine_series = evaluate(coefficients, fine);
  Numeric      lowest      = std::numeric_limits<Numeric>::infinity();
  for (Size i = 0; i < fine.size(); i++) lowest = std::min(lowest, fine_series[i, 0]);
  out.min_f11 = lowest * inv_scale;
  return out;
}
}  // namespace scattering::tro_legendre

namespace scattering {
LegendreReport::LegendreReport(Index n_temps, Index n_freqs)
    : reconstruction_error(n_temps, n_freqs, 6, 0.0),
      tail(n_temps, n_freqs, 6, 0.0),
      min_f11(n_temps, n_freqs, 0.0),
      asymmetry(n_temps, n_freqs, 0.0),
      normalisation_error(n_temps, n_freqs, std::numeric_limits<Numeric>::quiet_NaN()) {}

void LegendreReport::set(Index i_t, Index i_f, const tro_legendre::report& r) {
  reconstruction_error[i_t, i_f, joker] = r.reconstruction_error;
  tail[i_t, i_f, joker]                 = r.tail;
  min_f11[i_t, i_f]                     = r.min_f11;
  asymmetry[i_t, i_f]                   = r.asymmetry;
}
}  // namespace scattering
