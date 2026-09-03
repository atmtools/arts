#include "rtepack_transmission.h"

#include <array_algo.h>
#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_omp.h>
#include <debug.h>

#include <Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_spectral_matrix.h"

namespace rtepack {
static constexpr Numeric too_small = 1e-4;
// 8 * eps^(1/4) for IEEE binary64.  Below this relative separation the
// eigenvalue divided differences lose enough digits to warrant the rare
// scaling-and-squaring path.
static constexpr Numeric spectral_cluster_relative = 1.0 / 1024.0;

namespace {
constexpr propmat commutator(const propmat &lhs, const propmat &rhs) {
  return {0.0,
          rhs.C() * lhs.U() - lhs.C() * rhs.U() + rhs.D() * lhs.V() - lhs.D() * rhs.V(),
          lhs.B() * rhs.U() - rhs.B() * lhs.U() + rhs.D() * lhs.W() - lhs.D() * rhs.W(),
          lhs.B() * rhs.V() - rhs.B() * lhs.V() - rhs.C() * lhs.W() + lhs.C() * rhs.W(),
          lhs.B() * rhs.C() - rhs.B() * lhs.C() + rhs.V() * lhs.W() - lhs.V() * rhs.W(),
          lhs.B() * rhs.D() - rhs.B() * lhs.D() - rhs.U() * lhs.W() + lhs.U() * rhs.W(),
          lhs.C() * rhs.D() - rhs.C() * lhs.D() + rhs.U() * lhs.V() - lhs.U() * rhs.V()};
}

constexpr propmat magnus_exponent_deriv(
    const propmat &k1, const propmat &k2, const propmat &dk, const Numeric r, const Numeric dr, const bool k1_deriv) {
  const propmat dk1  = k1_deriv ? dk : propmat{};
  const propmat dk2  = k1_deriv ? propmat{} : dk;
  const propmat rk1  = r * k1;
  const propmat rk2  = r * k2;
  const propmat drk1 = dr * k1 + r * dk1;
  const propmat drk2 = dr * k2 + r * dk2;

  return -dr * avg(k1, k2) - r * avg(dk1, dk2) + (1.0 / 12.0) * (commutator(drk2, rk1) + commutator(rk2, drk1));
}

Numeric max_abs(const muelmat &m) {
  Numeric out{};
  for (Size i = 0; i < 4; ++i) {
    for (Size j = 0; j < 4; ++j) out = std::max(out, std::abs(m[i, j]));
  }
  return out;
}

// Return s such that ||2^-s m||_inf <= 1/2 without ever forming an
// overflowing norm or the (possibly overflowing) value 2^s.  Dividing by
// the largest entry before accumulating the row norm also keeps this defined
// for every finite binary64 matrix.
int taylor_scaling_steps(const muelmat &m) {
  const Numeric largest = max_abs(m);
  if (largest == 0.0 or not std::isfinite(largest)) return 0;

  Numeric scaled_row_norm{};
  for (Size i = 0; i < 4; ++i) {
    Numeric sum{};
    for (Size j = 0; j < 4; ++j) sum += std::abs(m[i, j]) / largest;
    scaled_row_norm = std::max(scaled_row_norm, sum);
  }

  const Numeric log2_norm = std::log2(largest) + std::log2(scaled_row_norm);
  return log2_norm > -1.0 ? std::max(0, static_cast<int>(std::ceil(log2_norm + 1.0))) : 0;
}

Numeric phi1_scalar(const Numeric z) { return z == 0.0 ? 1.0 : std::expm1(z) / z; }

Numeric phi1m1_scalar_from_expm1(const Numeric z, const Numeric expm1_z) {
  if (std::abs(z) >= 0.5) return expm1_z / z - 1.0;

  Numeric term = 0.5 * z;
  Numeric sum  = term;
  for (Size n = 2; n <= 18; ++n) {
    term *= z / static_cast<Numeric>(n + 1);
    sum  += term;
    if (std::abs(term) <= std::numeric_limits<Numeric>::epsilon() * std::abs(sum)) break;
  }
  return sum;
}

template <Size derivative_order> Numeric phi1_scalar_derivative(const Numeric z) {
  static_assert(derivative_order >= 1 and derivative_order <= 4);

  // phi_1^(k)(z) = sum_n z^n / (n! (n+k+1)).  The closed forms below
  // subtract nearly equal numbers near zero, increasingly severely with k.
  if (std::abs(z) < 0.5) {
    Numeric sum  = 1.0 / static_cast<Numeric>(derivative_order + 1);
    Numeric term = 1.0;
    for (Size n = 1; n <= 18; ++n) {
      term              *= z / static_cast<Numeric>(n);
      const Numeric add  = term / static_cast<Numeric>(n + derivative_order + 1);
      sum               += add;
      if (std::abs(add) <= std::numeric_limits<Numeric>::epsilon() * std::abs(sum)) break;
    }
    return sum;
  }

  const Numeric ez = std::exp(z);
  if constexpr (derivative_order == 1) {
    return (ez * (z - 1.0) + 1.0) / Math::pow2(z);
  } else if constexpr (derivative_order == 2) {
    return (ez * (Math::pow2(z) - 2.0 * z + 2.0) - 2.0) / Math::pow3(z);
  } else if constexpr (derivative_order == 3) {
    const Numeric z2 = Math::pow2(z);
    const Numeric z3 = z2 * z;
    return (ez * (z3 - 3.0 * z2 + 6.0 * z - 6.0) + 6.0) / Math::pow4(z);
  } else {
    const Numeric z2 = Math::pow2(z);
    const Numeric z3 = z2 * z;
    const Numeric z4 = z2 * z2;
    return (ez * (z4 - 4.0 * z3 + 12.0 * z2 - 24.0 * z + 24.0) - 24.0) / (z4 * z);
  }
}

Numeric phi1_scalar_derivative_from_exp(const Numeric z, const Numeric exp_z) {
  return std::abs(z) < 0.5 ? phi1_scalar_derivative<1>(z) : (exp_z * (z - 1.0) + 1.0) / Math::pow2(z);
}

Numeric polarized_row_norm(const propmat &g) {
  return std::max({std::abs(g.B()) + std::abs(g.C()) + std::abs(g.D()),
                   std::abs(g.B()) + std::abs(g.U()) + std::abs(g.V()),
                   std::abs(g.C()) + std::abs(g.U()) + std::abs(g.W()),
                   std::abs(g.D()) + std::abs(g.V()) + std::abs(g.W())});
}

muelmat phi1_cayley_hamilton_series(const propmat &g, Vector4 *diag_m1 = nullptr) {
  const Numeric a  = g.A();
  const Numeric b  = g.B();
  const Numeric c  = g.C();
  const Numeric d  = g.D();
  const Numeric u  = g.U();
  const Numeric v  = g.V();
  const Numeric w  = g.W();
  const Numeric B  = u * u + v * v + w * w - b * b - c * c - d * d;
  const Numeric C  = -Math::pow2(d * u - c * v + b * w);
  const Numeric n  = polarized_row_norm(g);
  const Numeric n2 = n * n;
  const Numeric n3 = n2 * n;

  Vector4 term{1.0, 0.0, 0.0, 0.0};
  Vector4 sum = term;
  Numeric sum0_m1{};

  for (Size k = 0; k < 80; ++k) {
    const Numeric inv_denominator = 1.0 / static_cast<Numeric>(k + 2);
    const auto    old             = term;
    term                          = {inv_denominator * (a * old[0] - C * old[3]),
                                     inv_denominator * (old[0] + a * old[1]),
                                     inv_denominator * (old[1] + a * old[2] - B * old[3]),
                                     inv_denominator * (old[2] + a * old[3])};
    for (Size j = 0; j < 4; ++j) sum[j] += term[j];
    sum0_m1 += term[0];

    const Numeric term_bound =
        std::abs(term[0]) + n * std::abs(term[1]) + n2 * std::abs(term[2]) + n3 * std::abs(term[3]);
    const Numeric sum_bound = std::abs(sum[0]) + n * std::abs(sum[1]) + n2 * std::abs(sum[2]) + n3 * std::abs(sum[3]);
    if (static_cast<Numeric>(k + 2) > std::abs(a) + n and
        term_bound <= std::numeric_limits<Numeric>::epsilon() * (1.0 + sum_bound))
      break;
  }

  const muelmat N{0, b, c, d, b, 0, u, v, c, -u, 0, w, d, -v, -w, 0};
  const muelmat N2 = N * N;
  const muelmat N3 = N * N2;
  if (diag_m1 != nullptr) {
    *diag_m1 = {sum0_m1 + sum[2] * N2[0, 0],
                sum0_m1 + sum[2] * N2[1, 1],
                sum0_m1 + sum[2] * N2[2, 2],
                sum0_m1 + sum[2] * N2[3, 3]};
  }
  return sum[0] * muelmat::id() + sum[1] * N + sum[2] * N2 + sum[3] * N3;
}

muelmat phi1_frechet_cayley_hamilton_series(const propmat &g, const propmat &dg) {
  const auto &[a, b, c, d, u, v, w]        = g.data;
  const auto &[da, db, dc, dd, du, dv, dw] = dg.data;
  const Numeric B                          = u * u + v * v + w * w - b * b - c * c - d * d;
  const Numeric C                          = -Math::pow2(d * u - c * v + b * w);
  const Numeric dB                         = 2.0 * (u * du + v * dv + w * dw - b * db - c * dc - d * dd);
  const Numeric dC = -2.0 * (d * u - c * v + b * w) * (dd * u + d * du - dc * v - c * dv + db * w + b * dw);
  const Numeric n  = polarized_row_norm(g);
  const Numeric e  = polarized_row_norm(propmat{0.0, db, dc, dd, du, dv, dw});
  const Numeric n2 = n * n;
  const Numeric n3 = n2 * n;

  Vector4 term{1.0, 0.0, 0.0, 0.0};
  Vector4 dterm{};
  Vector4 sum = term;
  Vector4 dsum{};

  for (Size k = 0; k < 80; ++k) {
    const Numeric inv_denominator = 1.0 / static_cast<Numeric>(k + 2);
    const auto    old             = term;
    const auto    dold            = dterm;
    term                          = {inv_denominator * (a * old[0] - C * old[3]),
                                     inv_denominator * (old[0] + a * old[1]),
                                     inv_denominator * (old[1] + a * old[2] - B * old[3]),
                                     inv_denominator * (old[2] + a * old[3])};
    dterm = {inv_denominator * (da * old[0] + a * dold[0] - dC * old[3] - C * dold[3]),
             inv_denominator * (dold[0] + da * old[1] + a * dold[1]),
             inv_denominator * (dold[1] + da * old[2] + a * dold[2] - dB * old[3] - B * dold[3]),
             inv_denominator * (dold[2] + da * old[3] + a * dold[3])};
    for (Size j = 0; j < 4; ++j) {
      sum[j]  += term[j];
      dsum[j] += dterm[j];
    }

    const Numeric term_bound =
        std::abs(term[0]) + n * std::abs(term[1]) + n2 * std::abs(term[2]) + n3 * std::abs(term[3]);
    const Numeric dterm_bound = std::abs(dterm[0]) + n * std::abs(dterm[1]) + n2 * std::abs(dterm[2]) +
                                n3 * std::abs(dterm[3]) +
                                e * (std::abs(term[1]) + 2.0 * n * std::abs(term[2]) + 3.0 * n2 * std::abs(term[3]));
    const Numeric sum_bound   = std::abs(sum[0]) + n * std::abs(sum[1]) + n2 * std::abs(sum[2]) + n3 * std::abs(sum[3]);
    const Numeric dsum_bound  = std::abs(dsum[0]) + n * std::abs(dsum[1]) + n2 * std::abs(dsum[2]) +
                                n3 * std::abs(dsum[3]) +
                                e * (std::abs(sum[1]) + 2.0 * n * std::abs(sum[2]) + 3.0 * n2 * std::abs(sum[3]));
    if (static_cast<Numeric>(k + 2) > std::abs(a) + n and
        term_bound <= std::numeric_limits<Numeric>::epsilon() * (1.0 + sum_bound) and
        dterm_bound <= std::numeric_limits<Numeric>::epsilon() * (1.0 + dsum_bound))
      break;
  }

  const muelmat N{0, b, c, d, b, 0, u, v, c, -u, 0, w, d, -v, -w, 0};
  const muelmat dN{0, db, dc, dd, db, 0, du, dv, dc, -du, 0, dw, dd, -dv, -dw, 0};
  const muelmat N2  = N * N;
  const muelmat dN2 = dN * N + N * dN;
  const muelmat N3  = N * N2;
  const muelmat dN3 = dN * N2 + N * dN2;
  return dsum[0] * muelmat::id() + dsum[1] * N + sum[1] * dN + dsum[2] * N2 + sum[2] * dN2 + dsum[3] * N3 +
         sum[3] * dN3;
}

// Scaling-and-squaring Taylor evaluation of the Frechet derivative of exp.
// This is used only for degenerate polarized generators, where the closed
// Cayley-Hamilton derivative has removable singularities in x, y, and S.
muelmat exp_frechet_degenerate(const propmat &g, const propmat &dg) {
  muelmat A = to_muelmat(g);
  muelmat E = to_muelmat(dg);

  const int     squarings  = taylor_scaling_steps(A);
  const Numeric scale      = std::scalbn(1.0, -squarings);
  A                       *= scale;
  E                       *= scale;

  muelmat value = muelmat::id();
  muelmat deriv = muelmat::constant(0.0);
  muelmat term  = muelmat::id();
  muelmat dterm = muelmat::constant(0.0);

  for (Size n = 1; n <= 32; ++n) {
    const muelmat previous  = term;
    dterm                   = (dterm * A + previous * E) / static_cast<Numeric>(n);
    term                    = (previous * A) / static_cast<Numeric>(n);
    value                  += term;
    deriv                  += dterm;

    if (max_abs(term) + max_abs(dterm) <=
        std::numeric_limits<Numeric>::epsilon() * (1.0 + max_abs(value) + max_abs(deriv)))
      break;
  }

  for (int i = 0; i < squarings; ++i) {
    deriv = value * deriv + deriv * value;
    value = value * value;
  }

  return deriv;
}

// Scaling-and-squaring Taylor evaluation of exp(G).  This is reserved for
// generators for which the faster spectral coefficients are ill-conditioned.
muelmat exp_degenerate(const propmat &g) {
  muelmat A = to_muelmat(g);

  const int     squarings  = taylor_scaling_steps(A);
  const Numeric scale      = std::scalbn(1.0, -squarings);
  A                       *= scale;

  muelmat value = muelmat::id();
  muelmat term  = muelmat::id();
  for (Size n = 1; n <= 32; ++n) {
    term   = (term * A) / static_cast<Numeric>(n);
    value += term;

    if (max_abs(term) <= std::numeric_limits<Numeric>::epsilon() * (1.0 + max_abs(value))) break;
  }

  for (int i = 0; i < squarings; ++i) value = value * value;
  return value;
}

// Scaling-and-squaring Taylor evaluation of exp(G)-I.  Squaring the remainder
// as E(2G) = 2E(G) + E(G)^2 avoids subtracting the identity at every scale.
muelmat expm1_degenerate(const propmat &g) {
  muelmat A = to_muelmat(g);

  const int     squarings  = taylor_scaling_steps(A);
  const Numeric scale      = std::scalbn(1.0, -squarings);
  A                       *= scale;

  muelmat value = A;
  muelmat term  = A;
  for (Size n = 2; n <= 32; ++n) {
    term   = (term * A) / static_cast<Numeric>(n);
    value += term;

    if (max_abs(term) <= std::numeric_limits<Numeric>::epsilon() * (1.0 + max_abs(value))) break;
  }

  for (int i = 0; i < squarings; ++i) value = 2.0 * value + value * value;
  return value;
}

// Scaling-and-squaring Taylor evaluation of phi_1(G).  Unlike the
// Cayley-Hamilton divided differences, this remains well conditioned when the
// polarized eigenvalues coalesce.
muelmat phi1_degenerate(const propmat &g) {
  muelmat A = to_muelmat(g);

  const int     squarings  = taylor_scaling_steps(A);
  const Numeric scale      = std::scalbn(1.0, -squarings);
  A                       *= scale;

  muelmat value = muelmat::id();
  muelmat term  = muelmat::id();
  for (Size n = 1; n <= 32; ++n) {
    term   = (term * A) / static_cast<Numeric>(n + 1);
    value += term;

    if (max_abs(term) <= std::numeric_limits<Numeric>::epsilon() * (1.0 + max_abs(value))) break;
  }

  // phi_1(2A) = (I + exp(A)) phi_1(A) / 2.
  muelmat exp_value = muelmat::id() + A * value;
  for (int i = 0; i < squarings; ++i) {
    value     = 0.5 * ((muelmat::id() + exp_value) * value);
    exp_value = exp_value * exp_value;
  }

  return value;
}

// Scaling-and-squaring Taylor evaluation of D phi_1(G)[dG].
muelmat phi1_frechet_degenerate(const propmat &g, const propmat &dg) {
  muelmat A = to_muelmat(g);
  muelmat E = to_muelmat(dg);

  const int     squarings  = taylor_scaling_steps(A);
  const Numeric scale      = std::scalbn(1.0, -squarings);
  A                       *= scale;
  E                       *= scale;

  muelmat value = muelmat::id();
  muelmat deriv = muelmat::constant(0.0);
  muelmat term  = muelmat::id();
  muelmat dterm = muelmat::constant(0.0);

  for (Size n = 1; n <= 32; ++n) {
    const muelmat previous  = term;
    dterm                   = (dterm * A + previous * E) / static_cast<Numeric>(n + 1);
    term                    = (previous * A) / static_cast<Numeric>(n + 1);
    value                  += term;
    deriv                  += dterm;

    if (max_abs(term) + max_abs(dterm) <=
        std::numeric_limits<Numeric>::epsilon() * (1.0 + max_abs(value) + max_abs(deriv)))
      break;
  }

  muelmat exp_value = muelmat::id() + A * value;
  muelmat exp_deriv = E * value + A * deriv;

  for (int i = 0; i < squarings; ++i) {
    const muelmat old_value     = value;
    const muelmat old_deriv     = deriv;
    const muelmat old_exp_value = exp_value;
    const muelmat old_exp_deriv = exp_deriv;

    value     = 0.5 * ((muelmat::id() + old_exp_value) * old_value);
    deriv     = 0.5 * (old_exp_deriv * old_value + (muelmat::id() + old_exp_value) * old_deriv);
    exp_value = old_exp_value * old_exp_value;
    exp_deriv = old_exp_deriv * old_exp_value + old_exp_value * old_exp_deriv;
  }

  return deriv;
}

struct scalar_linprop_result {
  Numeric value;
  Numeric value_m1;
  Numeric dtau;
  Numeric ddelta;
};

// Exact scalar linear-propagation source operator and its two partial
// derivatives.  With ks at the incoming/background end and ke at the
// outgoing/observer end,
//
//   L(tau, delta) = integral_0^1 exp(-tau*q-delta*q*(1-q)) dq,
//   tau = r*(ks+ke)/2, delta = r*(ke-ks)/2.
//
// The three branches avoid cancellation respectively for optically thick
// layers, nearly constant propagation, and well-separated endpoints.
scalar_linprop_result scalar_linprop(const Numeric start_optical,
                                     const Numeric end_optical,
                                     const Numeric minus_tau,
                                     const Numeric exp_minus_tau,
                                     const Numeric exp_minus_tau_m1) {
  const Numeric eps   = std::numeric_limits<Numeric>::epsilon();
  const Numeric tau   = -minus_tau;
  const Numeric delta = std::midpoint(end_optical, -start_optical);
  const Numeric beta  = end_optical;

  // Watson expansion about the observer-side endpoint.  Its finite-interval
  // remainder is below round-off once tau >= 40, and the ratio condition
  // keeps the asymptotic terms safely decreasing.
  if (tau >= 40.0 and beta > 0.0 and beta * beta >= 256.0 * std::abs(delta)) {
    scalar_linprop_result out{};
    Numeric               term = 1.0 / beta;

    for (Size n = 0; n < 64; ++n) {
      const Numeric order  = static_cast<Numeric>(2 * n + 1);
      out.value           += term;
      out.dtau            -= order * term / beta;
      out.ddelta          += order * term / beta * (-1.0 + static_cast<Numeric>(2 * n + 2) / beta);

      const Numeric next = term * (2.0 * delta * order) / (beta * beta);
      if (std::abs(next) > std::abs(term)) break;
      if (std::abs(next) <= eps * std::max(std::abs(out.value), 1.0 / beta)) break;
      term = next;
    }
    out.value_m1 = out.value - 1.0;
    return out;
  }

  // Expand only in the endpoint difference.  The exponential moments are
  // generated in their stable recurrence direction, so this branch is
  // continuous at delta=0 without sacrificing thick-layer accuracy.
  if (std::abs(delta) <= 4.0) {
    constexpr Size max_order  = 20;
    constexpr Size max_moment = 2 * (max_order + 1);

    const Numeric x = 0.25 * std::abs(delta);
    Size          order{};
    Numeric       omitted = x;
    // exp(x) <= e here.  A constant conservative bound avoids a libm call
    // in the common zero/small-gradient path.
    while (order < max_order and omitted > eps / 32.0) {
      ++order;
      omitted *= x / static_cast<Numeric>(order + 1);
    }

    const Size                          last_moment = 2 * (order + 1);
    std::array<Numeric, max_moment + 1> moment;

    if (tau == 0.0) {
      for (Size m = 0; m <= last_moment; ++m) moment[m] = 1.0 / static_cast<Numeric>(m + 1);
    } else if (tau > 0.0) {
      moment[0] = -exp_minus_tau_m1 / tau;

      const Size upward_last = tau >= static_cast<Numeric>(last_moment) ? last_moment : static_cast<Size>(tau);
      for (Size m = 1; m <= upward_last; ++m)
        moment[m] = (static_cast<Numeric>(m) * moment[m - 1] - exp_minus_tau) / tau;

      if (upward_last < last_moment) {
        // Positive endpoint series for I_M.  This is an independent stable
        // anchor for the downward recurrence when M is larger than tau.
        Numeric series = 1.0 / static_cast<Numeric>(last_moment + 1);
        Numeric term   = series;
        for (Size j = 1; j < 256; ++j) {
          term   *= tau / static_cast<Numeric>(last_moment + j + 1);
          series += term;
          if (term <= eps * series) break;
        }
        moment[last_moment] = exp_minus_tau * series;

        for (Size m = last_moment; m > upward_last + 1; --m)
          moment[m - 1] = (tau * moment[m] + exp_minus_tau) / static_cast<Numeric>(m);
        if (upward_last == 0)
          moment[0] = -exp_minus_tau_m1 / tau;
        else
          moment[upward_last] = (tau * moment[upward_last + 1] + exp_minus_tau) / static_cast<Numeric>(upward_last + 1);
      }
    } else {
      // Negative scalar propagation is nonphysical but remains useful to the
      // public low-level API.  A direct entire series is stable in the range
      // where the result itself is representable and this branch is selected.
      for (Size m = 0; m <= last_moment; ++m) {
        Numeric sum  = 1.0 / static_cast<Numeric>(m + 1);
        Numeric term = 1.0;
        for (Size j = 1; j < 256; ++j) {
          term              *= -tau / static_cast<Numeric>(j);
          const Numeric add  = term / static_cast<Numeric>(m + j + 1);
          sum               += add;
          if (std::abs(add) <= eps * std::max(1.0, std::abs(sum))) break;
        }
        moment[m] = sum;
      }
    }

    // Keep this path explicitly binary64 on every platform.  In particular,
    // long double is binary64 on arm64 macOS but x87 extended precision on
    // x86 Linux, which otherwise makes this hot retrieval path platform
    // dependent.  Compensated sums retain the cancellation protection.
    std::array<Numeric, max_order + 2> A;
    std::array<Numeric, max_order + 1> Bm;
    for (Size n = 0; n <= order + 1; ++n) {
      Numeric sum_a{}, correction_a{}, sum_b{}, correction_b{}, choose{1.0};
      for (Size j = 0; j <= n; ++j) {
        const Numeric term_value = (j & 1U ? -choose : choose) * moment[n + j];
        const Numeric yk         = term_value - correction_a;
        const Numeric next       = sum_a + yk;
        correction_a             = (next - sum_a) - yk;
        sum_a                    = next;

        if (n <= order) {
          const Numeric b_value = (j & 1U ? -choose : choose) * moment[n + j + 1];
          const Numeric b_yk    = b_value - correction_b;
          const Numeric b_next  = sum_b + b_yk;
          correction_b          = (b_next - sum_b) - b_yk;
          sum_b                 = b_next;
        }

        if (j != n) choose *= static_cast<Numeric>(n - j) / static_cast<Numeric>(j + 1);
      }
      A[n] = sum_a;
      if (n <= order) Bm[n] = sum_b;
    }

    Numeric    value{}, value_m1{}, dtau{}, ddelta{}, coefficient{1.0};
    Numeric    value_correction{}, value_m1_correction{}, dtau_correction{}, ddelta_correction{};
    const auto compensated_add = [](Numeric &sum, Numeric &correction, const Numeric addend) {
      const Numeric yk   = addend - correction;
      const Numeric next = sum + yk;
      correction         = (next - sum) - yk;
      sum                = next;
    };
    for (Size n = 0; n <= order; ++n) {
      compensated_add(value, value_correction, coefficient * A[n]);
      compensated_add(value_m1,
                      value_m1_correction,
                      n == 0 ? phi1m1_scalar_from_expm1(minus_tau, exp_minus_tau_m1) : coefficient * A[n]);
      compensated_add(dtau, dtau_correction, -coefficient * Bm[n]);
      compensated_add(ddelta, ddelta_correction, -coefficient * A[n + 1]);
      coefficient *= -delta / static_cast<Numeric>(n + 1);
    }
    return {value, value_m1, dtau, ddelta};
  }

  Numeric value;
  if (delta > 0.0) {
    const Numeric root = std::sqrt(delta);
    value              = (Faddeeva::Dawson(end_optical / (2.0 * root)) -
                          exp_minus_tau * Faddeeva::Dawson(start_optical / (2.0 * root))) /
                         root;
  } else {
    const Numeric gamma = -delta;
    const Numeric root  = std::sqrt(gamma);
    value =
        0.5 * Constant::sqrt_pi / root *
        (Faddeeva::erfcx(end_optical / (2.0 * root)) - exp_minus_tau * Faddeeva::erfcx(start_optical / (2.0 * root)));
  }

  const Numeric first_moment  = (exp_minus_tau_m1 + beta * value) / (2.0 * delta);
  const Numeric second_moment = (exp_minus_tau - value + beta * first_moment) / (2.0 * delta);
  return {value, value - 1.0, -first_moment, second_moment - first_moment};
}

Numeric scalar_linprop_deriv(const scalar_linprop_result &state,
                             const Numeric                ks,
                             const Numeric                ke,
                             const Numeric                dk,
                             const Numeric                r,
                             const Numeric                dr,
                             const bool                   start_deriv) {
  // Scale endpoint combinations before adding them.  Besides avoiding the
  // obvious ks+ke overflow, this preserves finite derivatives when r or dr
  // is the small factor that makes very large endpoint matrices physical.
  const Numeric dtau   = dr * std::midpoint(ks, ke) + 0.5 * r * dk;
  const Numeric ddelta = dr * std::midpoint(ke, -ks) + 0.5 * r * (start_deriv ? -dk : dk);
  return state.dtau * dtau + state.ddelta * ddelta;
}

using diagonal_remainder = Vector4;

void set_diagonal_remainder(stokvec_matrix &out, const Size iv, const Size ip, const diagonal_remainder &diag) {
  out[iv, ip] = diag.data;
}

muelmat augmented_linprop_value(const tran                  &tr,
                                const scalar_linprop_result &state,
                                const propmat               &k1,
                                const propmat               &k2,
                                const Numeric                r,
                                diagonal_remainder          *diag_m1   = nullptr,
                                std::optional<muelmat>      *base_phi  = nullptr,
                                std::optional<muelmat>      *base_q_m1 = nullptr) {
  const Numeric      phi_m1              = phi1m1_scalar_from_expm1(tr.a, tr.expm1_a);
  const propmat      delta_r             = r * k2 - r * k1;
  const muelmat      q_m1                = -(1.0 / 12.0) * to_muelmat(delta_r);
  const Numeric      q_s_m1              = -delta_r.A() / 12.0;
  const Numeric      augmented_scalar_m1 = phi_m1 + q_s_m1 + phi_m1 * q_s_m1;
  const Numeric      scalar_correction   = state.value_m1 - augmented_scalar_m1;
  diagonal_remainder f_diag_m1;
  const muelmat      f          = diag_m1 == nullptr ? tr.linsrc() : tr.linsrc(f_diag_m1);
  const muelmat      correction = f * q_m1;

  if (base_phi != nullptr) base_phi->emplace(f);
  if (base_q_m1 != nullptr) base_q_m1->emplace(q_m1);

  if (diag_m1 != nullptr) {
    *diag_m1 = f_diag_m1;
    for (Size is = 0; is < 4; ++is) (*diag_m1)[is] += correction[is, is] + scalar_correction;
  }

  return f + correction + scalar_correction * muelmat::id();
}

muelmat augmented_magnus_source(const tran             &tr,
                                const propmat          &k1,
                                const propmat          &k2,
                                const Numeric           r,
                                diagonal_remainder     *diag_m1   = nullptr,
                                std::optional<muelmat> *base_phi  = nullptr,
                                std::optional<muelmat> *base_q_m1 = nullptr) {
  const propmat      delta_r = r * k2 - r * k1;
  const muelmat      q_m1    = -(1.0 / 12.0) * to_muelmat(delta_r);
  diagonal_remainder f_diag_m1;
  const muelmat      f          = diag_m1 == nullptr ? tr.linsrc() : tr.linsrc(f_diag_m1);
  const muelmat      correction = f * q_m1;

  if (base_phi != nullptr) base_phi->emplace(f);
  if (base_q_m1 != nullptr) base_q_m1->emplace(q_m1);

  if (diag_m1 != nullptr) {
    *diag_m1 = f_diag_m1;
    for (Size is = 0; is < 4; ++is) (*diag_m1)[is] += correction[is, is];
  }

  return f + correction;
}

muelmat augmented_linprop_deriv(const tran                  &tr,
                                const scalar_linprop_result &state,
                                const muelmat               &f,
                                const muelmat               &q_m1,
                                const propmat               &kavg,
                                const Numeric                phi,
                                const Numeric                phi_prime,
                                const propmat               &k1,
                                const propmat               &k2,
                                const propmat               &dk,
                                const Numeric                r,
                                const Numeric                dr,
                                const bool                   k1_deriv) {
  const Numeric dvalue   = scalar_linprop_deriv(state, k1.A(), k2.A(), dk.A(), r, dr, k1_deriv);
  const propmat dg       = -dr * kavg - 0.5 * r * dk;
  const Numeric dphi     = phi_prime * dg.A();
  const propmat ddelta_r = dr * k2 - dr * k1 + (k1_deriv ? -r * dk : r * dk);
  const muelmat dq       = -(1.0 / 12.0) * to_muelmat(ddelta_r);
  const Numeric q_s      = 1.0 + q_m1[0, 0];
  const Numeric dq_s     = -ddelta_r.A() / 12.0;
  const Numeric dc       = dvalue - (dphi * q_s + phi * dq_s);
  const muelmat df       = tr.linsrc_deriv(dg);
  return df + df * q_m1 + f * dq + dc * muelmat::id();
}

muelmat augmented_magnus_source_deriv(const tran    &tr,
                                      const muelmat &f,
                                      const muelmat &q_m1,
                                      const propmat &k1,
                                      const propmat &k2,
                                      const propmat &dk,
                                      const Numeric  r,
                                      const Numeric  dr,
                                      const bool     k1_deriv) {
  const propmat ddelta_r = dr * k2 - dr * k1 + (k1_deriv ? -r * dk : r * dk);
  const propmat dg       = magnus_exponent_deriv(k1, k2, dk, r, dr, k1_deriv);
  const muelmat dq       = -(1.0 / 12.0) * to_muelmat(ddelta_r);
  const muelmat df       = tr.linsrc_deriv(dg);
  return df + df * q_m1 + f * dq;
}
}  // namespace

tran::tran(const propmat &k1, const propmat &k2, const Numeric r)
    : a{-r * std::midpoint(k1.A(), k2.A())}, polarized(k1.is_polarized() or k2.is_polarized()) {
  if (std::abs(a) < 0.5) {
    expm1_a = std::expm1(a);
    exp_a   = 1.0 + expm1_a;
  } else {
    exp_a   = std::exp(a);
    expm1_a = exp_a - 1.0;
  }
  if (not polarized) return;

  b = -r * std::midpoint(k1.B(), k2.B());
  c = -r * std::midpoint(k1.C(), k2.C());
  d = -r * std::midpoint(k1.D(), k2.D());
  u = -r * std::midpoint(k1.U(), k2.U());
  v = -r * std::midpoint(k1.V(), k2.V());
  w = -r * std::midpoint(k1.W(), k2.W());

  init_polarized();
}

tran::tran(const propmat &k1, const propmat &k2, const Numeric r, MagnusOperator)
    : a{-r * std::midpoint(k1.A(), k2.A())}, polarized(k1.is_polarized() or k2.is_polarized()) {
  if (std::abs(a) < 0.5) {
    expm1_a = std::expm1(a);
    exp_a   = 1.0 + expm1_a;
  } else {
    exp_a   = std::exp(a);
    expm1_a = exp_a - 1.0;
  }
  if (not polarized) return;

  // Scale before forming the commutator.  This evaluates r^2 [K2,K1]
  // without overflowing products that the path-length factors later reduce.
  const propmat correction = (1.0 / 12.0) * commutator(r * k2, r * k1);
  b                        = -r * std::midpoint(k1.B(), k2.B()) + correction.B();
  c                        = -r * std::midpoint(k1.C(), k2.C()) + correction.C();
  d                        = -r * std::midpoint(k1.D(), k2.D()) + correction.D();
  u                        = -r * std::midpoint(k1.U(), k2.U()) + correction.U();
  v                        = -r * std::midpoint(k1.V(), k2.V()) + correction.V();
  w                        = -r * std::midpoint(k1.W(), k2.W()) + correction.W();

  init_polarized();
}

void tran::init_polarized() {
  b2 = b * b;
  c2 = c * c;
  d2 = d * d;
  u2 = u * u;
  v2 = v * v;
  w2 = w * w;

  /** Solve:
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
  B               = u2 + v2 + w2 - b2 - c2 - d2;
  const Numeric q = d * u - c * v + b * w;
  C               = -Math::pow2(q);
  S               = std::hypot(B, 2.0 * q);

  /**
        We define:
            x2 and y2 are the squares of x and y
            x and y are real and positive
            x is from the real part of the Eigenvalues
            y is from the imag part of the Eigenvalues
        Notes:
             S  >=  0
            |S| >= |B|
            S-B >=  0
            S+B >=  0
            y2 stores the magnitude of the negative squared eigenvalue
    */
  // The characteristic equation gives lambda^2 = (-B +/- S) / 2.
  // Therefore these quantities already are x^2 and -y^2.  In particular,
  // no additional square root belongs here; x and y are taken below.
  if (B >= 0.0) {
    y2 = 0.5 * (S + B);
    x2 = y2 == 0.0 ? 0.0 : (std::abs(q) / y2) * std::abs(q);
  } else {
    x2 = 0.5 * (S - B);
    y2 = x2 == 0.0 ? 0.0 : (std::abs(q) / x2) * std::abs(q);
  }
  x = std::sqrt(x2);
  y = std::sqrt(y2);

  x_zero      = x < too_small;
  y_zero      = y < too_small;
  both_zero   = y_zero and x_zero;
  either_zero = y_zero or x_zero;

  /* Using:
     *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
     *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
     *    inv_x2 := 1 for x == 0,
     *    -i sin(ix) → sinh(x)
     *    cos(ix) → cosh(x)
     *    C0, C1, C2 ∝ [1/x^2]
     */
  ix = x_zero ? 0.0 : 1.0 / x;
  iy = y_zero ? 0.0 : 1.0 / y;

  // The first "1.0" is the trick for above limits
  inv_x2y2 = both_zero ? 1.0 : 1.0 / (x2 + y2);

  if (S < 1e-4) {
    const Numeric x4 = x2 * x2;
    const Numeric y4 = y2 * y2;
    const Numeric p4 = x4 - x2 * y2 + y4;
    const Numeric p6 = x4 * x2 - x4 * y2 + x2 * y4 - y4 * y2;

    // At this point x,y < 0.01.  Polynomial evaluation is both faster than
    // four library calls and accurate through round-off at this scale.
    cx = 1.0 + x2 * (0.5 + x2 * (1.0 / 24.0 + x2 * (1.0 / 720.0 + x2 / 40320.0)));
    sx = x * (1.0 + x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 * (1.0 / 5040.0 + x2 / 362880.0))));
    cy = 1.0 + y2 * (-0.5 + y2 * (1.0 / 24.0 + y2 * (-1.0 / 720.0 + y2 / 40320.0)));
    sy = y * (1.0 + y2 * (-1.0 / 6.0 + y2 * (1.0 / 120.0 + y2 * (-1.0 / 5040.0 + y2 / 362880.0))));

    C0_m1 = x2 * y2 * (1.0 / 24.0 + (x2 - y2) / 720.0 + p4 / 40320.0);
    C0    = 1.0 + C0_m1;
    C1    = 1.0 + x2 * y2 / 120.0 + x2 * y2 * (x2 - y2) / 5040.0;
    C2    = 0.5 + (x2 - y2) / 24.0 + p4 / 720.0 + p6 / 40320.0;
    C3    = 1.0 / 6.0 + (x2 - y2) / 120.0 + p4 / 5040.0 + p6 / 362880.0;
  } else {
    cy = std::cos(y);
    sy = std::sin(y);
    cx = std::cosh(x);
    sx = std::sinh(x);

    const auto coshm1 = [](const Numeric value, const Numeric z2) {
      return z2 < 0.01 ? z2 * (0.5 + z2 * (1.0 / 24.0 + z2 * (1.0 / 720.0 + z2 / 40320.0))) : value - 1.0;
    };
    const auto cosm1 = [](const Numeric value, const Numeric z2) {
      return z2 < 0.01 ? z2 * (-0.5 + z2 * (1.0 / 24.0 + z2 * (-1.0 / 720.0 + z2 / 40320.0))) : value - 1.0;
    };
    const auto sinhcm1 = [](const Numeric value, const Numeric z, const Numeric z2) {
      return z2 < 0.01 ? z2 * (1.0 / 6.0 + z2 * (1.0 / 120.0 + z2 * (1.0 / 5040.0 + z2 / 362880.0))) : value / z - 1.0;
    };
    const auto sincm1 = [](const Numeric value, const Numeric z, const Numeric z2) {
      return z2 < 0.01 ? z2 * (-1.0 / 6.0 + z2 * (1.0 / 120.0 + z2 * (-1.0 / 5040.0 + z2 / 362880.0)))
                       : value / z - 1.0;
    };

    const Numeric coshm1_x  = coshm1(cx, x2);
    const Numeric cosm1_y   = cosm1(cy, y2);
    const Numeric sinhcm1_x = sinhcm1(sx, x, x2);
    const Numeric sincm1_y  = sincm1(sy, y, y2);
    C0_m1                   = (x2 * cosm1_y + y2 * coshm1_x) * inv_x2y2;
    C0                      = 1.0 + C0_m1;
    C1                      = 1.0 + (x2 * sincm1_y + y2 * sinhcm1_x) * inv_x2y2;
    C2                      = (coshm1_x - cosm1_y) * inv_x2y2;
    C3                      = (sinhcm1_x - sincm1_y) * inv_x2y2;
  }

  stable_coefficients =
      std::isfinite(C0_m1) and std::isfinite(C0) and std::isfinite(C1) and std::isfinite(C2) and std::isfinite(C3);
  const Numeric polarization_norm =
      std::max({std::abs(b), std::abs(c), std::abs(d), std::abs(u), std::abs(v), std::abs(w)});
  spectrally_clustered = S > 0.0 and S <= spectral_cluster_relative * polarization_norm * polarization_norm;
}

namespace {
muelmat scaled_spectral_exponential(const tran &tr) {
  // Form exp(a)cosh(x) and exp(a)sinh(x) from the actual real eigenmodes.
  // This preserves a finite a+x mode when exp(a) underflows and cosh(x)
  // overflows, and avoids the loss accumulated by many matrix squarings.
  const Numeric exp_plus  = std::exp(tr.a + tr.x);
  const Numeric exp_minus = std::exp(tr.a - tr.x);
  const Numeric hc        = std::midpoint(exp_plus, exp_minus);
  const Numeric hs        = 0.5 * exp_plus - 0.5 * exp_minus;
  const Numeric tc        = tr.exp_a * tr.cy;
  const Numeric ts        = tr.exp_a * tr.sy;

  const Numeric wx        = tr.x2 / tr.S;
  const Numeric wy        = tr.y2 / tr.S;
  const Numeric hs_over_x = tr.x == 0.0 ? tr.exp_a : hs / tr.x;
  const Numeric ts_over_y = tr.y == 0.0 ? tr.exp_a : ts / tr.y;
  const Numeric D0        = wx * tc + wy * hc;
  const Numeric D1        = wx * ts_over_y + wy * hs_over_x;
  const Numeric D2        = (hc - tc) / tr.S;
  const Numeric D3        = (hs_over_x - ts_over_y) / tr.S;

  const Numeric D2b = D2 * (tr.c * tr.u + tr.d * tr.v);
  const Numeric D2c = D2 * (tr.b * tr.u - tr.d * tr.w);
  const Numeric D2d = D2 * (tr.b * tr.v + tr.c * tr.w);
  const Numeric D2u = D2 * (tr.b * tr.c - tr.v * tr.w);
  const Numeric D2v = D2 * (tr.b * tr.d + tr.u * tr.w);
  const Numeric D2w = D2 * (tr.c * tr.d - tr.u * tr.v);

  const Numeric D3b = D3 * (tr.b * (tr.B - tr.w2) + tr.w * (tr.c * tr.v - tr.d * tr.u));
  const Numeric D3c = D3 * (tr.c * (tr.v2 - tr.B) - tr.v * (tr.d * tr.u + tr.b * tr.w));
  const Numeric D3d = D3 * (tr.d * (tr.u2 - tr.B) - tr.u * (tr.c * tr.v - tr.b * tr.w));
  const Numeric D3u = D3 * (tr.d * (tr.c * tr.v - tr.b * tr.w) - tr.u * (tr.B + tr.d2));
  const Numeric D3v = D3 * (tr.c * (tr.d * tr.u + tr.b * tr.w) - tr.v * (tr.B + tr.c2));
  const Numeric D3w = D3 * (tr.b * (tr.c * tr.v - tr.d * tr.u) - tr.w * (tr.B + tr.b2));

  const Numeric M00 = D0 + D2 * (tr.b2 + tr.c2 + tr.d2);
  const Numeric M11 = D0 + D2 * (tr.b2 - tr.u2 - tr.v2);
  const Numeric M22 = D0 + D2 * (tr.c2 - tr.u2 - tr.w2);
  const Numeric M33 = D0 + D2 * (tr.d2 - tr.v2 - tr.w2);

  // clang-format off
  return muelmat{
    M00,                        D1 * tr.b - D2b - D3b,     D1 * tr.c + D2c + D3c,     D1 * tr.d + D2d + D3d,
    D1 * tr.b + D2b - D3b,     M11,                        D1 * tr.u + D2u + D3u,     D1 * tr.v + D2v + D3v,
    D1 * tr.c - D2c + D3c,   - D1 * tr.u + D2u - D3u,     M22,                        D1 * tr.w + D2w + D3w,
    D1 * tr.d - D2d + D3d,   - D1 * tr.v + D2v - D3v,   - D1 * tr.w + D2w - D3w,     M33};
  // clang-format on
}

bool use_scaled_spectral_exponential(const tran &tr) {
  return tr.S > 0.0 and std::isfinite(tr.S) and std::isfinite(tr.x) and std::isfinite(tr.y) and
         (tr.exp_a < std::numeric_limits<Numeric>::min() or tr.x > 500.0);
}

bool use_degenerate_exponential(const tran &tr) {
  // The explicit small-S polynomials are already the well-conditioned
  // clustered limit, and are substantially cheaper than matrix Taylor
  // scaling.  Larger clustered spectra still need the generic fallback.
  return not tr.stable_coefficients or (tr.spectrally_clustered and tr.S >= 1e-4) or
         tr.exp_a < std::numeric_limits<Numeric>::min();
}

bool use_real_pair_frechet_limit(const tran &tr) {
  if (tr.y == 0.0) return true;
  if (tr.x <= 2.0 or not std::isfinite(tr.S) or not std::isfinite(tr.x2) or not std::isfinite(tr.y2)) return false;

  // The near-real expansion below retains every term linear in y/x, so its
  // neglected terms are at round-off once (y/x)^2 <= epsilon.  Comparing
  // the already available squared roots adds no work to the common path.
  return tr.y2 <= std::numeric_limits<Numeric>::epsilon() * tr.x2;
}

enum class frechet_function { exponential, phi1 };

muelmat separated_real_frechet(const tran &tr, const propmat &dg, const frechet_function function) {
  // Separate the +/-x eigenspaces from the nearly repeated imaginary pair.
  // The real projectors are exact also for nonzero y; the imaginary block is
  // retained through first order in y/x below.
  const muelmat N{0, tr.b, tr.c, tr.d, tr.b, 0, tr.u, tr.v, tr.c, -tr.u, 0, tr.w, tr.d, -tr.v, -tr.w, 0};
  const muelmat I   = muelmat::id();
  const muelmat R   = N / tr.x;
  const Numeric rho = tr.y2 / tr.x2;
  const muelmat Pr  = (R * R + rho * I) / (1.0 + rho);
  const muelmat NPr = tr.y == 0.0 ? N : N * Pr;
  const muelmat Rr  = NPr / tr.x;
  const muelmat Pp  = 0.5 * (Pr + Rr);
  const muelmat Pm  = 0.5 * (Pr - Rr);
  const muelmat P0  = I - Pr;

  const Numeric lp  = tr.a + tr.x;
  const Numeric lm  = tr.a - tr.x;
  const Numeric l0  = tr.a;
  const Numeric vp  = function == frechet_function::exponential ? std::exp(lp) : phi1_scalar(lp);
  const Numeric vm  = function == frechet_function::exponential ? std::exp(lm) : phi1_scalar(lm);
  const Numeric v0  = function == frechet_function::exponential ? tr.exp_a : (l0 == 0.0 ? 1.0 : tr.expm1_a / l0);
  const Numeric dp  = function == frechet_function::exponential ? vp : phi1_scalar_derivative<1>(lp);
  const Numeric dm  = function == frechet_function::exponential ? vm : phi1_scalar_derivative<1>(lm);
  const Numeric d0  = function == frechet_function::exponential ? v0 : phi1_scalar_derivative_from_exp(l0, tr.exp_a);
  const Numeric dd0 = function == frechet_function::exponential ? v0 : phi1_scalar_derivative<2>(l0);
  const auto    divided_difference = [](const Numeric lhs,
                                        const Numeric rhs,
                                        const Numeric lhs_value,
                                        const Numeric rhs_value,
                                        const Numeric common_derivative) {
    if (lhs == rhs) return common_derivative;
    return (lhs_value - rhs_value) / (lhs - rhs);
  };

  const muelmat E         = to_muelmat(dg);
  const Numeric e_largest = max_abs(E);
  if (e_largest == 0.0) return muelmat::constant(0.0);

  // A power-of-two normalization is exact.  Unlike division by the largest
  // entry, it preserves a small eigenvalue perturbation such as da+db=-10
  // when da and db are O(1e12).
  int e_exponent{};
  std::frexp(e_largest, &e_exponent);
  const Numeric e_scale = std::scalbn(1.0, e_exponent - 1);
  const muelmat En      = E / e_scale;

  const muelmat PpE = Pp * En;
  const muelmat PmE = Pm * En;
  const muelmat P0E = P0 * En;
  const muelmat Epp = PpE * Pp;
  const muelmat Epm = PpE * Pm;
  const muelmat Ep0 = PpE * P0;
  const muelmat Emp = PmE * Pp;
  const muelmat Emm = PmE * Pm;
  const muelmat Em0 = PmE * P0;
  const muelmat E0p = P0E * Pp;
  const muelmat E0m = P0E * Pm;
  const muelmat E00 = P0E * P0;

  const Numeric cpm    = divided_difference(lp, lm, vp, vm, dp);
  const Numeric cp0    = divided_difference(lp, l0, vp, v0, dp);
  const Numeric cm0    = divided_difference(lm, l0, vm, v0, dm);
  muelmat       result = dp * Epp + dm * Emm + d0 * E00 + cpm * (Epm + Emp) + cp0 * (Ep0 + E0p) + cm0 * (Em0 + E0m);

  if (tr.y != 0.0) {
    const muelmat Z    = N - NPr;
    const muelmat ZE   = Z * En;
    const Numeric hp1  = (cp0 - d0) / (lp - l0);
    const Numeric hm1  = (cm0 - d0) / (lm - l0);
    result            += hp1 * (PpE * Z + ZE * Pp) + hm1 * (PmE * Z + ZE * Pm) + 0.5 * dd0 * (Z * E00 + E00 * Z);
  }

  return e_scale * result;
}
}  // namespace

muelmat tran::operator()() const noexcept {
  if (not polarized) return exp_a;

  if (use_scaled_spectral_exponential(*this)) return scaled_spectral_exponential(*this);

  if (use_degenerate_exponential(*this)) { return exp_degenerate(propmat{a, b, c, d, u, v, w}); }

  // Do the calculation exp(a) * (C0 * I + C1 * K + C2 * K^2 + C3 * K^3)

  // Scale the coefficients before forming powers of the polarized part.  This
  // avoids overflowing an intermediate Cj*K^j when exp(a) later brings the
  // physical transmission back into range.
  const Numeric D0 = exp_a * C0;
  const Numeric D1 = exp_a * C1;
  const Numeric D2 = exp_a * C2;
  const Numeric D3 = exp_a * C3;

  const Numeric D2b = D2 * (c * u + d * v);
  const Numeric D2c = D2 * (b * u - d * w);
  const Numeric D2d = D2 * (b * v + c * w);
  const Numeric D2u = D2 * (b * c - v * w);
  const Numeric D2v = D2 * (b * d + u * w);
  const Numeric D2w = D2 * (c * d - u * v);

  // B = u2 + v2 + w2 - b2 - c2 - d2
  const Numeric D3b = D3 * (b * (B - w2) + w * (c * v - d * u));
  const Numeric D3c = D3 * (c * (v2 - B) - v * (d * u + b * w));
  const Numeric D3d = D3 * (d * (u2 - B) - u * (c * v - b * w));
  const Numeric D3u = D3 * (d * (c * v - b * w) - u * (B + d2));
  const Numeric D3v = D3 * (c * (d * u + b * w) - v * (B + c2));
  const Numeric D3w = D3 * (b * (c * v - d * u) - w * (B + b2));

  const Numeric M00 = D0 + D2 * (b2 + c2 + d2);
  const Numeric M11 = D0 + D2 * (b2 - u2 - v2);
  const Numeric M22 = D0 + D2 * (c2 - u2 - w2);
  const Numeric M33 = D0 + D2 * (d2 - v2 - w2);

  // clang-format off
  return muelmat{
    M00,                  D1 * b - D2b - D3b,   D1 * c + D2c + D3c, D1 * d + D2d + D3d,
    D1 * b + D2b - D3b,   M11,                  D1 * u + D2u + D3u, D1 * v + D2v + D3v,
    D1 * c - D2c + D3c, - D1 * u + D2u - D3u,   M22,                D1 * w + D2w + D3w,
    D1 * d - D2d + D3d, - D1 * v + D2v - D3v, - D1 * w + D2w - D3w, M33};
  // clang-format on
}

muelmat tran::operator()(Vector4 &diag_m1) const noexcept {
  if (not polarized) {
    diag_m1.fill(expm1_a);
    return {exp_a};
  }

  if (use_scaled_spectral_exponential(*this)) {
    const muelmat value = scaled_spectral_exponential(*this);
    diag_m1             = {value[0, 0] - 1.0, value[1, 1] - 1.0, value[2, 2] - 1.0, value[3, 3] - 1.0};
    return value;
  }

  if (use_degenerate_exponential(*this)) {
    const propmat g{a, b, c, d, u, v, w};
    if (std::abs(a) + polarized_row_norm(g) < 0.5) {
      muelmat value = expm1_degenerate(g);
      diag_m1       = {value[0, 0], value[1, 1], value[2, 2], value[3, 3]};
      for (Size is = 0; is < 4; ++is) value[is, is] += 1.0;
      return value;
    }

    const muelmat value = exp_degenerate(g);
    diag_m1             = {value[0, 0] - 1.0, value[1, 1] - 1.0, value[2, 2] - 1.0, value[3, 3] - 1.0};
    return value;
  }

  const muelmat value = (*this)();
  diag_m1             = expm1_diagonal();
  return value;
}

muelmat tran::expm1() const noexcept {
  if (not polarized) return {expm1_a};

  if (use_scaled_spectral_exponential(*this)) {
    muelmat value = scaled_spectral_exponential(*this);
    for (Size is = 0; is < 4; ++is) value[is, is] -= 1.0;
    return value;
  }

  if (use_degenerate_exponential(*this)) { return expm1_degenerate(propmat{a, b, c, d, u, v, w}); }

  // Compute the diagonal offset: e^a * C0 - 1
  //    = (expm1(a) + 1) * (C0_m1 + 1) - 1
  //    = expm1(a) * (C0_m1 + 1) + C0_m1
  const Numeric em1_a       = expm1_a;
  const Numeric scaled_c0   = exp_a * C0;
  const Numeric direct_diag = scaled_c0 - 1.0;
  const Numeric robust_diag = em1_a * (C0_m1 + 1.0) + C0_m1;
  const Numeric diag_offset = std::abs(a) < 0.5 and std::abs(C0_m1) < 0.5 ? robust_diag : direct_diag;

  // Compute the rest of the matrix terms
  //    The result is: diag_offset * I + exp_a * (C1*K + C2*K^2 + C3*K^3)
  //    Note that K and K^3 have zero diagonals, but K^2 has diagonal terms.

  const Numeric D1 = exp_a * C1;
  const Numeric D2 = exp_a * C2;
  const Numeric D3 = exp_a * C3;

  const Numeric D2b = D2 * (c * u + d * v);
  const Numeric D2c = D2 * (b * u - d * w);
  const Numeric D2d = D2 * (b * v + c * w);
  const Numeric D2u = D2 * (b * c - v * w);
  const Numeric D2v = D2 * (b * d + u * w);
  const Numeric D2w = D2 * (c * d - u * v);

  const Numeric D3b = D3 * (b * (B - w2) + w * (c * v - d * u));
  const Numeric D3c = D3 * (c * (v2 - B) - v * (d * u + b * w));
  const Numeric D3d = D3 * (d * (u2 - B) - u * (c * v - b * w));
  const Numeric D3u = D3 * (d * (c * v - b * w) - u * (B + d2));
  const Numeric D3v = D3 * (c * (d * u + b * w) - v * (B + c2));
  const Numeric D3w = D3 * (b * (c * v - d * u) - w * (B + b2));

  // Diagonal contributions from K^2 term
  const Numeric M00_rest = D2 * (b2 + c2 + d2);
  const Numeric M11_rest = D2 * (b2 - u2 - v2);
  const Numeric M22_rest = D2 * (c2 - u2 - w2);
  const Numeric M33_rest = D2 * (d2 - v2 - w2);

  // clang-format off
  return muelmat{
    diag_offset + M00_rest,       D1 * b - D2b - D3b,         D1 * c + D2c + D3c,       D1 * d + D2d + D3d,
    D1 * b + D2b - D3b,           diag_offset + M11_rest,     D1 * u + D2u + D3u,       D1 * v + D2v + D3v,
    D1 * c - D2c + D3c,           -D1 * u + D2u - D3u,        diag_offset + M22_rest,     D1 * w + D2w + D3w,
    D1 * d - D2d + D3d,           -D1 * v + D2v - D3v,        -D1 * w + D2w - D3w,      diag_offset + M33_rest
  };
  // clang-format on
}

Vector4 tran::expm1_diagonal() const noexcept {
  if (not polarized) {
    const Numeric value = expm1_a;
    return {value, value, value, value};
  }

  if (use_scaled_spectral_exponential(*this)) {
    const muelmat value = scaled_spectral_exponential(*this);
    return {value[0, 0] - 1.0, value[1, 1] - 1.0, value[2, 2] - 1.0, value[3, 3] - 1.0};
  }

  if (use_degenerate_exponential(*this)) {
    const muelmat value = expm1_degenerate(propmat{a, b, c, d, u, v, w});
    return {value[0, 0], value[1, 1], value[2, 2], value[3, 3]};
  }

  // The odd powers of the polarized generator have zero diagonal.  Retain
  // exp(a) C0 - 1 without first rounding exp(a) C0 to one, then add the four
  // diagonal entries of exp(a) C2 N^2.
  const Numeric em1_a       = expm1_a;
  const Numeric direct_diag = exp_a * C0 - 1.0;
  const Numeric robust_diag = em1_a * (C0_m1 + 1.0) + C0_m1;
  const Numeric diag_offset = std::abs(a) < 0.5 and std::abs(C0_m1) < 0.5 ? robust_diag : direct_diag;
  const Numeric D2          = exp_a * C2;

  return {diag_offset + D2 * (b2 + c2 + d2),
          diag_offset + D2 * (b2 - u2 - v2),
          diag_offset + D2 * (c2 - u2 - w2),
          diag_offset + D2 * (d2 - v2 - w2)};
}

muelmat tran::linsrc() const noexcept { return linsrc_impl(nullptr); }

muelmat tran::linsrc(Vector4 &diag_m1) const noexcept { return linsrc_impl(&diag_m1); }

muelmat tran::linsrc_impl(Vector4 *diag_m1) const noexcept {
  const Numeric phi1_a = a == 0.0 ? 1.0 : expm1_a / a;
  if (not polarized) {
    if (diag_m1 != nullptr) diag_m1->fill(phi1m1_scalar_from_expm1(a, expm1_a));
    return {phi1_a};
  }

  const propmat g{a, b, c, d, u, v, w};
  const Numeric generator_norm        = std::abs(a) + polarized_row_norm(g);
  const bool    usable_spectral_roots = S > 0.0 and std::isfinite(S) and std::isfinite(x) and std::isfinite(y) and
                                        std::isfinite(x2) and std::isfinite(y2);
  if ((not stable_coefficients or exp_a < std::numeric_limits<Numeric>::min()) and not usable_spectral_roots) {
    const muelmat value = phi1_degenerate(g);
    if (diag_m1 != nullptr) *diag_m1 = linsrcm1_diagonal(value);
    return value;
  }

  bool phi1_ill_conditioned = x_zero or y_zero or spectrally_clustered;
  if (not phi1_ill_conditioned) {
    const Numeric polarization_scale =
        std::max({std::abs(b), std::abs(c), std::abs(d), std::abs(u), std::abs(v), std::abs(w)});
    phi1_ill_conditioned = S > 0.0 and S <= 0.01 * polarization_scale * polarization_scale;
  }
  if (not phi1_ill_conditioned)
    phi1_ill_conditioned = std::abs(a + x) < 0.01 or std::abs(a - x) < 0.01 or a * a + y * y < 1e-4;

  // The divided differences below subtract functions evaluated at clustered
  // eigenvalues and then divide by x^2+y^2 (= S).  The specialized
  // Cayley-Hamilton series is accurate and cheap for the small/medium
  // generators where this occurs; scaling and squaring covers larger cases.
  if (phi1_ill_conditioned) {
    if (generator_norm <= 8.0) return phi1_cayley_hamilton_series(g, diag_m1);

    // An exactly zero x or y is handled by explicit analytic limits below.
    // When the other root is well separated, those limits are preferable to
    // scaling-and-squaring, which can erase a weak eigenmode in an extremely
    // opaque layer.
    const Numeric polarization_scale =
        std::max({std::abs(b), std::abs(c), std::abs(d), std::abs(u), std::abs(v), std::abs(w)});
    const bool separated_zero_limit =
        (x_zero != y_zero) and S > 0.01 * polarization_scale * polarization_scale and usable_spectral_roots;
    if (not separated_zero_limit) {
      const muelmat value = phi1_degenerate(g);
      if (diag_m1 != nullptr) *diag_m1 = linsrcm1_diagonal(value);
      return value;
    }
  }

  const Numeric phi1p_a = (x_zero or y_zero) ? phi1_scalar_derivative_from_exp(a, exp_a) : 0.0;
  Numeric       l0, l1, l2, l3;

  if (both_zero) {
    l0 = phi1_a;
    l1 = phi1p_a;
    l2 = 0.5 * phi1_scalar_derivative<2>(a);
    l3 = phi1_scalar_derivative<3>(a) / 6.0;
  } else {
    Numeric Pp, Pm_div_x;
    if (x_zero) {
      Pp       = phi1_a;
      Pm_div_x = phi1p_a;
    } else {
      const Numeric f1 = phi1_scalar(a + x);
      const Numeric f2 = phi1_scalar(a - x);
      Pp               = 0.5 * (f1 + f2);
      Pm_div_x         = 0.5 * (f1 - f2) / x;
    }

    Numeric Qp, q_im;
    if (y_zero) {
      Qp   = phi1_a;
      q_im = phi1p_a;
    } else {
      const Numeric denom       = a * a + y * y;
      const Numeric ea_cy_m1    = exp_a * cy - 1.0;
      const Numeric ea_sy       = exp_a * sy;
      Qp                        = (a * ea_cy_m1 + y * ea_sy) / denom;
      const Numeric sin_y_div_y = (std::abs(y) < 1e-6) ? 1.0 - y * y / 6.0 : sy / y;
      q_im                      = (a * exp_a * sin_y_div_y - ea_cy_m1) / denom;
    }

    l2 = (Pp - Qp) * inv_x2y2;
    l0 = Pp - l2 * x2;
    l3 = (Pm_div_x - q_im) * inv_x2y2;
    l1 = Pm_div_x - l3 * x2;
  }

  const muelmat S_mat{0, b, c, d, b, 0, u, v, c, -u, 0, w, d, -v, -w, 0};

  const muelmat S2_mat = S_mat * S_mat;
  const muelmat S3_mat = S_mat * S2_mat;

  const muelmat value =
      muelmat{l0, 0, 0, 0, 0, l0, 0, 0, 0, 0, l0, 0, 0, 0, 0, l0} + S_mat * l1 + S2_mat * l2 + S3_mat * l3;
  if (diag_m1 != nullptr) *diag_m1 = {value[0, 0] - 1.0, value[1, 1] - 1.0, value[2, 2] - 1.0, value[3, 3] - 1.0};
  return value;
}

Vector4 tran::linsrcm1_diagonal(const muelmat &lambda) const noexcept {
  if (not polarized) {
    const Numeric value = phi1m1_scalar_from_expm1(a, expm1_a);
    return {value, value, value, value};
  }

  const propmat g{a, b, c, d, u, v, w};
  const Numeric n = polarized_row_norm(g);
  if (std::abs(a) + n > 0.5) {
    return {lambda[0, 0] - 1.0, lambda[1, 1] - 1.0, lambda[2, 2] - 1.0, lambda[3, 3] - 1.0};
  }

  // For a thin layer, form only the Cayley-Hamilton coefficients needed by
  // the diagonal of phi_1(G)-I.  This costs scalar arithmetic (no 4x4 matrix
  // products) and retains terms far below the spacing around one.
  const Numeric local_B = u * u + v * v + w * w - b * b - c * c - d * d;
  const Numeric local_C = -Math::pow2(d * u - c * v + b * w);
  const Numeric n2      = n * n;
  const Numeric n3      = n2 * n;

  Vector4 term{1.0, 0.0, 0.0, 0.0};
  Vector4 sum_m1{};
  for (Size k = 0; k < 80; ++k) {
    const Numeric inv_denominator = 1.0 / static_cast<Numeric>(k + 2);
    const auto    old             = term;
    term                          = {inv_denominator * (a * old[0] - local_C * old[3]),
                                     inv_denominator * (old[0] + a * old[1]),
                                     inv_denominator * (old[1] + a * old[2] - local_B * old[3]),
                                     inv_denominator * (old[2] + a * old[3])};
    for (Size j = 0; j < 4; ++j) sum_m1[j] += term[j];

    const Numeric term_bound =
        std::abs(term[0]) + n * std::abs(term[1]) + n2 * std::abs(term[2]) + n3 * std::abs(term[3]);
    const Numeric sum_bound =
        std::abs(sum_m1[0]) + n * std::abs(sum_m1[1]) + n2 * std::abs(sum_m1[2]) + n3 * std::abs(sum_m1[3]);
    if (static_cast<Numeric>(k + 2) > std::abs(a) + n and
        term_bound <= std::numeric_limits<Numeric>::epsilon() * (1.0 + sum_bound))
      break;
  }

  const Numeric c2_diag = sum_m1[2];
  return {sum_m1[0] + c2_diag * (b2 + c2 + d2),
          sum_m1[0] + c2_diag * (b2 - u2 - v2),
          sum_m1[0] + c2_diag * (c2 - u2 - w2),
          sum_m1[0] + c2_diag * (d2 - v2 - w2)};
}

muelmat tran::linsrc_deriv(const propmat &dk, const Numeric r, const Numeric dr) const {
  ARTS_USER_ERROR_IF(r == 0.0 and dr != 0.0,
                     "A nonzero layer-length derivative at r=0 requires the endpoint propagation matrices.");
  const Numeric dr_r = r == 0.0 ? 0.0 : dr / r;
  return linsrc_deriv(dr_r * propmat{a, b, c, d, u, v, w} - 0.5 * r * dk);
}

muelmat tran::linsrc_deriv(const propmat &dg) const {
  const auto &[da, db, dc, dd, du, dv, dw] = dg.data;

  if (not polarized) return phi1_scalar_derivative_from_exp(a, exp_a) * to_muelmat(dg);

  const Numeric phi1_a = a == 0.0 ? 1.0 : expm1_a / a;
  const propmat g{a, b, c, d, u, v, w};
  const Numeric n              = polarized_row_norm(g);
  const Numeric generator_norm = std::abs(a) + n;
  const Numeric polarization_scale =
      std::max({std::abs(b), std::abs(c), std::abs(d), std::abs(u), std::abs(v), std::abs(w)});
  const bool phi1_ill_conditioned = x_zero or y_zero or spectrally_clustered or
                                    (S > 0.0 and S <= 0.01 * polarization_scale * polarization_scale) or
                                    std::abs(a + x) < 0.01 or std::abs(a - x) < 0.01 or a * a + y * y < 1e-4;
  const bool separated_real_limit = use_real_pair_frechet_limit(*this) and
                                    S > 0.01 * polarization_scale * polarization_scale and std::isfinite(S) and
                                    std::isfinite(x2) and std::isfinite(y2);

  if (generator_norm <= 2.0 or (phi1_ill_conditioned and generator_norm <= 8.0))
    return phi1_frechet_cayley_hamilton_series(g, dg);

  if (separated_real_limit) return separated_real_frechet(*this, dg, frechet_function::phi1);

  if ((not stable_coefficients or phi1_ill_conditioned or exp_a < std::numeric_limits<Numeric>::min() or x > 500.0) and
      not separated_real_limit) {
    return phi1_frechet_degenerate(g, dg);
  }

  const Numeric phi1p_a = (x_zero or y_zero) ? phi1_scalar_derivative_from_exp(a, exp_a) : 0.0;
  const Numeric db2     = 2.0 * db * b;
  const Numeric dc2     = 2.0 * dc * c;
  const Numeric dd2     = 2.0 * dd * d;
  const Numeric du2     = 2.0 * du * u;
  const Numeric dv2     = 2.0 * dv * v;
  const Numeric dw2     = 2.0 * dw * w;
  const Numeric dB      = du2 + dv2 + dw2 - db2 - dc2 - dd2;
  const Numeric dC      = -2.0 * (d * u - c * v + b * w) * (dd * u + d * du - dc * v - c * dv + db * w + b * dw);
  const Numeric dS_val  = (S > 1e-9) ? (B * dB - 2.0 * dC) / S : 0.0;
  const Numeric dx2     = 0.5 * (dS_val - dB);
  const Numeric dy2     = 0.5 * (dS_val + dB);
  const Numeric dx      = (x > 1e-9) ? 0.5 * dx2 / x : 0.0;
  const Numeric dy      = (y > 1e-9) ? 0.5 * dy2 / y : 0.0;

  Numeric l1, l2, l3;
  Numeric dl0, dl1, dl2, dl3;

  if (both_zero) {
    const Numeric fpa  = phi1p_a;
    const Numeric fppa = phi1_scalar_derivative<2>(a);
    dl0                = fpa * da;
    l1                 = fpa;
    dl1                = fppa * da;

    const Numeric f3pa = phi1_scalar_derivative<3>(a);
    const Numeric f4pa = phi1_scalar_derivative<4>(a);
    l2                 = 0.5 * fppa;
    dl2                = 0.5 * f3pa * da;
    l3                 = f3pa / 6.0;
    dl3                = f4pa / 6.0 * da;
  } else {
    Numeric Pp = 0.0, Pm_div_x = 0.0, Qp = 0.0, q_im = 0.0;
    Numeric dPp = 0.0, dPm_div_x = 0.0, dQp = 0.0, dq_im = 0.0;

    if (x_zero) {
      Pp  = phi1_a;
      dPp = phi1p_a * da + 0.5 * phi1_scalar_derivative<2>(a) * dx2;

      Pm_div_x  = phi1p_a;
      dPm_div_x = phi1_scalar_derivative<2>(a) * da + (phi1_scalar_derivative<3>(a) / 6.0) * dx2;

    } else {
      const Numeric f_apx   = phi1_scalar(a + x);
      const Numeric f_amx   = phi1_scalar(a - x);
      const Numeric fp_apx  = phi1_scalar_derivative<1>(a + x);
      const Numeric fp_amx  = phi1_scalar_derivative<1>(a - x);
      Pp                    = 0.5 * (f_apx + f_amx);
      const Numeric sum_fp  = fp_apx + fp_amx;
      const Numeric diff_fp = fp_apx - fp_amx;
      dPp                   = 0.5 * sum_fp * da + 0.5 * diff_fp * dx;
      Pm_div_x              = 0.5 * (f_apx - f_amx) / x;
      const Numeric dPm_da  = 0.5 * diff_fp / x;

      Numeric dPm_dx;
      if (x < 1e-3) {
        dPm_dx = phi1_scalar_derivative<3>(a) * x / 3.0;
      } else {
        const Numeric diff_f = f_apx - f_amx;
        dPm_dx               = (x * sum_fp - diff_f) / (2.0 * x * x);
      }

      dPm_div_x = dPm_da * da + dPm_dx * dx;
    }

    if (y_zero) {
      Qp    = phi1_a;
      dQp   = phi1p_a * da - 0.5 * phi1_scalar_derivative<2>(a) * dy2;
      q_im  = phi1p_a;
      dq_im = phi1_scalar_derivative<2>(a) * da - (phi1_scalar_derivative<3>(a) / 6.0) * dy2;
    } else {
      const Numeric denom       = a * a + y * y;
      const Numeric ea_cy_m1    = exp_a * cy - 1.0;
      const Numeric ea_sy       = exp_a * sy;
      Qp                        = (a * ea_cy_m1 + y * ea_sy) / denom;
      const Numeric sin_y_div_y = (std::abs(y) < 1e-6) ? 1.0 - y * y / 6.0 : sy / y;
      q_im                      = (a * exp_a * sin_y_div_y - ea_cy_m1) / denom;
      const Numeric ImF         = (a * exp_a * sy - y * ea_cy_m1) / denom;
      const Numeric A_val       = exp_a * ((a - 1.0) * cy - y * sy) + 1.0;
      const Numeric B_val       = exp_a * ((a - 1.0) * sy + y * cy);
      const Numeric C_val       = a * a - y * y;
      const Numeric D_val       = 2.0 * a * y;
      const Numeric denom2      = C_val * C_val + D_val * D_val;
      const Numeric Fp_re       = (A_val * C_val + B_val * D_val) / denom2;
      const Numeric Fp_im       = (B_val * C_val - A_val * D_val) / denom2;
      dQp                       = Fp_re * da - Fp_im * dy;
      const Numeric dImF        = Fp_im * da + Fp_re * dy;
      dq_im                     = (y * dImF - ImF * dy) / (y * y);
    }

    const Numeric inv  = inv_x2y2;
    const Numeric dinv = -inv * inv * (dx2 + dy2);
    l2                 = (Pp - Qp) * inv;
    dl2                = (dPp - dQp) * inv + (Pp - Qp) * dinv;
    dl0                = dPp - dl2 * x2 - l2 * dx2;
    l3                 = (Pm_div_x - q_im) * inv;
    dl3                = (dPm_div_x - dq_im) * inv + (Pm_div_x - q_im) * dinv;
    l1                 = Pm_div_x - l3 * x2;
    dl1                = dPm_div_x - dl3 * x2 - l3 * dx2;
  }

  const muelmat S_mat{0, b, c, d, b, 0, u, v, c, -u, 0, w, d, -v, -w, 0};
  const muelmat dS_mat{0, db, dc, dd, db, 0, du, dv, dc, -du, 0, dw, dd, -dv, -dw, 0};

  const muelmat S2_mat  = S_mat * S_mat;
  const muelmat dS2_mat = dS_mat * S_mat + S_mat * dS_mat;
  const muelmat S3_mat  = S_mat * S2_mat;
  const muelmat dS3_mat = dS_mat * S2_mat + S_mat * dS2_mat;

  return muelmat{dl0, 0, 0, 0, 0, dl0, 0, 0, 0, 0, dl0, 0, 0, 0, 0, dl0} + S_mat * dl1 + dS_mat * l1 + S2_mat * dl2 +
         dS2_mat * l2 + S3_mat * dl3 + dS3_mat * l3;
}

muelmat tran::linsrc_linprop(const muelmat &t, const propmat &k1, const propmat &k2, const Numeric r) const noexcept {
  static_cast<void>(t);

  const scalar_linprop_result state = scalar_linprop(r * k1.A(), r * k2.A(), a, exp_a, expm1_a);
  if (not polarized) return {state.value};

  // The scalar completion-of-the-square formula does not extend by applying
  // Dawson element by element.  Treat the isotropic gradient with its exact
  // scalar integral and the polarized gradient with the commutator-free
  // augmented-source term.  This also makes the result continuous as the
  // polarization tends to zero.
  return augmented_linprop_value(*this, state, k1, k2, r);
}

muelmat tran::linsrc_linprop_deriv(const muelmat &lambda,
                                   const muelmat &t,
                                   const propmat &k1,
                                   const propmat &k2,
                                   const propmat &dk_in,
                                   const muelmat &dt,
                                   const Numeric  r,
                                   const Numeric  dr,
                                   bool           k1_deriv) const {
  const scalar_linprop_result state  = scalar_linprop(r * k1.A(), r * k2.A(), a, exp_a, expm1_a);
  const Numeric               dvalue = scalar_linprop_deriv(state, k1.A(), k2.A(), dk_in.A(), r, dr, k1_deriv);

  if (not polarized and not dk_in.is_polarized()) {
    static_cast<void>(lambda);
    static_cast<void>(t);
    static_cast<void>(dt);
    return {dvalue};
  }

  static_cast<void>(lambda);
  static_cast<void>(t);
  static_cast<void>(dt);
  const muelmat f         = linsrc();
  const muelmat q_m1      = -(1.0 / 12.0) * to_muelmat(r * k2 - r * k1);
  const propmat kavg      = avg(k1, k2);
  const Numeric phi       = a == 0.0 ? 1.0 : expm1_a / a;
  const Numeric phi_prime = phi1_scalar_derivative_from_exp(a, exp_a);
  return augmented_linprop_deriv(*this, state, f, q_m1, kavg, phi, phi_prime, k1, k2, dk_in, r, dr, k1_deriv);
}

muelmat tran::deriv(const muelmat &t,
                    const propmat &k1,
                    const propmat &k2,
                    const propmat &dk,
                    const Numeric  r,
                    const Numeric  dr) const {
  return deriv(t, -0.5 * r * dk - dr * avg(k1, k2));
}

muelmat tran::deriv(const muelmat &t, const propmat &dg) const {
  const auto &[da, db, dc, dd, du, dv, dw] = dg.data;

  if (not polarized) return exp_a * to_muelmat(dg);

  if (x > 2.0 and std::isfinite(S) and use_real_pair_frechet_limit(*this))
    return separated_real_frechet(*this, dg, frechet_function::exponential);

  if (not stable_coefficients or spectrally_clustered or x_zero or y_zero or
      exp_a < std::numeric_limits<Numeric>::min() or x > 500.0) {
    return exp_frechet_degenerate(propmat{a, b, c, d, u, v, w}, dg);
  }

  const Numeric db2 = 2 * db * b;
  const Numeric dc2 = 2 * dc * c;
  const Numeric dd2 = 2 * dd * d;
  const Numeric du2 = 2 * du * u;
  const Numeric dv2 = 2 * dv * v;
  const Numeric dw2 = 2 * dw * w;

  /* Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
  const Numeric dB = du2 + dv2 + dw2 - db2 - dc2 - dd2;
  const Numeric dC = -2 * (d * u - c * v + b * w) * (dd * u + d * du - dc * v - c * dv + db * w + b * dw);
  const Numeric dS = (B * dB - 2 * dC) / S;

  const Numeric dx2 = 0.5 * (dS - dB);
  const Numeric dy2 = 0.5 * (dS + dB);
  Numeric       dC0, dC1, dC2, dC3;

  if (S < 1e-4) {
    // Differentiate the same coefficient polynomials used by
    // init_polarized().  Re-entering the trigonometric divided differences
    // here would reintroduce cancellation into an otherwise stable value.
    const Numeric x4     = x2 * x2;
    const Numeric y4     = y2 * y2;
    const Numeric dx4    = 2.0 * x2 * dx2;
    const Numeric dy4    = 2.0 * y2 * dy2;
    const Numeric p4     = x4 - x2 * y2 + y4;
    const Numeric dp4    = dx4 - dx2 * y2 - x2 * dy2 + dy4;
    const Numeric dp6    = dx4 * x2 + x4 * dx2 - dx4 * y2 - x4 * dy2 + dx2 * y4 + x2 * dy4 - dy4 * y2 - y4 * dy2;
    const Numeric xy2    = x2 * y2;
    const Numeric dxy2   = dx2 * y2 + x2 * dy2;
    const Numeric delta  = x2 - y2;
    const Numeric ddelta = dx2 - dy2;

    dC0 = dxy2 * (1.0 / 24.0 + delta / 720.0 + p4 / 40320.0) + xy2 * (ddelta / 720.0 + dp4 / 40320.0);
    dC1 = dxy2 / 120.0 + (dxy2 * delta + xy2 * ddelta) / 5040.0;
    dC2 = ddelta / 24.0 + dp4 / 720.0 + dp6 / 40320.0;
    dC3 = ddelta / 120.0 + dp4 / 5040.0 + dp6 / 362880.0;
  } else {
    const Numeric dx     = 0.5 * dx2 / x;
    const Numeric dy     = 0.5 * dy2 / y;
    const Numeric dcy    = -sy * dy;
    const Numeric dsy    = cy * dy;
    const Numeric dcx    = sx * dx;
    const Numeric dsx    = cx * dx;
    const Numeric dix    = -dx * ix * ix;
    const Numeric diy    = -dy * iy * iy;
    const Numeric dx2dy2 = dx2 + dy2;

    dC0 = (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0 * dx2dy2) * inv_x2y2;
    dC1 =
        (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy + dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix - C1 * dx2dy2) *
        inv_x2y2;
    dC2 = ((dcx - C2 * dx2) - (dcy + C2 * dy2)) * inv_x2y2;
    dC3 = ((dsx * ix + sx * dix - C3 * dx2) - (dsy * iy + sy * diy + C3 * dy2)) * inv_x2y2;
  }

  const Numeric dC2b = dC2 * (c * u + d * v) + C2 * (dc * u + c * du + dd * v + d * dv);
  const Numeric dC2c = dC2 * (b * u - d * w) + C2 * (db * u + b * du - dd * w - d * dw);
  const Numeric dC2d = dC2 * (b * v + c * w) + C2 * (db * v + b * dv + dc * w + c * dw);
  const Numeric dC2u = dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw);
  const Numeric dC2v = dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw);
  const Numeric dC2w = dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv);

  const Numeric dC3b =
      dC3 * (b * (B - w2) + w * (c * v - d * u)) +
      C3 * (db * (B - w2) + b * (dB - dw2) + dw * (c * v - d * u) + w * (dc * v + c * dv - dd * u - d * du));
  const Numeric dC3c =
      dC3 * (c * (v2 - B) - v * (d * u + b * w)) +
      C3 * (dc * (v2 - B) + c * (dv2 - dB) - dv * (d * u + b * w) - v * (dd * u + d * du + db * w + b * dw));
  const Numeric dC3d =
      dC3 * (d * (u2 - B) - u * (c * v - b * w)) +
      C3 * (dd * (u2 - B) + d * (du2 - dB) - du * (c * v - b * w) - u * (dc * v + c * dv - db * w - b * dw));
  const Numeric dC3u =
      dC3 * (d * (c * v - b * w) - u * (B + d2)) +
      C3 * (dd * (c * v - b * w) + d * (dc * v + c * dv - db * w - b * dw) - du * (B + d2) - u * (dB + dd2));
  const Numeric dC3v =
      dC3 * (c * (d * u + b * w) - v * (B + c2)) +
      C3 * (dc * (d * u + b * w) + c * (dd * u + d * du + db * w + b * dw) - dv * (B + c2) - v * (dB + dc2));
  const Numeric dC3w =
      dC3 * (b * (c * v - d * u) - w * (B + b2)) +
      C3 * (db * (c * v - d * u) + b * (dc * v + c * dv - dd * u - d * du) - dw * (B + b2) - w * (dB + db2));

  const Numeric dM00 = dC0 + dC2 * (b2 + c2 + d2) + C2 * (db2 + dc2 + dd2);
  const Numeric dM11 = dC0 + dC2 * (b2 - u2 - v2) + C2 * (db2 - du2 - dv2);
  const Numeric dM22 = dC0 + dC2 * (c2 - u2 - w2) + C2 * (dc2 - du2 - dw2);
  const Numeric dM33 = dC0 + dC2 * (d2 - v2 - w2) + C2 * (dd2 - dv2 - dw2);

  // clang-format off
  return da * t + exp_a * muelmat{
    dM00,                              dC1 * b + C1 * db - dC2b - dC3b,   dC1 * c + C1 * dc + dC2c + dC3c, dC1 * d + C1 * dd + dC2d + dC3d,
    dC1 * b + C1 * db + dC2b - dC3b,   dM11,                              dC1 * u + C1 * du + dC2u + dC3u, dC1 * v + C1 * dv + dC2v + dC3v,
    dC1 * c + C1 * dc - dC2c + dC3c, - dC1 * u - C1 * du + dC2u - dC3u,   dM22,                            dC1 * w + C1 * dw + dC2w + dC3w,
    dC1 * d + C1 * dd - dC2d + dC3d, - dC1 * v - C1 * dv + dC2v - dC3v, - dC1 * w - C1 * dw + dC2w - dC3w, dM33};
  // clang-format on
}

muelmat tran::magnus_deriv(const propmat &k1,
                           const propmat &k2,
                           const propmat &dk,
                           const Numeric  r,
                           const Numeric  dr,
                           const bool     k1_deriv) const {
  return magnus_deriv((*this)(), k1, k2, dk, r, dr, k1_deriv);
}

muelmat tran::magnus_deriv(const muelmat &t,
                           const propmat &k1,
                           const propmat &k2,
                           const propmat &dk,
                           const Numeric  r,
                           const Numeric  dr,
                           const bool     k1_deriv) const {
  return deriv(t, magnus_exponent_deriv(k1, k2, dk, r, dr, k1_deriv));
}

muelmat tran::magnus_linsrc(const propmat &k1, const propmat &k2, const Numeric r) const {
  return augmented_magnus_source(*this, k1, k2, r);
}

muelmat tran::magnus_linsrc_deriv(const propmat &k1,
                                  const propmat &k2,
                                  const propmat &dk,
                                  const Numeric  r,
                                  const Numeric  dr,
                                  const bool     k1_deriv) const {
  const muelmat f    = linsrc();
  const muelmat q_m1 = -(1.0 / 12.0) * to_muelmat(r * k2 - r * k1);
  return augmented_magnus_source_deriv(*this, f, q_m1, k1, k2, dk, r, dr, k1_deriv);
}

muelmat exp(propmat k, Numeric r) { return tran(k, k, r)(); }

propmat logK(const muelmat &m) {
  if (not m.is_polarized()) return std::log(midtr(m));

  /**
    The code tries to retrieve K from exp(K) = M input as muelmat m,
    where K is a propmat.

    K is unknown but, as a reminder, it looks like:
      K = [ a    b    c    d
            b    a    u    v
            c   -u    a    w
            d   -v   -w    a ],

    where a, b, c, d, u, v, w are real numbers.

    Since det(e^A) = e^tr(A), we can find a from the determinant of m,
    using tr(K) = 4 a.
   */

  // Small fix: in cases this matter, logK is not useful
  const Numeric det_m = det(m);
  const Numeric a     = (det_m > std::numeric_limits<Numeric>::min()) ? 0.25 * std::log(det_m) : 0.0;

  /**
    We use this knowing
        m = exp(K)
          = exp(a) * (C0 I + C1 K + C2 K^2 + C3 K^3)

    where C0, C1, C2, C3 are the Cayley Hamilton coefficients derived in tran::tran.
    
    This means that we can extract
        T = m * exp(-a)
          = C0 I + C1 K + C2 K^2 + C3 K^3
    using the same logic and notation.
   */

  const muelmat T = m * std::exp(-a);

  /**
    These coefficients C0, C1, C2, C3 are found from the eigenvalues of K, which
    have the form: lambda = exp(x), exp(-x), exp(iy), exp(-iy), where x and y are real.

    Since the trace of a matrix is the sum of its eigenvalues, and the eigenvalues
    of a squared matrix are the squares of the eigenvalues of the original matrix,
    we can form these two relations from our known T matrix:
      S1 = tr(T)
         = e^x + e^-x + e^iy + e^-iy
         = 2 cosh(x) + 2 cos(y)
    and
      S2 = tr(T^2)
         = (e^x + e^-x)^2 + (e^iy + e^-iy)^2
         = 2 cosh(2x) + 2 cos(2y)
         = 2 (2 cosh^2(x) - 1) + 2 (2 cos^2(y) - 1)

    Letting u = cosh(x), v = cos(y), we can write:
        S1 = 2 u + 2 v  ->  
            u + v = S1 / 2
    and
        S2 = 2 (2 u^2 - 1) + 2 (2 v^2 - 1)
           = 4 u^2 + 4 v^2 - 4 ->
            u^2 + v^2 = (S2 + 4) / 4

    Since
        (u - v)^2 = 2 (u^2 + v^2) - (u + v)^2
    we can solve for u and v:
        (u - v)^2 = 2 * (S2 + 4) / 4 - (S1 / 2)^2
                  = (S2 + 8 - S1^2) / 4
    and write
        D = (u - v)^2
          = (S2 + 8 - S1^2) / 4
   */

  const Numeric S1 = tr(T);
  const Numeric D  = 0.25 * (2 * tr(T * T) + 8 - S1 * S1);

  /**
    We can use D to find u and v:
        u = (u + v + u - v) / 2
          = (S1 / 2 + sqrt(D)) / 2
    and
        v = (u + v - (u - v)) / 2
          = (S1 / 2 - sqrt(D)) / 2
    Since u - v = sqrt(D)
   */

  const Numeric sqrtD = D > 0 ? std::sqrt(D) : 0.0;
  const Numeric u     = std::max((S1 * 0.5 + sqrtD) * 0.5, 1.0);
  const Numeric v     = std::clamp((S1 * 0.5 - sqrtD) * 0.5, -1.0, 1.0);

  /**
    From the definition of u and v we can now define x and y:
        x = acosh(u)
        y = acos(v)
   */

  const Numeric x = std::acosh(u);
  const Numeric y = std::acos(v);

  // From here, we simply compute the coefficients C0, C1, C2, C3 as in tran::tran
  const Numeric x2          = x * x;
  const Numeric y2          = y * y;
  const bool    x_zero      = x < too_small;
  const bool    y_zero      = y < too_small;
  const bool    both_zero   = x_zero && y_zero;
  const bool    either_zero = x_zero || y_zero;
  const Numeric cy          = v;
  const Numeric sy          = std::sin(y);
  const Numeric cx          = u;
  const Numeric sx          = std::sinh(x);
  const Numeric ix          = x_zero ? 0.0 : 1.0 / x;
  const Numeric iy          = y_zero ? 0.0 : 1.0 / y;
  const Numeric inv_x2y2    = both_zero ? 1.0 : 1.0 / (x2 + y2);
  const Numeric C0          = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
  const Numeric C1          = either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
  // Skipping C2 and C3

  /**
    Now we are just one step away from reconstructing log(T) and thus K from exp(K).

    We have:
      T = C0 I + C1 K + C2 K^2 + C3 K^3

    We write:
        T02 = C0 I + C2 K^2  ->
                K^2 = (T02 - C0 I) / C2
    and
        T13 = C1 K + C3 K^3  ->
                T13 = K (C1 I + C3 K^2)

    So that:
        K = T13 * inv(C1 I + C3 K^2)
          = T13 * inv(C1 I + C3 (T02 - C0 I) / C2)
    
    OK.  We just have to be able to construct T02 and T13 from T.
    This requires some separation of terms in T.

    For T02, 
      [   T[0,0]             (T[0,1]-T[1,0])/2   (T[0,2]-T[2,0])/2   (T[0,3]-T[3,0])/2
        -(T[1,0]-T[0,1])/2    T[1,1]             (T[1,2]+T[2,1])/2   (T[1,3]+T[3,1])/2
        -(T[2,0]-T[0,2])/2   (T[2,1]+T[1,2])/2    T[2,2]             (T[2,3]+T[3,2])/2
        -(T[3,0]-T[0,3])/2   (T[3,1]+T[1,3])/2   (T[3,2]+T[2,3])/2    T[3,3] ]
    For T13,
      [  0                   (T[0,1]+T[1,0])/2    (T[0,2]+T[2,0])/2   (T[0,3]+T[3,0])/2
        (T[1,0]+T[0,1])/2     0                   (T[1,2]-T[2,1])/2   (T[1,3]-T[3,1])/2
        (T[2,0]+T[0,2])/2   -(T[2,1]-T[1,2])/2     0                  (T[2,3]-T[3,2])/2
        (T[3,0]+T[0,3])/2   -(T[3,1]-T[1,3])/2   -(T[3,2]-T[2,3])/2    0 ]
    
    (To make the above symmetries and anti-symmetries clear, the return value of
     tran::operator() has been formatted to demonstrate them.)
   */
  const Numeric factor =
      (both_zero or sqrtD <= 0) ? 1.0 / 3.0 : ((x_zero ? 1.0 : sx * ix) - (y_zero ? 1.0 : sy * iy)) / sqrtD;
  const Numeric z01 = factor * std::midpoint(T[0, 1], -T[1, 0]);
  const Numeric z02 = factor * std::midpoint(T[0, 2], -T[2, 0]);
  const Numeric z03 = factor * std::midpoint(T[0, 3], -T[3, 0]);
  const Numeric z12 = factor * std::midpoint(T[1, 2], T[2, 1]);
  const Numeric z13 = factor * std::midpoint(T[1, 3], T[3, 1]);
  const Numeric z23 = factor * std::midpoint(T[2, 3], T[3, 2]);
  const Numeric z00 = factor * (T[0, 0] - C0) + C1;
  const Numeric z11 = factor * (T[1, 1] - C0) + C1;
  const Numeric z22 = factor * (T[2, 2] - C0) + C1;
  const Numeric z33 = factor * (T[3, 3] - C0) + C1;
  const Numeric b01 = std::midpoint(T[0, 1], T[1, 0]);
  const Numeric b02 = std::midpoint(T[0, 2], T[2, 0]);
  const Numeric b03 = std::midpoint(T[0, 3], T[3, 0]);
  const Numeric b12 = std::midpoint(T[1, 2], -T[2, 1]);
  const Numeric b13 = std::midpoint(T[1, 3], -T[3, 1]);
  const Numeric b23 = std::midpoint(T[2, 3], -T[3, 2]);

  // clang-format off
  const muelmat Kp = muelmat{ 0.0,  b01,  b02, b03,
                              b01,  0.0,  b12, b13,
                              b02, -b12,  0.0, b23,
                              b03, -b13, -b23, 0.0} *
              inv(muelmat{ z00,  z01,  z02, z03,
                             -z01,  z11,  z12, z13,
                             -z02,  z12,  z22, z23,
                             -z03,  z13,  z23, z33});
  // clang-format on

  return {a,
          std::midpoint(Kp[0, 1], Kp[1, 0]),
          std::midpoint(Kp[0, 2], Kp[2, 0]),
          std::midpoint(Kp[0, 3], Kp[3, 0]),
          std::midpoint(Kp[1, 2], -Kp[2, 1]),
          std::midpoint(Kp[1, 3], -Kp[3, 1]),
          std::midpoint(Kp[2, 3], -Kp[3, 2])};
}

specmat sqrt(const propmat &pm) {
  const Numeric a      = pm.A();
  const Complex sqrt_a = std::sqrt(Complex(a));

  if (not pm.is_polarized()) return sqrt_a;

  const Numeric b = pm.B();
  const Numeric c = pm.C();
  const Numeric d = pm.D();
  const Numeric u = pm.U();
  const Numeric v = pm.V();
  const Numeric w = pm.W();

  const Numeric b2 = b * b;
  const Numeric c2 = c * c;
  const Numeric d2 = d * d;
  const Numeric u2 = u * u;
  const Numeric v2 = v * v;
  const Numeric w2 = w * w;

  Complex d0c{}, d1c{}, d2c{}, d3c{};

  if (pm.is_rotational()) {
    const Numeric rho = std::hypot(u, v, w);
    if (rho == 0.0) return {0.0};

    // sqrt(N) = sqrt(rho/2) (N/rho - (N/rho)^2).  Normalizing first
    // preserves every nonzero binary64 rotation and avoids the singular
    // 1/(rho*sqrt(rho)) coefficient of the equivalent unscaled formula.
    const Numeric nh_u = u / rho;
    const Numeric nh_v = v / rho;
    const Numeric nh_w = w / rho;
    const Numeric s    = Constant::inv_sqrt_2 * std::sqrt(rho);
    return {0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            s * (nh_u * nh_u + nh_v * nh_v),
            s * (nh_u + nh_v * nh_w),
            s * (nh_v - nh_u * nh_w),
            0.0,
            s * (-nh_u + nh_v * nh_w),
            s * (nh_u * nh_u + nh_w * nh_w),
            s * (nh_w + nh_u * nh_v),
            0.0,
            s * (-nh_v - nh_u * nh_w),
            s * (-nh_w + nh_u * nh_v),
            s * (nh_v * nh_v + nh_w * nh_w)};
  } else {
    /** Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
    const Numeric B = u2 + v2 + w2 - b2 - c2 - d2;
    const Numeric q = d * u - c * v + b * w;
    const Numeric S = std::hypot(B, 2.0 * q);

    Numeric x2, abs_y2;
    if (B >= 0.0) {
      abs_y2 = 0.5 * (S + B);
      x2     = abs_y2 == 0.0 ? 0.0 : (std::abs(q) / abs_y2) * std::abs(q);
    } else {
      x2     = 0.5 * (S - B);
      abs_y2 = x2 == 0.0 ? 0.0 : (std::abs(q) / x2) * std::abs(q);
    }
    const Numeric x  = std::sqrt(x2);
    const Complex y  = Complex(0, std::sqrt(abs_y2));
    const Complex sx = std::sqrt(Complex(a + x));
    const Complex dx = std::sqrt(Complex(a - x));
    const Complex sy = std::sqrt(a + y);
    const Complex dy = std::sqrt(a - y);
    const Complex Sx = sx + dx;
    const Complex Sy = sy + dy;

    if (S == 0.0) {
      ARTS_USER_ERROR_IF(a == 0.0,
                         "The principal square root is not defined for a nonzero nilpotent propagation matrix.");

      // sqrt(a I + N), with N^4=0.  This exact finite Taylor polynomial
      // handles the coalesced-eigenvalue case without divided differences.
      d0c = sqrt_a;
      d1c = 0.5 / sqrt_a;
      d2c = -0.125 / (a * sqrt_a);
      d3c = 0.0625 / (a * a * sqrt_a);
    } else {
      const Numeric inv_S = 1.0 / S;

      // Rationalized divided differences.  1/Sx and 1/Sy are exactly the
      // odd first divided differences, but avoid subtracting nearby roots.
      const Complex term_x = Sx == 0.0 ? 0.0 : 1.0 / Sx;
      const Complex term_y = Sy == 0.0 ? 0.0 : 1.0 / Sy;
      d0c                  = (abs_y2 * Sx + x2 * Sy) * (0.5 * inv_S);
      d1c                  = (abs_y2 * term_x + x2 * term_y) * inv_S;

      // Rationalize Sx-Sy twice.  This remains conditioned as x,y -> 0:
      // (Sx-Sy)/(2S) = -1 / ((sx*dx+sy*dy)(Sx+Sy)).
      const Complex common = (sx * dx + sy * dy) * (Sx + Sy);
      d2c                  = common == 0.0 ? 0.5 * (Sx - Sy) * inv_S : -1.0 / common;
      d3c                  = Sx == 0.0 or Sy == 0.0 ? (term_x - term_y) * inv_S : -2.0 * d2c / (Sx * Sy);
    }
  }

  // Terms for K0^2
  const Numeric k2_00 = b2 + c2 + d2;
  const Numeric k2_11 = b2 - u2 - v2;
  const Numeric k2_22 = c2 - u2 - w2;
  const Numeric k2_33 = d2 - v2 - w2;

  const Numeric k2_01_val = -(c * u + d * v);
  const Numeric k2_02_val = b * u - d * w;
  const Numeric k2_03_val = b * v + c * w;
  const Numeric k2_12_val = b * c - v * w;
  const Numeric k2_13_val = b * d + u * w;
  const Numeric k2_23_val = c * d - u * v;

  // Terms for K0^3
  const Numeric k3_01_val = b * k2_00 - u * k2_02_val - v * k2_03_val;
  const Numeric k3_02_val = c * k2_00 + u * k2_01_val - w * k2_03_val;
  const Numeric k3_03_val = d * k2_00 + v * k2_01_val + w * k2_02_val;
  const Numeric k3_12_val = -c * k2_01_val + u * k2_11 - w * k2_13_val;
  const Numeric k3_13_val = -d * k2_01_val + v * k2_11 + w * k2_12_val;
  const Numeric k3_23_val = -d * k2_02_val + v * k2_12_val + w * k2_22;

  specmat K;

  K[0, 0] = d0c + d2c * k2_00;
  K[1, 1] = d0c + d2c * k2_11;
  K[2, 2] = d0c + d2c * k2_22;
  K[3, 3] = d0c + d2c * k2_33;

  K[0, 1] = d1c * b + d2c * k2_01_val + d3c * k3_01_val;
  K[1, 0] = d1c * b - d2c * k2_01_val + d3c * k3_01_val;

  K[0, 2] = d1c * c + d2c * k2_02_val + d3c * k3_02_val;
  K[2, 0] = d1c * c - d2c * k2_02_val + d3c * k3_02_val;

  K[0, 3] = d1c * d + d2c * k2_03_val + d3c * k3_03_val;
  K[3, 0] = d1c * d - d2c * k2_03_val + d3c * k3_03_val;

  K[1, 2] = d1c * u + d2c * k2_12_val + d3c * k3_12_val;
  K[2, 1] = -d1c * u + d2c * k2_12_val - d3c * k3_12_val;

  K[1, 3] = d1c * v + d2c * k2_13_val + d3c * k3_13_val;
  K[3, 1] = -d1c * v + d2c * k2_13_val - d3c * k3_13_val;

  K[2, 3] = d1c * w + d2c * k2_23_val + d3c * k3_23_val;
  K[3, 2] = -d1c * w + d2c * k2_23_val - d3c * k3_23_val;

  return K;
}

void TransmittanceMatrix::constant(const std::span<const propmat>        &K,
                                   const std::span<const propmat_vector> &dK,
                                   const ConstVectorView                 &r,
                                   const ConstTensor3View                &dr) {
  const Size N  = K.size();
  const Size nq = dT.ncols();

  auto dT0 = dT[0, 0];
  auto dT1 = dT[1, 0];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    const tran         tr{K[i - 1], K[i], r[i - 1]};
    const propmat      kavg = avg(K[i - 1], K[i]);
    diagonal_remainder t_diag_m1;
    T[0, i] = tr(t_diag_m1);
    set_diagonal_remainder(T_diag_m1, 0, i, t_diag_m1);

    for (Size iq = 0; iq < nq; iq++) {
      dT0[i - 1, iq] = tr.deriv(T[0, i], -0.5 * r[i - 1] * dK[i - 1][iq] - dr0[i - 1, iq] * kavg);
      dT1[i, iq]     = tr.deriv(T[0, i], -0.5 * r[i - 1] * dK[i][iq] - dr1[i - 1, iq] * kavg);
    }
  }
}

void TransmittanceMatrix::linsrc(const std::span<const propmat>        &K,
                                 const std::span<const propmat_vector> &dK,
                                 const ConstVectorView                 &r,
                                 const ConstTensor3View                &dr) {
  const Size N  = K.size();
  const Size nq = dT.ncols();

  auto dT0 = dT[0, 0];
  auto dT1 = dT[1, 0];
  auto dL0 = dL[0, 0];
  auto dL1 = dL[1, 0];
  auto dr0 = dr[0];
  auto dr1 = dr[1];
  auto T_  = T[0];
  auto L_  = L[0];

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    auto              &k1 = K[i - 1];
    auto              &k2 = K[i];
    const tran         tr{k1, k2, r[i - 1]};
    const propmat      kavg = avg(k1, k2);
    diagonal_remainder t_diag_m1;
    T_[i] = tr(t_diag_m1);
    set_diagonal_remainder(T_diag_m1, 0, i, t_diag_m1);
    diagonal_remainder l_diag_m1;
    L_[i] = tr.linsrc(l_diag_m1);
    set_diagonal_remainder(L_diag_m1, 0, i, l_diag_m1);

    for (Size iq = 0; iq < nq; iq++) {
      dT0[i - 1, iq] = tr.deriv(T_[i], -0.5 * r[i - 1] * dK[i - 1][iq] - dr0[i - 1, iq] * kavg);
      dT1[i, iq]     = tr.deriv(T_[i], -0.5 * r[i - 1] * dK[i][iq] - dr1[i - 1, iq] * kavg);
      dL0[i - 1, iq] = tr.linsrc_deriv(-0.5 * r[i - 1] * dK[i - 1][iq] - dr0[i - 1, iq] * kavg);
      dL1[i, iq]     = tr.linsrc_deriv(-0.5 * r[i - 1] * dK[i][iq] - dr1[i - 1, iq] * kavg);
    }
  }
}

void TransmittanceMatrix::linprop(const std::span<const propmat>        &K,
                                  const std::span<const propmat_vector> &dK,
                                  const ConstVectorView                 &r,
                                  const ConstTensor3View                &dr) {
  const Size N  = K.size();
  const Size nq = dT.ncols();

  auto dT0 = dT[0, 0];
  auto dT1 = dT[1, 0];
  auto dL0 = dL[0, 0];
  auto dL1 = dL[1, 0];
  auto dr0 = dr[0];
  auto dr1 = dr[1];
  auto T_  = T[0];
  auto L_  = L[0];

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    const propmat     &start = K[i];
    const propmat     &end   = K[i - 1];
    const Numeric      ri    = r[i - 1];
    const tran         tr{start, end, ri};
    diagonal_remainder t_diag_m1;
    T_[i] = tr(t_diag_m1);
    set_diagonal_remainder(T_diag_m1, 0, i, t_diag_m1);
    const scalar_linprop_result state = scalar_linprop(ri * start.A(), ri * end.A(), tr.a, tr.exp_a, tr.expm1_a);
    const propmat               kavg  = avg(start, end);

    diagonal_remainder     l_diag_m1;
    std::optional<muelmat> base_phi;
    std::optional<muelmat> base_q_m1;
    Numeric                base_phi_scalar;
    Numeric                base_phi_prime;
    if (not tr.polarized) {
      L_[i] = state.value;
      l_diag_m1.fill(state.value_m1);
    } else {
      L_[i] = augmented_linprop_value(
          tr, state, start, end, ri, &l_diag_m1, nq == 0 ? nullptr : &base_phi, nq == 0 ? nullptr : &base_q_m1);
      if (nq != 0) {
        base_phi_scalar = tr.a == 0.0 ? 1.0 : tr.expm1_a / tr.a;
        base_phi_prime  = phi1_scalar_derivative_from_exp(tr.a, tr.exp_a);
      }
    }
    set_diagonal_remainder(L_diag_m1, 0, i, l_diag_m1);

    for (Size iq = 0; iq < nq; iq++) {
      const Numeric  dr_end   = dr0[i - 1, iq];
      const Numeric  dr_start = dr1[i - 1, iq];
      const propmat &dk_end   = dK[i - 1][iq];
      const propmat &dk_start = dK[i][iq];
      dT0[i - 1, iq]          = tr.deriv(T_[i], -0.5 * ri * dk_end - dr_end * kavg);
      dT1[i, iq]              = tr.deriv(T_[i], -0.5 * ri * dk_start - dr_start * kavg);

      if (not tr.polarized and not dk_end.is_polarized()) {
        dL0[i - 1, iq] = muelmat{scalar_linprop_deriv(state, start.A(), end.A(), dk_end.A(), ri, dr_end, false)};
      } else {
        if (not base_phi.has_value()) {
          base_phi.emplace(tr.linsrc());
          base_q_m1.emplace(-(1.0 / 12.0) * to_muelmat(ri * end - ri * start));
          base_phi_scalar = tr.a == 0.0 ? 1.0 : tr.expm1_a / tr.a;
          base_phi_prime  = phi1_scalar_derivative_from_exp(tr.a, tr.exp_a);
        }
        dL0[i - 1, iq] = augmented_linprop_deriv(tr,
                                                 state,
                                                 *base_phi,
                                                 *base_q_m1,
                                                 kavg,
                                                 base_phi_scalar,
                                                 base_phi_prime,
                                                 start,
                                                 end,
                                                 dk_end,
                                                 ri,
                                                 dr_end,
                                                 false);
      }

      if (not tr.polarized and not dk_start.is_polarized()) {
        dL1[i, iq] = muelmat{scalar_linprop_deriv(state, start.A(), end.A(), dk_start.A(), ri, dr_start, true)};
      } else {
        if (not base_phi.has_value()) {
          base_phi.emplace(tr.linsrc());
          base_q_m1.emplace(-(1.0 / 12.0) * to_muelmat(ri * end - ri * start));
          base_phi_scalar = tr.a == 0.0 ? 1.0 : tr.expm1_a / tr.a;
          base_phi_prime  = phi1_scalar_derivative_from_exp(tr.a, tr.exp_a);
        }
        dL1[i, iq] = augmented_linprop_deriv(tr,
                                             state,
                                             *base_phi,
                                             *base_q_m1,
                                             kavg,
                                             base_phi_scalar,
                                             base_phi_prime,
                                             start,
                                             end,
                                             dk_start,
                                             ri,
                                             dr_start,
                                             true);
      }
    }
  }
}

void TransmittanceMatrix::magop(const std::span<const propmat>        &K,
                                const std::span<const propmat_vector> &dK,
                                const ConstVectorView                 &r,
                                const ConstTensor3View                &dr) {
  const Size N  = K.size();
  const Size nq = dT.ncols();

  auto dT0 = dT[0, 0];
  auto dT1 = dT[1, 0];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    const propmat     &start = K[i];
    const propmat     &end   = K[i - 1];
    const tran         tr{start, end, r[i - 1], MagnusOperator::magnus};
    diagonal_remainder t_diag_m1;
    T[0, i] = tr(t_diag_m1);
    set_diagonal_remainder(T_diag_m1, 0, i, t_diag_m1);

    for (Size iq = 0; iq < nq; iq++) {
      dT0[i - 1, iq] = tr.magnus_deriv(T[0, i], start, end, dK[i - 1][iq], r[i - 1], dr0[i - 1, iq], false);
      dT1[i, iq]     = tr.magnus_deriv(T[0, i], start, end, dK[i][iq], r[i - 1], dr1[i - 1, iq], true);
    }
  }
}

void TransmittanceMatrix::magop_linsrc(const std::span<const propmat>        &K,
                                       const std::span<const propmat_vector> &dK,
                                       const ConstVectorView                 &r,
                                       const ConstTensor3View                &dr) {
  const Size N  = K.size();
  const Size nq = dT.ncols();

  auto dT0 = dT[0, 0];
  auto dT1 = dT[1, 0];
  auto dL0 = dL[0, 0];
  auto dL1 = dL[1, 0];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; ++i) {
    const propmat     &start = K[i];
    const propmat     &end   = K[i - 1];
    const Numeric      ri    = r[i - 1];
    const tran         tr{start, end, ri, MagnusOperator::magnus};
    diagonal_remainder t_diag_m1;
    T[0, i] = tr(t_diag_m1);
    set_diagonal_remainder(T_diag_m1, 0, i, t_diag_m1);
    diagonal_remainder     l_diag_m1;
    std::optional<muelmat> base_phi;
    std::optional<muelmat> base_q_m1;
    L[0, i] = augmented_magnus_source(
        tr, start, end, ri, &l_diag_m1, nq == 0 ? nullptr : &base_phi, nq == 0 ? nullptr : &base_q_m1);
    set_diagonal_remainder(L_diag_m1, 0, i, l_diag_m1);

    for (Size iq = 0; iq < nq; ++iq) {
      dT0[i - 1, iq] = tr.magnus_deriv(T[0, i], start, end, dK[i - 1][iq], ri, dr0[i - 1, iq], false);
      dT1[i, iq]     = tr.magnus_deriv(T[0, i], start, end, dK[i][iq], ri, dr1[i - 1, iq], true);
      dL0[i - 1, iq] = augmented_magnus_source_deriv(
          tr, *base_phi, *base_q_m1, start, end, dK[i - 1][iq], ri, dr0[i - 1, iq], false);
      dL1[i, iq] =
          augmented_magnus_source_deriv(tr, *base_phi, *base_q_m1, start, end, dK[i][iq], ri, dr1[i - 1, iq], true);
    }
  }
}

void TransmittanceMatrix::constant(const std::span<const propmat_vector> &K,
                                   const std::span<const propmat_matrix> &dK,
                                   const ConstVectorView                 &r,
                                   const ConstTensor3View                &dr) {
  const Size nf = dT.npages();
  const Size np = dT.nrows();
  const Size nq = dT.ncols();

  auto dT0 = dT[0];
  auto dT1 = dT[1];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 1; i < np; i++) {
    for (Size iv = 0; iv < nf; ++iv) {
      const tran         tran_state{K[i - 1][iv], K[i][iv], r[i - 1]};
      const propmat      kavg = avg(K[i - 1][iv], K[i][iv]);
      diagonal_remainder t_diag_m1;
      T[iv, i] = tran_state(t_diag_m1);
      set_diagonal_remainder(T_diag_m1, iv, i, t_diag_m1);

      for (Size j = 0; j < nq; j++) {
        dT0[iv, i - 1, j] = tran_state.deriv(T[iv, i], -0.5 * r[i - 1] * dK[i - 1][j, iv] - dr0[i - 1, j] * kavg);
        dT1[iv, i, j]     = tran_state.deriv(T[iv, i], -0.5 * r[i - 1] * dK[i][j, iv] - dr1[i - 1, j] * kavg);
      }
    }
  }
}

void TransmittanceMatrix::linsrc(const std::span<const propmat_vector> &K,
                                 const std::span<const propmat_matrix> &dK,
                                 const ConstVectorView                 &r,
                                 const ConstTensor3View                &dr) {
  const Size nf = dT.npages();
  const Size np = dT.nrows();
  const Size nq = dT.ncols();

  auto dT0 = dT[0];
  auto dT1 = dT[1];
  auto dL0 = dL[0];
  auto dL1 = dL[1];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 1; i < np; i++) {
    for (Size iv = 0; iv < nf; ++iv) {
      const tran         tran_state{K[i - 1][iv], K[i][iv], r[i - 1]};
      diagonal_remainder t_diag_m1;
      T[iv, i] = tran_state(t_diag_m1);
      set_diagonal_remainder(T_diag_m1, iv, i, t_diag_m1);
      diagonal_remainder l_diag_m1;
      L[iv, i] = tran_state.linsrc(l_diag_m1);
      set_diagonal_remainder(L_diag_m1, iv, i, l_diag_m1);
      const propmat kavg = avg(K[i - 1][iv], K[i][iv]);

      for (Size j = 0; j < nq; j++) {
        dT0[iv, i - 1, j] = tran_state.deriv(T[iv, i], -0.5 * r[i - 1] * dK[i - 1][j, iv] - dr0[i - 1, j] * kavg);
        dT1[iv, i, j]     = tran_state.deriv(T[iv, i], -0.5 * r[i - 1] * dK[i][j, iv] - dr1[i - 1, j] * kavg);
        dL0[iv, i - 1, j] = tran_state.linsrc_deriv(-0.5 * r[i - 1] * dK[i - 1][j, iv] - dr0[i - 1, j] * kavg);
        dL1[iv, i, j]     = tran_state.linsrc_deriv(-0.5 * r[i - 1] * dK[i][j, iv] - dr1[i - 1, j] * kavg);
      }
    }
  }
}

void TransmittanceMatrix::linprop(const std::span<const propmat_vector> &K,
                                  const std::span<const propmat_matrix> &dK,
                                  const ConstVectorView                 &r,
                                  const ConstTensor3View                &dr) {
  const Size nf = dT.npages();
  const Size np = dT.nrows();
  const Size nq = dT.ncols();

  auto dT0 = dT[0];
  auto dT1 = dT[1];
  auto dL0 = dL[0];
  auto dL1 = dL[1];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 1; i < np; i++) {
    for (Size iv = 0; iv < nf; ++iv) {
      const propmat     &start = K[i][iv];
      const propmat     &end   = K[i - 1][iv];
      const Numeric      ri    = r[i - 1];
      const tran         tran_state{start, end, ri};
      diagonal_remainder t_diag_m1;
      T[iv, i] = tran_state(t_diag_m1);
      set_diagonal_remainder(T_diag_m1, iv, i, t_diag_m1);
      const scalar_linprop_result state =
          scalar_linprop(ri * start.A(), ri * end.A(), tran_state.a, tran_state.exp_a, tran_state.expm1_a);
      const propmat kavg = avg(start, end);

      diagonal_remainder     l_diag_m1;
      std::optional<muelmat> base_phi;
      std::optional<muelmat> base_q_m1;
      Numeric                base_phi_scalar;
      Numeric                base_phi_prime;
      if (not tran_state.polarized) {
        L[iv, i] = state.value;
        l_diag_m1.fill(state.value_m1);
      } else {
        L[iv, i] = augmented_linprop_value(tran_state,
                                           state,
                                           start,
                                           end,
                                           ri,
                                           &l_diag_m1,
                                           nq == 0 ? nullptr : &base_phi,
                                           nq == 0 ? nullptr : &base_q_m1);
        if (nq != 0) {
          base_phi_scalar = tran_state.a == 0.0 ? 1.0 : tran_state.expm1_a / tran_state.a;
          base_phi_prime  = phi1_scalar_derivative_from_exp(tran_state.a, tran_state.exp_a);
        }
      }
      set_diagonal_remainder(L_diag_m1, iv, i, l_diag_m1);

      for (Size j = 0; j < nq; j++) {
        const Numeric  dr_end   = dr0[i - 1, j];
        const Numeric  dr_start = dr1[i - 1, j];
        const propmat &dk_end   = dK[i - 1][j, iv];
        const propmat &dk_start = dK[i][j, iv];
        dT0[iv, i - 1, j]       = tran_state.deriv(T[iv, i], -0.5 * ri * dk_end - dr_end * kavg);
        dT1[iv, i, j]           = tran_state.deriv(T[iv, i], -0.5 * ri * dk_start - dr_start * kavg);

        if (not tran_state.polarized and not dk_end.is_polarized()) {
          dL0[iv, i - 1, j] = muelmat{scalar_linprop_deriv(state, start.A(), end.A(), dk_end.A(), ri, dr_end, false)};
        } else {
          if (not base_phi.has_value()) {
            base_phi.emplace(tran_state.linsrc());
            base_q_m1.emplace(-(1.0 / 12.0) * to_muelmat(ri * end - ri * start));
            base_phi_scalar = tran_state.a == 0.0 ? 1.0 : tran_state.expm1_a / tran_state.a;
            base_phi_prime  = phi1_scalar_derivative_from_exp(tran_state.a, tran_state.exp_a);
          }
          dL0[iv, i - 1, j] = augmented_linprop_deriv(tran_state,
                                                      state,
                                                      *base_phi,
                                                      *base_q_m1,
                                                      kavg,
                                                      base_phi_scalar,
                                                      base_phi_prime,
                                                      start,
                                                      end,
                                                      dk_end,
                                                      ri,
                                                      dr_end,
                                                      false);
        }

        if (not tran_state.polarized and not dk_start.is_polarized()) {
          dL1[iv, i, j] = muelmat{scalar_linprop_deriv(state, start.A(), end.A(), dk_start.A(), ri, dr_start, true)};
        } else {
          if (not base_phi.has_value()) {
            base_phi.emplace(tran_state.linsrc());
            base_q_m1.emplace(-(1.0 / 12.0) * to_muelmat(ri * end - ri * start));
            base_phi_scalar = tran_state.a == 0.0 ? 1.0 : tran_state.expm1_a / tran_state.a;
            base_phi_prime  = phi1_scalar_derivative_from_exp(tran_state.a, tran_state.exp_a);
          }
          dL1[iv, i, j] = augmented_linprop_deriv(tran_state,
                                                  state,
                                                  *base_phi,
                                                  *base_q_m1,
                                                  kavg,
                                                  base_phi_scalar,
                                                  base_phi_prime,
                                                  start,
                                                  end,
                                                  dk_start,
                                                  ri,
                                                  dr_start,
                                                  true);
        }
      }
    }
  }
}

void TransmittanceMatrix::magop(const std::span<const propmat_vector> &K,
                                const std::span<const propmat_matrix> &dK,
                                const ConstVectorView                 &r,
                                const ConstTensor3View                &dr) {
  const Size nf = dT.npages();
  const Size np = dT.nrows();
  const Size nq = dT.ncols();

  auto dT0 = dT[0];
  auto dT1 = dT[1];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 1; i < np; i++) {
    for (Size iv = 0; iv < nf; ++iv) {
      const propmat     &start = K[i][iv];
      const propmat     &end   = K[i - 1][iv];
      const tran         tr{start, end, r[i - 1], MagnusOperator::magnus};
      diagonal_remainder t_diag_m1;
      T[iv, i] = tr(t_diag_m1);
      set_diagonal_remainder(T_diag_m1, iv, i, t_diag_m1);

      for (Size iq = 0; iq < nq; iq++) {
        dT0[iv, i - 1, iq] = tr.magnus_deriv(T[iv, i], start, end, dK[i - 1][iq, iv], r[i - 1], dr0[i - 1, iq], false);
        dT1[iv, i, iq]     = tr.magnus_deriv(T[iv, i], start, end, dK[i][iq, iv], r[i - 1], dr1[i - 1, iq], true);
      }
    }
  }
}

void TransmittanceMatrix::magop_linsrc(const std::span<const propmat_vector> &K,
                                       const std::span<const propmat_matrix> &dK,
                                       const ConstVectorView                 &r,
                                       const ConstTensor3View                &dr) {
  const Size nf = dT.npages();
  const Size np = dT.nrows();
  const Size nq = dT.ncols();

  auto dT0 = dT[0];
  auto dT1 = dT[1];
  auto dL0 = dL[0];
  auto dL1 = dL[1];
  auto dr0 = dr[0];
  auto dr1 = dr[1];

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 1; i < np; ++i) {
    for (Size iv = 0; iv < nf; ++iv) {
      const propmat     &start = K[i][iv];
      const propmat     &end   = K[i - 1][iv];
      const Numeric      ri    = r[i - 1];
      const tran         tr{start, end, ri, MagnusOperator::magnus};
      diagonal_remainder t_diag_m1;
      T[iv, i] = tr(t_diag_m1);
      set_diagonal_remainder(T_diag_m1, iv, i, t_diag_m1);
      diagonal_remainder     l_diag_m1;
      std::optional<muelmat> base_phi;
      std::optional<muelmat> base_q_m1;
      L[iv, i] = augmented_magnus_source(
          tr, start, end, ri, &l_diag_m1, nq == 0 ? nullptr : &base_phi, nq == 0 ? nullptr : &base_q_m1);
      set_diagonal_remainder(L_diag_m1, iv, i, l_diag_m1);

      for (Size iq = 0; iq < nq; ++iq) {
        dT0[iv, i - 1, iq] = tr.magnus_deriv(T[iv, i], start, end, dK[i - 1][iq, iv], ri, dr0[i - 1, iq], false);
        dT1[iv, i, iq]     = tr.magnus_deriv(T[iv, i], start, end, dK[i][iq, iv], ri, dr1[i - 1, iq], true);
        dL0[iv, i - 1, iq] = augmented_magnus_source_deriv(
            tr, *base_phi, *base_q_m1, start, end, dK[i - 1][iq, iv], ri, dr0[i - 1, iq], false);
        dL1[iv, i, iq] = augmented_magnus_source_deriv(
            tr, *base_phi, *base_q_m1, start, end, dK[i][iq, iv], ri, dr1[i - 1, iq], true);
      }
    }
  }
}

void TransmittanceMatrix::init(const std::span<const propmat_vector> &K,
                               const std::span<const propmat_matrix> &dK,
                               const ConstVectorView                 &r,
                               const ConstTensor3View                &dr,
                               const TransmittanceOption              opt) {
  option = opt;

  ARTS_USER_ERROR_IF(not arr::same_size(K, dK),
                     "K and dK must have the same size: K: {}, dK: {}, r: {}",
                     K.size(),
                     dK.size(),
                     r.size());

  constexpr Size nt = 2;
  const Size     np = K.size();

  if (np == 0) {
    ARTS_USER_ERROR_IF(r.size() != 0 or dr.npages() != nt or dr.nrows() != 0,
                       "Empty K requires empty r and dr shape (2, 0, nq), got r: {:B,}, dr: {:B,}",
                       r.shape(),
                       dr.shape());
    const Size nf = T.nrows();
    const Size nq = dr.ncols();
    T.resize(nf, 0);
    T_diag_m1.resize(nf, 0);
    L.resize(nf, 0);
    L_diag_m1.resize(nf, 0);
    P.resize(nf, 0);
    dT.resize(nt, nf, 0, nq);
    dL.resize(nt, nf, 0, nq);
    return;
  }

  const Size nf                = K.front().size();
  const Size nq                = dK.front().nrows();
  const bool has_linear_source = option == TransmittanceOption::linsrc or option == TransmittanceOption::linprop or
                                 option == TransmittanceOption::magop_linsrc;

  ARTS_USER_ERROR_IF(
      dr.npages() != nt or dr.nrows() != static_cast<Index>(np - 1) or dr.ncols() != static_cast<Index>(nq) or
          r.size() != np - 1,
      "dr and r must have compatible sizes. dr: (nt, np-1, nq)->{:B,}, r: (np-1)->{:B,}, nt: 2, np-1: {}, nq: {}",
      dr.shape(),
      r.shape(),
      np - 1,
      nq);

  ARTS_USER_ERROR_IF(not all_same_shape({nf}, K), "All propmats in K must have same size ({}).", nf);

  ARTS_USER_ERROR_IF(not all_same_shape({nq, nf}, dK), "All propmats in dK must have same shape ({}, {}).", nq, nf);

  switch (option) {
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:
    case TransmittanceOption::magop_linsrc:
      L.resize(nf, np);
      L_diag_m1.resize(nf, np);
      dL.resize(nt, nf, np, nq);
      dL = muelmat::constant(0);
      [[fallthrough]];
    case TransmittanceOption::magop:
    case TransmittanceOption::constant:
      T.resize(nf, np);
      T_diag_m1.resize(nf, np);
      dT.resize(nt, nf, np, nq);
      dT = muelmat::constant(0);
      P.resize(nf, np);
  }

  // Kernels overwrite every layer (i>=1); only the boundary column needs a
  // seed.  Avoid a full extra cache/matrix write over this bandwidth-heavy
  // initialization path.
  for (Size iv = 0; iv < nf; ++iv) {
    T[iv, 0]         = muelmat::id();
    T_diag_m1[iv, 0] = stokvec{};
    if (has_linear_source) {
      L[iv, 0]         = muelmat::id();
      L_diag_m1[iv, 0] = stokvec{};
    }
  }

  switch (option) {
    case TransmittanceOption::constant:     constant(K, dK, r, dr); break;
    case TransmittanceOption::linsrc:       linsrc(K, dK, r, dr); break;
    case TransmittanceOption::linprop:      linprop(K, dK, r, dr); break;
    case TransmittanceOption::magop:        magop(K, dK, r, dr); break;
    case TransmittanceOption::magop_linsrc: magop_linsrc(K, dK, r, dr); break;
  }

  for (Size i = 0; i < nf; i++) {
    P[i, 0] = muelmat::id();
    for (Size j = 1; j < np; j++) { P[i, j] = P[i, j - 1] * T[i, j]; }
  }
}

void TransmittanceMatrix::init(const std::span<const propmat>        &K,
                               const std::span<const propmat_vector> &dK,
                               const ConstVectorView                 &r,
                               const ConstTensor3View                &dr,
                               const TransmittanceOption              opt) {
  option = opt;

  constexpr Size nt = 2;
  constexpr Size nf = 1;
  const Size     np = K.size();

  ARTS_USER_ERROR_IF(not arr::same_size(K, dK),
                     "K and dK must have the same size: K: {}, dK: {}, r: {}",
                     K.size(),
                     dK.size(),
                     r.size());

  if (np == 0) {
    ARTS_USER_ERROR_IF(r.size() != 0 or dr.npages() != nt or dr.nrows() != 0,
                       "Empty K requires empty r and dr shape (2, 0, nq), got r: {:B,}, dr: {:B,}",
                       r.shape(),
                       dr.shape());
    const Size nq = dr.ncols();
    T.resize(nf, 0);
    T_diag_m1.resize(nf, 0);
    L.resize(nf, 0);
    L_diag_m1.resize(nf, 0);
    P.resize(nf, 0);
    dT.resize(nt, nf, 0, nq);
    dL.resize(nt, nf, 0, nq);
    return;
  }

  const Size nq                = dK.front().size();
  const bool has_linear_source = option == TransmittanceOption::linsrc or option == TransmittanceOption::linprop or
                                 option == TransmittanceOption::magop_linsrc;

  ARTS_USER_ERROR_IF(
      dr.npages() != nt or dr.nrows() != static_cast<Index>(np - 1) or dr.ncols() != static_cast<Index>(nq) or
          r.size() != np - 1,
      "dr and r must have compatible sizes. dr: (nt, np-1, nq)->{:B,}, r: (np-1)->{:B,}, nt: {}, np-1: {}, nq: {}",
      dr.shape(),
      r.shape(),
      nt,
      np - 1,
      nq);

  ARTS_USER_ERROR_IF(not all_same_shape({nq}, dK), "All propmats in dK must have same shape ({}).", nq);

  switch (option) {
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:
    case TransmittanceOption::magop_linsrc:
      L.resize(nf, np);
      L_diag_m1.resize(nf, np);
      dL.resize(nt, nf, np, nq);
      dL = muelmat::constant(0);
      [[fallthrough]];
    case TransmittanceOption::magop:
    case TransmittanceOption::constant:
      T.resize(nf, np);
      T_diag_m1.resize(nf, np);
      dT.resize(nt, nf, np, nq);
      dT = muelmat::constant(0);
      P.resize(nf, np);
  }

  T[0, 0]         = muelmat::id();
  T_diag_m1[0, 0] = stokvec{};
  if (has_linear_source) {
    L[0, 0]         = muelmat::id();
    L_diag_m1[0, 0] = stokvec{};
  }

  switch (option) {
    case TransmittanceOption::constant:     constant(K, dK, r, dr); break;
    case TransmittanceOption::linsrc:       linsrc(K, dK, r, dr); break;
    case TransmittanceOption::linprop:      linprop(K, dK, r, dr); break;
    case TransmittanceOption::magop:        magop(K, dK, r, dr); break;
    case TransmittanceOption::magop_linsrc: magop_linsrc(K, dK, r, dr); break;
  }

  P[0, 0] = muelmat::id();
  for (Size j = 1; j < np; j++) { P[0, j] = P[0, j - 1] * T[0, j]; }
}

void TransmittanceMatrix::check(Size np, Size nq, Size nf, const std::string_view caller) const {
  switch (option) {
    case TransmittanceOption::magop:
    case TransmittanceOption::constant:
      ARTS_USER_ERROR_IF(
          not same_shape({nf, np}, T, P) or (T_diag_m1.size() != 0 and not same_shape({nf, np}, T_diag_m1)),
          R"(Mismatched shapes in Transmittance in {6}:

T: {0:B,}
T_diag_m1: {1:B,}
P: {2:B,} - optional on nq > 0, nq = {5}

Expected shapes: ({3}, {4})
)",
          T.shape(),
          T_diag_m1.shape(),
          P.shape(),
          nf,
          np,
          nq,
          caller);

      ARTS_USER_ERROR_IF(not same_shape({2, nf, np, nq}, dT),
                         R"(Mismatched shapes in Transmittance in {5}:

dT: {0:B,}

Expected shape: (2, {2}, {3}, {4}
)",
                         dT.shape(),
                         dL.shape(),
                         nf,
                         np,
                         nq,
                         caller)
      break;
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:
    case TransmittanceOption::magop_linsrc:
      ARTS_USER_ERROR_IF(not same_shape({nf, np}, T, L, P) or
                             (T_diag_m1.size() != 0 and not same_shape({nf, np}, T_diag_m1)) or
                             (L_diag_m1.size() != 0 and not same_shape({nf, np}, L_diag_m1)),
                         R"(Mismatched shapes in Transmittance in {8}:

T: {0:B,}
T_diag_m1: {1:B,}
L: {2:B,}
L_diag_m1: {3:B,}
P: {4:B,} - optional on nq > 0, nq = {7}

Expected shapes: ({5}, {6})
)",
                         T.shape(),
                         T_diag_m1.shape(),
                         L.shape(),
                         L_diag_m1.shape(),
                         P.shape(),
                         nf,
                         np,
                         nq,
                         caller);

      ARTS_USER_ERROR_IF(not same_shape({2, nf, np, nq}, dT, dL),
                         R"(Mismatched shapes in Transmittance in {5}:

dT: {0:B,}
dL: {1:B,}

Expected shape: (2, {2}, {3}, {4}
)",
                         dT.shape(),
                         dL.shape(),
                         nf,
                         np,
                         nq,
                         caller)
      break;
  }
}

[[nodiscard]] std::array<Size, 3> TransmittanceMatrix::shape() const noexcept {
  return {static_cast<Size>(dT.npages()), static_cast<Size>(dT.nrows()), static_cast<Size>(dT.ncols())};
}
}  // namespace rtepack
