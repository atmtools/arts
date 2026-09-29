#include "lbl_lineshape_voigt_ecs_hadded.h"

#include <arts_constants.h>
#include <wigner_functions.h>

#include <array>
#include <cerrno>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>

namespace lbl::voigt::ecs::hadded {
namespace {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

Numeric parity(Index n) { return n % 2 == 0 ? 1.0 : -1.0; }

void validate_state(Index J, Index K) {
  ARTS_USER_ERROR_IF(J < 0 or K < 0 or K > J or J > (std::numeric_limits<int>::max() - 1) / 24,
                     "NH3 ECS requires integer 0 <= K <= J within the Wigner range; got J={}, K={}",
                     J,
                     K)
}

void validate_line(const rotational_line& line) {
  const auto& u = line.upper;
  const auto& l = line.lower;
  validate_state(u.J, u.K);
  validate_state(l.J, l.K);
  const auto valid_inversion = [](inversion x) { return x == inversion::symmetric or x == inversion::antisymmetric; };
  ARTS_USER_ERROR_IF(not valid_inversion(u.symmetry) or not valid_inversion(l.symmetry) or u.symmetry == l.symmetry or
                         u.K != l.K or std::abs(u.J - l.J) > 1 or (u.J == l.J and u.K == 0),
                     "NH3 ECS requires an allowed rank-1 parallel transition with equal K, "
                     "opposite inversion symmetry, and no K=0 Q line")
}

// A stable tie-break for equal lower-state energies, independent of matrix order.
auto state_key(const rotational_line& line) {
  return std::tuple{line.lower.J, line.lower.K, line.lower.symmetry, line.upper.J, line.upper.K, line.upper.symmetry};
}

Numeric epsilon(const rotational_state& state) {
  if (state.K == 0) return 0;
  return parity(state.J) * (state.symmetry == inversion::antisymmetric ? 1 : -1);
}

Numeric tangent(ConstMatrixView values, Index q, Index i) { return values.empty() ? 0.0 : values[q, i]; }

// Only the shapes are checked here.  A nonfinite tangent propagates into dW
// and is reported once for the result rather than scanned for per element.
void validate_derivatives(ConstMatrixView values, Index targets, Index inputs) {
  ARTS_USER_ERROR_IF(not values.empty() and (values.nrows() != targets or values.ncols() != inputs),
                     "Inconsistent NH3 ECS derivative dimensions")
}

void validate_derivatives(ConstVectorView values, Index targets) {
  ARTS_USER_ERROR_IF(not values.empty() and values.size() != static_cast<Size>(targets),
                     "Inconsistent NH3 ECS derivative dimensions")
}

void validate_output(MatrixView W, Tensor3View dW, Index n){
    ARTS_USER_ERROR_IF(W.nrows() != n or W.ncols() != n or (dW.npages() != 0 and (dW.nrows() != n or dW.ncols() != n)),
                       "Inconsistent NH3 ECS matrix dimensions")}

Numeric wig3(Index a, Index b, Index c, Index d, Index e, Index f) {
  return WIGNER3(static_cast<int>(2 * a),
                 static_cast<int>(2 * b),
                 static_cast<int>(2 * c),
                 static_cast<int>(2 * d),
                 static_cast<int>(2 * e),
                 static_cast<int>(2 * f));
}
Numeric wig6(Index a, Index b, Index c, Index d, Index e, Index f) {
  return WIGNER6(static_cast<int>(2 * a),
                 static_cast<int>(2 * b),
                 static_cast<int>(2 * c),
                 static_cast<int>(2 * d),
                 static_cast<int>(2 * e),
                 static_cast<int>(2 * f));
}

using channel_index = std::map<std::array<Index, 3>, Index>;

// Hadded et al. (2002): the IOS angular kernel (Eq. 10), with detailed balance
// (Eq. 16) and ECS energy corrections (Eqs. 17-18). This call owns Wigner scratch for all pairs
// and all Jacobian targets; no public Wigner wrappers are called in the loop.
void coupling_kernel(MatrixView                       W,
                     std::span<const rotational_line> lines,
                     const channel_index&             channels,
                     ConstVectorView                  corrected_Q,
                     ConstMatrixView                  dcorrected_Q,
                     ConstVectorView                  e0,
                     ConstVectorView                  Omega_line,
                     Numeric                          T,
                     Index                            max_rank,
                     int                              scratch_size,
                     Tensor3View                      dW,
                     const derivative_data&           derivatives) {
  const Index n  = static_cast<Index>(lines.size());
  const Index nq = dW.npages();
  Vector      dsum(nq);
  errno = 0;
  {
    arts_wigner_thread_init(scratch_size);
    struct scratch_guard {
      ~scratch_guard() { arts_wigner_thread_free(); }
    } guard;

    for (Index a = 0; a < n; ++a) {
      for (Index b = a + 1; b < n; ++b) {
        // NH3 nuclear-spin species cannot be connected by these collisions.
        if ((lines[a].lower.K % 3 == 0) != (lines[b].lower.K % 3 == 0)) continue;
        Index from = a, to = b;
        if (e0[from] < e0[to] or (e0[from] == e0[to] and state_key(lines[from]) < state_key(lines[to])))
          std::swap(from, to);
        const auto&   i      = lines[from].lower;
        const auto&   f      = lines[from].upper;
        const auto&   ip     = lines[to].lower;
        const auto&   fp     = lines[to].upper;
        const Numeric NiNf   = i.K == 0 ? 1.0 : 0.5;
        const Numeric NipNfp = ip.K == 0 ? 1.0 : 0.5;
        const Numeric scale =
            -NiNf * NipNfp * Numeric(2 * ip.J + 1) * std::sqrt(Numeric(2 * fp.J + 1) * Numeric(2 * f.J + 1));
        const std::array Mi{i.K - ip.K, i.K + ip.K};
        const std::array Mf{f.K - fp.K, f.K + fp.K};
        const Index      L0  = std::max(std::abs(i.J - ip.J), std::abs(f.J - fp.J));
        const Index      L1  = std::min({i.J + ip.J, f.J + fp.J, max_rank});
        Numeric          sum = 0;
        dsum                 = 0;
        for (Index L = L0; L <= L1; ++L) {
          const Numeric    Pi = parity(ip.J + ip.K + i.J + i.K + L);
          const Numeric    Pf = parity(fp.J + fp.K + f.J + f.K + L);
          const std::array Ci{1 + epsilon(i) * epsilon(ip) * Pi, epsilon(ip) + epsilon(i) * Pi};
          const std::array Cf{1 + epsilon(f) * epsilon(fp) * Pf, epsilon(fp) + epsilon(f) * Pf};
          const Numeric    sixj = wig6(i.J, f.J, 1, fp.J, ip.J, L);
          if (sixj == 0) continue;
          const Numeric angular_scale = scale * parity(fp.J + f.J + fp.K + ip.K + 1 + L) * sixj;
          for (Index si = 0; si < 2; ++si) {
            if (Ci[si] == 0 or std::abs(Mi[si]) > L) continue;
            const Numeric initial = Ci[si] * wig3(ip.J, L, i.J, si == 0 ? ip.K : -ip.K, Mi[si], -i.K);
            for (Index sf = 0; sf < 2; ++sf) {
              if (Cf[sf] == 0 or std::abs(Mf[sf]) > L) continue;
              const auto entry = channels.find({L, Mi[si], Mf[sf]});
              if (entry == channels.end()) continue;
              const Numeric final        = Cf[sf] * wig3(fp.J, L, f.J, sf == 0 ? fp.K : -fp.K, Mf[sf], -f.K);
              const Numeric coefficient  = angular_scale * initial * final;
              const Index   c            = entry->second;
              sum                       += coefficient * corrected_Q[c];
              for (Index q = 0; q < nq; ++q) dsum[q] += coefficient * dcorrected_Q[q, c];
            }
          }
        }
        // Eqs. 17-18 use Omega >= 1, unlike the reciprocal CO2/O2 factor.
        const Numeric down     = sum / Omega_line[from];
        const Numeric exponent = (e0[to] - e0[from]) / (Constant::k * T);
        const Numeric balance  = Numeric(2 * i.J + 1) / Numeric(2 * ip.J + 1) * std::exp(exponent);
        W[to, from]            = down;
        W[from, to]            = down * balance;
        for (Index q = 0; q < nq; ++q) {
          const Numeric ddown = (dsum[q] - down * tangent(derivatives.dOmega_line, q, from)) / Omega_line[from];
          const Numeric dT    = derivatives.dT.empty() ? 0.0 : derivatives.dT[q];
          const Numeric dexponent =
              (tangent(derivatives.de0, q, to) - tangent(derivatives.de0, q, from)) / (Constant::k * T) -
              exponent * dT / T;
          dW[q, to, from] = ddown;
          dW[q, from, to] = (ddown + down * dexponent) * balance;
        }
      }
    }
  }
  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute NH3 ECS Wigner symbols")
}
}  // namespace

Numeric rotational_energy(Index J, Index K, Numeric B, Numeric C) {
  validate_state(J, K);
  ARTS_USER_ERROR_IF(not std::isfinite(B) or B <= 0 or not std::isfinite(C) or C <= 0,
                     "NH3 rotational constants must be positive and finite")
  const Numeric j = Numeric(J), k = Numeric(K);
  const Numeric energy = B * (j * (j + 1) - k * k) + C * k * k;
  ARTS_USER_ERROR_IF(not std::isfinite(energy), "NH3 rotational energy overflowed")
  return energy;
}

Numeric reduced_dipole(const rotational_line& line) {
  validate_line(line);
  const auto& u = line.upper;
  const auto& l = line.lower;
  return std::sqrt(Numeric(2 * u.J + 1)) * parity(u.J + u.K) *
         wigner3j(Rational{u.J}, Rational{1}, Rational{l.J}, Rational{u.K}, Rational{l.K - u.K}, Rational{-l.K});
}

void adiabatic_factors(VectorView      Omega,
                       ConstVectorView gap,
                       Numeric         duration,
                       MatrixView      dOmega,
                       ConstMatrixView dgap,
                       ConstVectorView dduration) {
  const Index n  = static_cast<Index>(gap.size());
  const Index nq = dOmega.nrows();
  ARTS_USER_ERROR_IF(Omega.size() != gap.size() or (nq != 0 and dOmega.ncols() != n),
                     "Inconsistent NH3 adiabatic-factor dimensions")
  ARTS_USER_ERROR_IF(not std::isfinite(duration) or duration < 0,
                     "NH3 collision duration must be finite and nonnegative")
  validate_derivatives(dgap, nq, n);
  validate_derivatives(dduration, nq);
  for (Index i = 0; i < n; ++i) {
    // A negative gap silently reverses the correction, so it is checked; an
    // overflow to a nonfinite factor shows up in the result instead.
    ARTS_USER_ERROR_IF(gap[i] < 0, "NH3 adiabatic energy gaps must be nonnegative")
    const Numeric x = duration * (gap[i] / Constant::h_bar);
    const Numeric a = 1 + x * x / 24;
    Omega[i]        = a * a;
    for (Index q = 0; q < nq; ++q) {
      const Numeric dt = dduration.empty() ? 0.0 : dduration[q];
      const Numeric dx = (duration * tangent(dgap, q, i) + gap[i] * dt) / Constant::h_bar;
      dOmega[q, i]     = a * x * dx / 6;
    }
  }
}

void relaxation_matrix_offdiagonal(MatrixView                       W,
                                   std::span<const rotational_line> lines,
                                   const basis_data&                basis,
                                   ConstVectorView                  e0,
                                   ConstVectorView                  Omega_line,
                                   Numeric                          T,
                                   Tensor3View                      dW,
                                   const derivative_data&           derivatives) {
  const Index n  = static_cast<Index>(lines.size());
  const Index nc = static_cast<Index>(basis.channels.size());
  const Index nq = dW.npages();
  validate_output(W, dW, n);
  ARTS_USER_ERROR_IF(e0.size() != lines.size() or Omega_line.size() != lines.size() or
                         basis.Q.size() != basis.channels.size() or basis.Omega.size() != basis.channels.size(),
                     "Inconsistent NH3 ECS input dimensions")
  ARTS_USER_ERROR_IF(not std::isfinite(T) or T <= 0, "NH3 ECS requires positive finite temperature")
  validate_derivatives(derivatives.dT, nq);
  validate_derivatives(derivatives.de0, nq, n);
  validate_derivatives(derivatives.dQ, nq, nc);
  validate_derivatives(derivatives.dOmega_basis, nq, nc);
  validate_derivatives(derivatives.dOmega_line, nq, n);
  Index maxJ = 1;
  for (Index i = 0; i < n; ++i) {
    validate_line(lines[i]);
    maxJ = std::max({maxJ, lines[i].upper.J, lines[i].lower.J});
    ARTS_USER_ERROR_IF(Omega_line[i] < 1, "NH3 ECS requires Omega >= 1")
    for (Index j = 0; j < i; ++j)
      ARTS_USER_ERROR_IF(state_key(lines[i]) == state_key(lines[j]),
                         "Duplicate rotational/inversion transition in NH3 parallel band")
  }
  channel_index channels;
  Vector        corrected_Q(nc);
  Matrix        dcorrected_Q(nq, nc);
  Index         max_rank = 0;
  for (Index c = 0; c < nc; ++c) {
    const auto& [L, Mi, Mf] = basis.channels[c];
    ARTS_USER_ERROR_IF(L < 0 or L > (std::numeric_limits<int>::max() - 1) / 12 or Mi < -L or Mi > L or Mf < -L or
                           Mf > L or Mi % 3 != 0 or Mf % 3 != 0,
                       "Invalid NH3 ECS channel (L,Mi,Mf)=({},{},{})",
                       L,
                       Mi,
                       Mf)
    ARTS_USER_ERROR_IF(not channels.emplace(std::array{L, Mi, Mf}, c).second, "Duplicate NH3 ECS channel")
    ARTS_USER_ERROR_IF(basis.Omega[c] < 1, "NH3 ECS requires basis Omega >= 1")
    max_rank       = std::max(max_rank, L);
    corrected_Q[c] = basis.Q[c] * basis.Omega[c];
    for (Index q = 0; q < nq; ++q) {
      dcorrected_Q[q, c] =
          tangent(derivatives.dQ, q, c) * basis.Omega[c] + basis.Q[c] * tangent(derivatives.dOmega_basis, q, c);
    }
  }
  for (Index i = 0; i < n; ++i) {
    for (Index j = 0; j < n; ++j) {
      if (i == j) continue;
      W[i, j] = 0;
      for (Index q = 0; q < nq; ++q) dW[q, i, j] = 0;
    }
  }
  if (n < 2 or nc == 0) return;
  const Index max_arg = std::max(maxJ, std::min(2 * maxJ, max_rank));
  ARTS_USER_ERROR_IF(not is_wigner3_ready(Rational{max_arg}) or not is_wigner6_ready(Rational{max_arg}),
                     "Initialize Wigner tables before calculating NH3 ECS")
  coupling_kernel(W,
                  lines,
                  channels,
                  corrected_Q,
                  dcorrected_Q,
                  e0,
                  Omega_line,
                  T,
                  max_rank,
                  static_cast<int>(2 * max_arg + 1),
                  dW,
                  derivatives);
}

void sum_rule_diagonal(MatrixView W, ConstVectorView dipr, Tensor3View dW) {
  const Index n = static_cast<Index>(dipr.size());
  validate_output(W, dW, n);
  for (Numeric dipole : dipr) ARTS_USER_ERROR_IF(dipole == 0, "NH3 sum rule requires nonzero dipoles")
  for (Index i = 0; i < n; ++i) {
    Numeric width = 0;
    for (Index j = 0; j < n; ++j)
      if (j != i) width -= dipr[j] * W[i, j] / dipr[i];
    W[i, i] = width;
    for (Index q = 0; q < dW.npages(); ++q) {
      Numeric derivative = 0;
      for (Index j = 0; j < n; ++j)
        if (j != i) derivative -= dipr[j] * dW[q, i, j] / dipr[i];
      dW[q, i, i] = derivative;
    }
  }
}
}  // namespace lbl::voigt::ecs::hadded
