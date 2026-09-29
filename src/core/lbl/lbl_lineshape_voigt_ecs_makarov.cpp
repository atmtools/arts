#include "lbl_lineshape_voigt_ecs_makarov.h"

#include <arts_conversions.h>
#include <wigner_functions.h>

#include <cmath>

#include "lbl_lineshape_voigt_ecs.h"

namespace lbl::voigt::ecs::makarov {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

namespace {
void validate_rotational_pair(const Rational Ju, const Rational Jl, const Rational N) {
  ARTS_USER_ERROR_IF(
      N.denom != 1 or N <= 0 or iseven(N), "Makarov O2-66 ECS requires a positive odd integer N, got {}", N);
  ARTS_USER_ERROR_IF(Ju.denom != 1 or Jl.denom != 1 or Ju < 0 or Jl < 0 or abs(Ju - N) > 1 or abs(Jl - N) > 1,
                     "Makarov O2-66 ECS requires integer J >= 0 with |J-N| <= 1; got Ju={}, Jl={}, N={}",
                     Ju,
                     Jl,
                     N);
  ARTS_USER_ERROR_IF(abs(Ju - Jl) > 1 or (Ju == 0 and Jl == 0),
                     "Makarov O2-66 ECS requires a rank-1 transition with |Ju-Jl| <= 1 and excludes J=0 to J=0; "
                     "got Ju={}, Jl={}",
                     Ju,
                     Jl);
}

void validate_band_id(const QuantumIdentifier& bnd_qid) {
  ARTS_USER_ERROR_IF(bnd_qid.isot != "O2-66"_isot, "Makarov ECS currently supports only O2-66, got {}", bnd_qid.isot);
  const auto& S = bnd_qid.state.at(QuantumNumberType::S);
  ARTS_USER_ERROR_IF(S.upper != 1 or S.lower != 1,
                     "Makarov O2-66 ECS requires electron spin S=1 in both states; got upper={}, lower={}",
                     S.upper,
                     S.lower);
}

void validate_rotational_line(const rotational_line& ln) {
  ARTS_USER_ERROR_IF(ln.Nu != ln.Nl,
                     "Makarov O2-66 ECS implements microwave transitions with unchanged N; got upper={}, lower={}",
                     ln.Nu,
                     ln.Nl);
  validate_rotational_pair(ln.Ju, ln.Jl, ln.Nu);
}

Numeric wig3(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER3(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(
    const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& e, const Rational& f) {
  return WIGNER6(a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}
// Fill all off-diagonal pairs and manage Wigner scratch for the complete loop.
void coupling_kernel(MatrixView                       W,
                     std::span<const rotational_line> lines,
                     const Rational                   Si,
                     const Rational                   Sf,
                     const basis_data&                basis,
                     const Vector&                    e0,
                     const Numeric                    T,
                     const int                        maxL,
                     Tensor3View                      dW,
                     ConstVectorView                  dT,
                     ConstMatrixView                  dQ,
                     ConstMatrixView                  dOmega) {
  using Conversion::kelvin2joule;
  const auto& Q       = basis.Q;
  const auto& Om      = basis.Omega;
  const Index nq      = dW.npages();
  const auto  tangent = [](ConstMatrixView values, Index q, Index i) -> Numeric {
    return values.empty() ? 0.0 : values[q, i];
  };
  Vector     dsum(nq);
  const auto bk = [](const Rational& r) -> Numeric { return sqrtr(2 * r + 1); };
  const Size n  = lines.size();

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    const auto& [Ji, Jf, Ni, Nf] = lines[i];

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;
      const auto& [Ji_p, Jf_p, Ni_p, Nf_p] = lines[j];

      if (Jf_p > Jf) continue;

      // Tran etal 2006 symbol with modifications:
      //    1) [Ji] * [Ji_p] instead of [Ji_p] ^ 2 in partial accordance with Makarov etal 2013
      Numeric sum       = 0;
      dsum              = 0;
      const Numeric scl = (iseven(Ji_p + Ji + 1) ? 1 : -1) * bk(Ni) * bk(Nf) * bk(Nf_p) * bk(Ni_p) * bk(Jf) * bk(Jf_p) *
                          bk(Ji) * bk(Ji_p);
      const auto [L0, L1] =
          wigner_limits(wigner3j_limits<3>(Ni_p, Ni), {Rational(2), Rational{std::numeric_limits<Index>::max()}});
      for (Rational L = L0; L <= L1; L += 2) {
        const Numeric a  = wig3(Ni_p, Ni, L, Rational{0}, Rational{0}, Rational{0});
        const Numeric b  = wig3(Nf_p, Nf, L, Rational{0}, Rational{0}, Rational{0});
        const Numeric c  = wig6(L, Ji, Ji_p, Si, Ni_p, Ni);
        const Numeric d  = wig6(L, Jf, Jf_p, Sf, Nf_p, Nf);
        const Numeric e  = wig6(L, Ji, Ji_p, Rational{1}, Jf_p, Jf);
        sum             += a * b * c * d * e * Numeric(2 * L + 1) * Q[L.toIndex()] / Om[L.toIndex()];
        for (Index q = 0; q < nq; ++q) {
          const Index idx  = L.toIndex();
          dsum[q]         += a * b * c * d * e * Numeric(2 * L + 1) *
                             (tangent(dQ, q, idx) / Om[idx] - Q[idx] / Om[idx] * (tangent(dOmega, q, idx) / Om[idx]));
        }
      }
      for (Index q = 0; q < nq; ++q) {
        dsum[q] = dsum[q] * (scl * Om[Ni.toIndex()]) + sum * (scl * tangent(dOmega, q, Ni.toIndex()));
      }
      sum *= scl * Om[Ni.toIndex()];

      // Add to W and rescale to upwards element by the populations.
      W[i, j] = sum;
      W[j, i] = sum * std::exp((e0[j] - e0[i]) / kelvin2joule(T));
      if (nq != 0) {
        const Numeric exponent = (e0[j] - e0[i]) / kelvin2joule(T);
        const Numeric balance  = std::exp(exponent);
        for (Index q = 0; q < nq; ++q) {
          const Numeric dexponent = dT.empty() ? 0.0 : -exponent * dT[q] / T;
          dW[q, i, j]             = dsum[q];
          dW[q, j, i]             = (dsum[q] + sum * dexponent) * balance;
        }
      }
    }
  }
  arts_wigner_thread_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")
}
}  // namespace

Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N) {
  validate_rotational_pair(Ju, Jl, N);
  return (iseven(Jl + N) ? 1 : -1) * sqrtr(6 * (2 * Jl + 1) * (2 * Ju + 1)) *
         wigner6j(Rational{1}, Rational{1}, Rational{1}, Jl, Ju, N);
};

namespace {
// O2-66 constants in MHz, from Tretyakov et al., JMS 231 (2005), Table 3.
constexpr Numeric B0 = 43100.4425;

constexpr Numeric rotational_energy_mhz(const Rational N) {
  const Numeric X = Numeric(N * (N + 1));
  return B0 * X - 0.145123 * Math::pow2(X) + 3.8e-08 * Math::pow3(X);
}

// Approximate spin-triplet expressions. These retain the existing treatment
// of centrifugal distortion and spin rotation; they are not a replacement
// for catalogue transition frequencies (low-N residuals reach about 24 MHz).
constexpr Numeric level_energy_mhz(const Rational N, const Rational J) {
  const Numeric XN         = Numeric(N);
  const Numeric X          = XN * (XN + 1);
  const Numeric lambda     = 59501.3435 + 0.058369 * X + 2.899e-07 * Math::pow2(X);
  const Numeric gamma      = -252.58633 - 2.4344e-04 * X - 1.45e-09 * Math::pow2(X);
  const Numeric rotational = rotational_energy_mhz(N);

  if (J < N) {
    // J=0 has no N=-1 mixing partner. Its energy is rotational - 2*lambda - gamma.
    // Van Vleck (1947), p. 414, footnote 3: doi:10.1103/PhysRev.71.413.
    // The old special case dropped the lambda-B0 contribution, putting the
    // N=1, J=1 <- J=0 splitting at 102.349 GHz instead of 118.750 GHz.
    if (N == 1) return rotational - 2 * lambda - gamma;
    return rotational - (lambda + B0 * (2 * XN - 1) + gamma * XN) +
           std::sqrt(Math::pow2(B0 * (2 * XN - 1)) + Math::pow2(lambda) - 2 * B0 * lambda);
  }
  if (J > N)
    return rotational - (lambda - B0 * (2 * XN + 3) - gamma * (XN + 1)) -
           std::sqrt(Math::pow2(B0 * (2 * XN + 3)) + Math::pow2(lambda) - 2 * B0 * lambda);
  return rotational;
}

constexpr Numeric ground_energy_mhz = level_energy_mhz(Rational{1}, Rational{0});
}  // namespace

Numeric rotational_energy(const Rational N) {
  // The ECS reference rotor has no resolved spin splitting. Preserve a common
  // ground-state zero with level_energy, including in Q's absolute energy.
  return Conversion::mhz2joule(rotational_energy_mhz(N) - ground_energy_mhz);
}

Numeric level_energy(const Rational N, const Rational J) {
  ARTS_USER_ERROR_IF(N.denom != 1 or N <= 0 or iseven(N) or J.denom != 1 or J < 0 or abs(J - N) > 1,
                     "O2-66 level energies require positive odd integer N and integer J >= 0 with |J-N| <= 1; "
                     "got N={}, J={}",
                     N,
                     J)
  return Conversion::mhz2joule(level_energy_mhz(N, J) - ground_energy_mhz);
}

void prepare_energies(energy_data& energies, const QuantumIdentifier& qid, std::span<const rotational_line> lines) {
  validate_band_id(qid);
  const auto& S = qid.state.at(QuantumNumberType::S);
  Rational    maxJ{0}, maxN{0};
  energies.e0.resize(lines.size());
  for (Size i = 0; i < lines.size(); ++i) {
    const auto& ln = lines[i];
    maxJ           = std::max({maxJ, ln.Ju, ln.Jl});
    maxN           = std::max({maxN, ln.Nu, ln.Nl});
    energies.e0[i] = level_energy(ln.Nl, ln.Jl);
  }
  const std::array rats{maxJ, maxN, Rational{S.upper}, Rational{S.lower}};
  const int        maxL = wigner_init_size(rats);
  prepare_rotational_ladder(energies, maxL, rotational_energy);
}

void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd) {
  validate_band_id(bnd_qid);
  for (const auto& ln : bnd) {
    const auto& J = ln.qn.at(QuantumNumberType::J);
    const auto& N = ln.qn.at(QuantumNumberType::N);
    validate_rotational_line({J.upper, J.lower, N.upper, N.lower});
  }
}

void relaxation_matrix_offdiagonal(MatrixView&                      W,
                                   const QuantumIdentifier&         bnd_qid,
                                   std::span<const rotational_line> lines,
                                   Numeric                          T0,
                                   const SpeciesEnum                broadening_species,
                                   const linemixing::species_data&  rovib_data,
                                   const Vector&                    dipr,
                                   const energy_data&               energies,
                                   const AtmPoint&                  atm,
                                   Tensor3View                      dW,
                                   ConstVectorView                  dT,
                                   ConstMatrixView                  dQ,
                                   ConstMatrixView                  dOmega) try {
  ARTS_USER_ERROR_IF((dW.npages() != 0 and (dW.nrows() != W.nrows() or dW.ncols() != W.ncols())) or
                         (not dT.empty() and dT.size() != static_cast<Size>(dW.npages())),
                     "Inconsistent Makarov ECS derivative dimensions")
  if (lines.empty()) return;
  validate_band_id(bnd_qid);
  const auto& e0 = energies.e0;

  const auto n = lines.size();
  ARTS_USER_ERROR_IF(
      e0.size() != n or dipr.size() != n or W.nrows() != static_cast<Index>(n) or W.ncols() != static_cast<Index>(n),
      "Inconsistent Makarov ECS kernel dimensions")

  auto&          S  = bnd_qid.state.at(QuantumNumberType::S);
  const Rational Si = S.upper;
  const Rational Sf = S.lower;

  Rational maxJ{0}, maxN{0};
  for (const auto& ln : lines) {
    maxJ = std::max({maxJ, ln.Ju, ln.Jl});
    maxN = std::max({maxN, ln.Nu, ln.Nl});
  }

  const std::array rats{maxJ, maxN, Si, Sf};
  const int        maxL  = wigner_init_size(rats);
  const auto       basis = prepare_basis(maxL, energies, rovib_data, T0, bnd_qid.isot, broadening_species, atm);
  ARTS_USER_ERROR_IF(
      (not dQ.empty() and (dQ.nrows() != dW.npages() or dQ.ncols() != static_cast<Index>(basis.Q.size()))) or
          (not dOmega.empty() and
           (dOmega.nrows() != dW.npages() or dOmega.ncols() != static_cast<Index>(basis.Omega.size()))),
      "Inconsistent Makarov ECS basis derivative dimensions")

  coupling_kernel(W, lines, Si, Sf, basis, e0, atm.temperature, maxL, dW, dT, dQ, dOmega);

  apply_sum_rule(W, dipr, e0, atm.temperature, dW, dT);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs::makarov
