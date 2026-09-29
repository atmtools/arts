#include "lbl_lineshape_voigt_ecs_hartmann.h"

#include <arts_conversions.h>
#include <atm.h>
#include <wigner_functions.h>

#include <cmath>

#include "lbl_lineshape_voigt_ecs.h"

namespace lbl::voigt::ecs::hartmann {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

namespace {
void validate_rotational_state(const Rational J, const Rational l) {
  ARTS_USER_ERROR_IF(J.denom != 1 or l.denom != 1 or J < 0 or abs(l) > J,
                     "Hartmann ECS requires integer J >= 0 and integer |l| <= J; got J={}, l={}",
                     J,
                     l);
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
// Band labels are already swapped; only the line angular labels are swapped here.
void coupling_kernel(MatrixView                       W,
                     std::span<const rotational_line> lines,
                     const Rational                   li,
                     const Rational                   lf,
                     const bool                       swap_order,
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
  const Size n = lines.size();

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    Rational Ji = lines[i].Ju;
    Rational Jf = lines[i].Jl;
    if (swap_order) std::swap(Ji, Jf);

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;
      Rational Ji_p = lines[j].Ju;
      Rational Jf_p = lines[j].Jl;
      if (swap_order) std::swap(Ji_p, Jf_p);

      // Select the direction in the kernel's angular convention.
      if (Jf_p > Jf) continue;

      Index L         = std::max(std::abs((Ji - Ji_p).toIndex()), std::abs((Jf - Jf_p).toIndex()));
      L              += L % 2;
      const Index Lf  = std::min((Ji + Ji_p).toIndex(), (Jf + Jf_p).toIndex());

      Numeric sum = 0;
      dsum        = 0;
      for (; L <= Lf; L += 2) {
        const Numeric a  = wig3(Ji, Ji_p, Rational{L}, li, -li, Rational{0});
        const Numeric b  = wig3(Jf, Jf_p, Rational{L}, lf, -lf, Rational{0});
        const Numeric c  = wig6(Ji, Jf, Rational{1}, Jf_p, Ji_p, Rational{L});
        sum             += a * b * c * Numeric(2 * L + 1) * Q[L] / Om[L];
        for (Index q = 0; q < nq; ++q) {
          dsum[q] += a * b * c * Numeric(2 * L + 1) *
                     (tangent(dQ, q, L) / Om[L] - Q[L] / Om[L] * (tangent(dOmega, q, L) / Om[L]));
        }
      }
      const Numeric ECS = Om[Ji.toIndex()];
      const Numeric scl = ECS * Numeric(2 * Ji_p + 1) * sqrtr((2 * Jf + 1) * (2 * Jf_p + 1));
      for (Index q = 0; q < nq; ++q) {
        const Numeric dscl =
            tangent(dOmega, q, Ji.toIndex()) * Numeric(2 * Ji_p + 1) * sqrtr((2 * Jf + 1) * (2 * Jf_p + 1));
        dsum[q] = dsum[q] * scl + sum * dscl;
      }
      sum *= scl;

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

void validate_isotopologue(const SpeciesIsotope& isot) {
  ARTS_USER_ERROR_IF(
      isot != "CO2-626"_isot, "Hartmann ECS currently supports rotational energies only for CO2-626, got {}", isot)
}
}  // namespace

Numeric reduced_dipole(const Rational Jf, const Rational Ji, const Rational lf, const Rational li, const Rational k) {
  validate_rotational_state(Jf, lf);
  validate_rotational_state(Ji, li);
  if (not iseven(Jf + lf + 1)) return -sqrtr(2 * Jf + 1) * wigner3j(Jf, k, Ji, lf, li - lf, -li);
  return +sqrtr(2 * Jf + 1) * wigner3j(Jf, k, Ji, lf, li - lf, -li);
}

Numeric rotational_energy(const Rational J) { return Conversion::kaycm2joule(0.39021) * Numeric(J * (J + 1)); }

Numeric level_energy(const Rational J) {
  validate_rotational_state(J, Rational{0});
  return rotational_energy(J);
}

void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd) {
  validate_isotopologue(bnd_qid.isot);
  const auto& l2 = bnd_qid.state.at(QuantumNumberType::l2);

  for (const auto& ln : bnd) {
    const auto& J = ln.qn.at(QuantumNumberType::J);
    validate_rotational_state(J.upper, l2.upper);
    validate_rotational_state(J.lower, l2.lower);
  }

  // This linear-rotor kernel has no elastic L=0 basis rate.  Distinct
  // transitions with the same rotational pair require additional state labels
  // and collision dynamics; evaluating the power law at L=0 is singular.
  for (Size i = 0; i < bnd.size(); ++i) {
    const auto& Ji = bnd.lines[i].qn.at(QuantumNumberType::J);
    for (Size j = i + 1; j < bnd.size(); ++j) {
      const auto& Jj = bnd.lines[j].qn.at(QuantumNumberType::J);
      ARTS_USER_ERROR_IF(Ji.upper == Jj.upper and Ji.lower == Jj.lower,
                         "Hartmann ECS does not support distinct lines with the same rotational pair "
                         "(J upper={}, J lower={}); an elastic L=0 collision model is required.",
                         Ji.upper,
                         Ji.lower)
    }
  }
}

void prepare_energies(energy_data& energies, const QuantumIdentifier& qid, std::span<const rotational_line> lines) {
  validate_isotopologue(qid.isot);
  Rational maxJ{0};
  energies.e0.resize(lines.size());
  for (Size i = 0; i < lines.size(); ++i) {
    const auto& ln = lines[i];
    maxJ           = std::max({maxJ, ln.Ju, ln.Jl});
    // Keep the original lower state even when the angular kernel swaps J roles.
    energies.e0[i] = level_energy(ln.Jl);
  }
  const auto&      l2 = qid.state.at(QuantumNumberType::l2);
  const std::array rats{maxJ, Rational{l2.upper}, Rational{l2.lower}};
  const int        maxL = wigner_init_size(rats);
  prepare_rotational_ladder(energies, maxL, rotational_energy);
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
                                   ConstMatrixView                  dOmega) {
  ARTS_USER_ERROR_IF((dW.npages() != 0 and (dW.nrows() != W.nrows() or dW.ncols() != W.ncols())) or
                         (not dT.empty() and dT.size() != static_cast<Size>(dW.npages())),
                     "Inconsistent Hartmann ECS derivative dimensions")
  const Size n = lines.size();
  if (not n) return;
  validate_isotopologue(bnd_qid.isot);
  const auto& e0 = energies.e0;
  ARTS_USER_ERROR_IF(
      e0.size() != n or dipr.size() != n or W.nrows() != static_cast<Index>(n) or W.ncols() != static_cast<Index>(n),
      "Inconsistent Hartmann ECS kernel dimensions")

  // These are constant for a band
  auto&    l2 = bnd_qid.state.at(QuantumNumberType::l2);
  Rational li = l2.upper;
  Rational lf = l2.lower;

  Rational maxJ{0};
  for (const auto& ln : lines) maxJ = std::max({maxJ, ln.Ju, ln.Jl});

  using std::swap;
  const bool swap_order = li > lf;
  if (swap_order) swap(li, lf);
  if (abs(li - lf) > 1) return;

  const Numeric T = atm.temperature;

  const std::array rats{maxJ, li, lf};
  const int        maxL  = wigner_init_size(rats);
  const auto       basis = prepare_basis(maxL, energies, rovib_data, T0, bnd_qid.isot, broadening_species, atm);
  ARTS_USER_ERROR_IF(
      (not dQ.empty() and (dQ.nrows() != dW.npages() or dQ.ncols() != static_cast<Index>(basis.Q.size()))) or
          (not dOmega.empty() and
           (dOmega.nrows() != dW.npages() or dOmega.ncols() != static_cast<Index>(basis.Omega.size()))),
      "Inconsistent Hartmann ECS basis derivative dimensions")

  coupling_kernel(W, lines, li, lf, swap_order, basis, e0, T, maxL, dW, dT, dQ, dOmega);

  apply_sum_rule(W, dipr, e0, T, dW, dT);
}
}  // namespace lbl::voigt::ecs::hartmann
