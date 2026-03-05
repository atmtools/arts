#include "lbl_lineshape_voigt_ecs_stotop.h"

#include <arts_conversions.h>
#include <atm.h>
#include <wigner_functions.h>

namespace lbl::voigt::ecs::stotop {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

namespace {
Numeric wig3(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) {
  return WIGNER3(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) {
  return WIGNER6(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

/*! Compute rotational energy for a symmetric top molecule
 *
 * E(J, K) = B*J*(J+1) + (A-B)*K^2 - D_J*[J*(J+1)]^2
 *            - D_JK*J*(J+1)*K^2 - D_K*K^4
 *
 * For the ECS basis rates we only need E(L) for the collisional
 * angular momentum transfer channel, where L has no K-dependence
 * (the basis rates are for the atom-like IOS limit).  So we use
 * the simple rigid rotor formula E = B*J*(J+1).
 *
 * Species-specific rotational constants:
 *   NH3-4111 (14NH3):  B0 = 9.9402 cm^{-1}
 *   PH3-1111 (31PH3):  B0 = 4.4522 cm^{-1}
 *
 * @param[in] isot  The isotopologue
 * @return A function J -> E(J) in Joule
 */
std::function<Numeric(Rational)> erot_selection(const SpeciesIsotope& isot) {
  // NH3 main isotopologue (14N-1H3)
  if (isot == "NH3-4111"_isot) {
    return [](const Rational J) -> Numeric {
      return Conversion::kaycm2joule(9.9402) * Numeric(J * (J + 1));
    };
  }

  // PH3 main isotopologue (31P-1H3)
  if (isot == "PH3-1111"_isot) {
    return [](const Rational J) -> Numeric {
      return Conversion::kaycm2joule(4.4522) * Numeric(J * (J + 1));
    };
  }

  ARTS_USER_ERROR("{} has no rotational energies for symmetric top ECS in ARTS",
                  isot.FullName())
  return [](const Rational J) -> Numeric {
    return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN();
  };
}
}  // namespace

Numeric reduced_dipole(const Rational Jf, const Rational Ji, const Rational K) {
  // d(Jf, Ji, K) = (-1)^{K+Jf} sqrt(2*Jf+1) * 3j(Jf, 1, Ji; K, 0, -K)
  // This is identical in structure to Eq. (10) of Rodrigues et al. 1997
  // with K replacing the vibrational angular momentum l.
  if (not iseven(Jf + K + 1))
    return -sqrtr(2 * Jf + 1) *
           wigner3j(Jf, Rational{1}, Ji, K, Rational{0}, -K);
  return +sqrtr(2 * Jf + 1) * wigner3j(Jf, Rational{1}, Ji, K, Rational{0}, -K);
}

void relaxation_matrix_offdiagonal(MatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesEnum broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm) {
  using Conversion::kelvin2joule;

  const Size n = bnd.size();
  if (not n) return;

  // K is a line-level quantum number for symmetric tops.
  // Lines with different K are not coupled (ΔK=0 collisions only).
  // Read K for each line and find the maximum K for Wigner init.
  std::vector<Rational> Kvec(n);
  for (Size i = 0; i < n; i++) {
    Kvec[i] = bnd.lines[sorting[i]].qn.at(QuantumNumberType::K).lower;
  }

  const Numeric T = atm.temperature;

  const auto erot = erot_selection(bnd_qid.isot);

  const std::array rats{bnd.max(QuantumNumberType::J),
                        bnd.max(QuantumNumberType::K)};
  const int maxL = wigner_init_size(rats);

  const auto Om = [&]() {
    Vector out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Omega(atm.temperature,
                                bnd.front().ls.T0,
                                broadening_species == SpeciesEnum::Bath
                                    ? atm.mean_mass()
                                    : atm.mean_mass(broadening_species),
                                bnd_qid.isot.mass,
                                erot(Rational{i}),
                                erot(Rational{i - 2}));
    return out;
  }();

  const auto Q = [&]() {
    Vector out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Q(
          Rational{i}, atm.temperature, bnd.front().ls.T0, erot(Rational{i}));
    return out;
  }();

  // The coupling for symmetric tops with ΔK=0 is identical to the
  // linear molecule (Hartmann) case with K replacing l:
  //
  //   R(J→J'; K) = (-1)^{J+J'+1} (2J'+1) Ω(J)
  //     × Σ_L (2L+1) 3j(J, J', L; K, -K, 0) × 3j(Jf, Jf', L; K, -K, 0)
  //                   × 6j(Ji, Jf, 1; Jf', Ji', L) × Q(L)/Ω(L)
  //
  // For pure rotational ΔK=0 transitions: Ji = J (upper), Jf = J±1 (lower).
  // Ki = Kf = K is constant for the entire sub-band.

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    auto& J           = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J);
    const Rational Ji = J.upper;
    const Rational Jf = J.lower;
    const Rational Ki = Kvec[i];

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;

      // Only couple lines within the same K sub-band
      if (Kvec[j] != Ki) continue;

      auto& J_p           = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J);
      const Rational Ji_p = J_p.upper;
      const Rational Jf_p = J_p.lower;

      // Only compute the downward element (J' ≤ J)
      if (Jf_p > Jf) continue;

      Index L         = std::max(std::abs((Ji - Ji_p).toIndex()),
                                 std::abs((Jf - Jf_p).toIndex()));
      L              += L % 2;
      const Index Lf  = std::min((Ji + Ji_p).toIndex(), (Jf + Jf_p).toIndex());

      Numeric sum = 0;
      for (; L <= Lf; L += 2) {
        const Numeric a  = wig3(Ji, Ji_p, Rational{L}, Ki, -Ki, Rational{0});
        const Numeric b  = wig3(Jf, Jf_p, Rational{L}, Ki, -Ki, Rational{0});
        const Numeric c  = wig6(Ji, Jf, Rational{1}, Jf_p, Ji_p, Rational{L});
        sum             += a * b * c * Numeric(2 * L + 1) * Q[L] / Om[L];
      }
      const Numeric ECS = Om[Ji.toIndex()];
      const Numeric scl =
          ECS * Numeric(2 * Ji_p + 1) * sqrtr((2 * Jf + 1) * (2 * Jf_p + 1));
      sum *= scl;

      // Downward element and detailed balance for upward
      W[j, i] = sum;
      W[i, j] = sum * std::exp((erot(Jf_p) - erot(Jf)) / kelvin2joule(T));
    }
  }
  arts_wigner_thread_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")

  // Sum rule correction (same as Hartmann, but within K sub-bands)
  for (Size i = 0; i < n; i++) {
    Numeric sumlw     = 0.0;
    Numeric sumup     = 0.0;
    const Rational Ki = Kvec[i];

    for (Size j = 0; j < n; j++) {
      if (Kvec[j] != Ki) continue;  // Only within same K sub-band

      if (j > i) {
        sumlw += dipr[j] * W[j, i];
      } else {
        sumup += dipr[j] * W[j, i];
      }
    }

    const Rational Ji = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J).lower;
    for (Size j = i + 1; j < n; j++) {
      if (Kvec[j] != Ki) continue;  // Only within same K sub-band

      const Rational Jj =
          bnd.lines[sorting[j]].qn.at(QuantumNumberType::J).lower;
      if (std::abs(sumlw) <= std::abs(sumup) * 1e-6) {
        // The sum rule correction would amplify off-diagonal elements
        // excessively.  This can happen for Q-branch (ΔJ=0) transitions
        // where the coupling structure differs from P/R branches.
        // Zero the elements to prevent NaN propagation.
        W[j, i] = 0.0;
        W[i, j] = 0.0;
      } else {
        W[j, i] *= -sumup / sumlw;
        W[i, j]  = W[j, i] * std::exp((erot(Ji) - erot(Jj)) / kelvin2joule(T));
      }
    }
  }
}
}  // namespace lbl::voigt::ecs::stotop
