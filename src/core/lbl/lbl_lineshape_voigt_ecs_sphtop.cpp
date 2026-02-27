#include "lbl_lineshape_voigt_ecs_sphtop.h"

#include <arts_conversions.h>
#include <atm.h>
#include <wigner_functions.h>

namespace lbl::voigt::ecs::sphtop {
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

/*! Compute rotational energy for a spherical top molecule
 *
 * E(J) = B*J*(J+1)
 *
 * Species-specific rotational constants:
 *   CH4-211 (12CH4):  B0 = 5.2410 cm^{-1}
 *
 * @param[in] isot  The isotopologue
 * @return A function J -> E(J) in Joule
 */
std::function<Numeric(Rational)> erot_selection(const SpeciesIsotope& isot) {
  // CH4 main isotopologue (12C-1H4)
  if (isot.spec == SpeciesEnum::Methane and isot.isotname == "211") {
    return [](const Rational J) -> Numeric {
      return Conversion::kaycm2joule(5.2410) * Numeric(J * (J + 1));
    };
  }

  ARTS_USER_ERROR(
      "{} has no rotational energies for spherical top ECS in ARTS",
      isot.FullName())
  return [](const Rational J) -> Numeric {
    return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN();
  };
}
}  // namespace

Numeric reduced_dipole(const Rational Jf, const Rational Ji) {
  // d(Jf, Ji) = (-1)^{Jf+1} sqrt(2*Jf+1) * 3j(Jf, 1, Ji; 0, 0, 0)
  // This is the l=0 limit of the Hartmann formula.
  if (not iseven(Jf + 1))
    return -sqrtr(2 * Jf + 1) *
           wigner3j(Jf, Rational{1}, Ji, Rational{0}, Rational{0}, Rational{0});
  return +sqrtr(2 * Jf + 1) *
         wigner3j(Jf, Rational{1}, Ji, Rational{0}, Rational{0}, Rational{0});
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

  const Numeric T = atm.temperature;

  const auto erot = erot_selection(bnd_qid.isot);

  const std::array rats{bnd.max(QuantumNumberType::J)};
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

  // The coupling for spherical tops with l=0 is:
  //
  //   R(J→J'; 0) = (-1)^{J+J'+1} (2J'+1) Ω(J)
  //     × Σ_L (2L+1) 3j(J, J', L; 0, 0, 0) × 3j(Jf, Jf', L; 0, 0, 0)
  //                   × 6j(Ji, Jf, 1; Jf', Ji', L) × Q(L)/Ω(L)
  //
  // This is identical to the Hartmann (linear molecule) formula with
  // l_i = l_f = 0.

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    auto& J = bnd.lines[sorting[i]].qn.at(QuantumNumberType::J);
    const Rational Ji = J.upper;
    const Rational Jf = J.lower;

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;

      auto& J_p = bnd.lines[sorting[j]].qn.at(QuantumNumberType::J);
      const Rational Ji_p = J_p.upper;
      const Rational Jf_p = J_p.lower;

      // Only compute the downward element (J' ≤ J)
      if (Jf_p > Jf) continue;

      Index L  = std::max(std::abs((Ji - Ji_p).toIndex()),
                         std::abs((Jf - Jf_p).toIndex()));
      L       += L % 2;
      const Index Lf =
          std::min((Ji + Ji_p).toIndex(), (Jf + Jf_p).toIndex());

      Numeric sum = 0;
      for (; L <= Lf; L += 2) {
        const Numeric a =
            wig3(Ji, Ji_p, Rational{L}, Rational{0}, Rational{0}, Rational{0});
        const Numeric b =
            wig3(Jf, Jf_p, Rational{L}, Rational{0}, Rational{0}, Rational{0});
        const Numeric c = wig6(Ji, Jf, Rational{1}, Jf_p, Ji_p, Rational{L});
        sum += a * b * c * Numeric(2 * L + 1) * Q[L] / Om[L];
      }
      const Numeric ECS = Om[Ji.toIndex()];
      const Numeric scl = ECS * Numeric(2 * Ji_p + 1) *
                          sqrtr((2 * Jf + 1) * (2 * Jf_p + 1));
      sum *= scl;

      // Downward element and detailed balance for upward
      W[j, i] = sum;
      W[i, j] = sum * std::exp((erot(Jf_p) - erot(Jf)) / kelvin2joule(T));
    }
  }
  arts_wigner_thread_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")

  // Sum rule correction
  for (Size i = 0; i < n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;

    for (Size j = 0; j < n; j++) {
      if (j > i) {
        sumlw += dipr[j] * W[j, i];
      } else {
        sumup += dipr[j] * W[j, i];
      }
    }

    const Rational Ji =
        bnd.lines[sorting[i]].qn.at(QuantumNumberType::J).lower;
    for (Size j = i + 1; j < n; j++) {
      const Rational Jj =
          bnd.lines[sorting[j]].qn.at(QuantumNumberType::J).lower;
      if (std::abs(sumlw) <= std::abs(sumup) * 1e-6) {
        W[j, i] = 0.0;
        W[i, j] = 0.0;
      } else {
        W[j, i] *= -sumup / sumlw;
        W[i, j] =
            W[j, i] * std::exp((erot(Ji) - erot(Jj)) / kelvin2joule(T));
      }
    }
  }
}
}  // namespace lbl::voigt::ecs::sphtop
