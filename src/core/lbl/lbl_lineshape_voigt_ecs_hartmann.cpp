#include <wigner_functions.h>

#include "atm.h"
#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

namespace lbl::voigt::ecs::hartmann {
std::function<Numeric(Rational)> erot_selection(
    const SpeciesIsotopeRecord& isot) {
  if (isot.spec == Species::Species::CarbonDioxide and isot.isotname == "626") {
    return [](const Rational J) -> Numeric {
      return Conversion::kaycm2joule(0.39021) * Numeric(J * (J + 1));
    };
  }

  ARTS_USER_ERROR(isot.FullName(), " has no rotational energies in ARTS")
  return [](const Rational J) -> Numeric {
    return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN();
  };
}

Numeric reduced_dipole(const Rational Jf,
                       const Rational Ji,
                       const Rational lf,
                       const Rational li,
                       const Rational k = 1) {
  if (not iseven(Jf + lf + 1))
    return -sqrt(2 * Jf + 1) * wigner3j(Jf, k, Ji, li, lf - li, -lf);
  return +sqrt(2 * Jf + 1) * wigner3j(Jf, k, Ji, li, lf - li, -lf);
}

void relaxation_matrix_offdiagonal(ExhaustiveMatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const Species::Species broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm) {
  using Conversion::kelvin2joule;

  const Size n = bnd.size();
  if (not n) return;

  // These are constant for a band
  auto& l2 = bnd_qid.val[QuantumNumberType::l2];
  Rational li = l2.upp();
  Rational lf = l2.low();

  using std::swap;
  const bool swap_order = li > lf;
  if (swap_order) swap(li, lf);
  const int sgn = iseven(li + lf + 1) ? -1 : 1;
  if (abs(li - lf) > 1) return;

  const Numeric T0 = bnd.front().ls.T0;
  const Numeric T = atm.temperature;

  const Numeric mass = broadening_species == SpeciesEnum::Bath
                           ? atm.mean_mass()
                           : atm.mean_mass(broadening_species);

  const auto erot = erot_selection(bnd_qid.Isotopologue());

  for (Size i = 0; i < n; i++) {
    auto& J = bnd.lines[sorting[i]].qn.val[QuantumNumberType::J];
    Rational Ji = J.upp();
    Rational Jf = J.low();
    if (swap_order) swap(Ji, Jf);

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;
      auto& J_p = bnd.lines[sorting[j]].qn.val[QuantumNumberType::J];
      Rational Ji_p = J_p.upp();
      Rational Jf_p = J_p.low();
      if (swap_order) swap(Ji_p, Jf_p);

      // Select upper quantum number
      if (Jf_p > Jf) continue;

      Index L = std::max(std::abs((Ji - Ji_p).toIndex()),
                         std::abs((Jf - Jf_p).toIndex()));
      L += L % 2;
      const Index Lf = std::min((Ji + Ji_p).toIndex(), (Jf + Jf_p).toIndex());

      Numeric sum = 0;
      for (; L <= Lf; L += 2) {
        const Numeric a = wigner3j(Ji_p, L, Ji, li, 0, -li);
        const Numeric b = wigner3j(Jf_p, L, Jf, lf, 0, -lf);
        const Numeric c = wigner6j(Ji, Jf, 1, Jf_p, Ji_p, L);
        const Numeric QL = rovib_data.Q(L, T, T0, erot(L));
        const Numeric ECS = rovib_data.Omega(
            T, T0, mass, bnd_qid.Isotopologue().mass, erot(L), erot(L - 2));
        sum += a * b * c * Numeric(2 * L + 1) * QL / ECS;
      }
      const Numeric ECS = rovib_data.Omega(
          T, T0, mass, bnd_qid.Isotopologue().mass, erot(Ji), erot(Ji - 2));
      const Numeric scl = sgn * ECS * Numeric(2 * Ji_p + 1) *
                          sqrt((2 * Jf + 1) * (2 * Jf_p + 1));
      sum *= scl;

      // Add to W and rescale to upwards element by the populations
      W(j, i) = sum;
      W(i, j) = sum * std::exp((erot(Jf_p) - erot(Jf)) / kelvin2joule(T));
    }
  }

  // Undocumented negative absolute sign
  for (Size i = 0; i < n; i++)
    for (Size j = 0; j < n; j++)
      if (j not_eq i and W(i, j) > 0) W(i, j) *= -1;

  // Sum rule correction
  for (Size i = 0; i < n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;

    for (Size j = 0; j < n; j++) {
      if (j > i) {
        sumlw += std::abs(dipr[j]) * W(j, i);  // Undocumented abs-sign
      } else {
        sumup += std::abs(dipr[j]) * W(j, i);  // Undocumented abs-sign
      }
    }

    const Rational Ji =
        bnd.lines[sorting[i]].qn.val[QuantumNumberType::J].low();
    for (Size j = i + 1; j < n; j++) {
      const Rational Jj =
          bnd.lines[sorting[j]].qn.val[QuantumNumberType::J].low();
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= -sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(Ji) - erot(Jj)) /
                                     kelvin2joule(T));  // This gives LTE
      }
    }
  }
}
}  // namespace lbl::voigt::ecs::hartmann
