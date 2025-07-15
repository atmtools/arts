#include "lbl_lineshape_voigt_ecs_makarov.h"

#include <wigner_functions.h>

#include "debug.h"

namespace lbl::voigt::ecs::makarov {
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
}  // namespace

Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N) {
  return (iseven(Jl + N) ? 1 : -1) * sqrt(6 * (2 * Jl + 1) * (2 * Ju + 1)) *
         wigner6j(1, 1, 1, Jl, Ju, N);
};

namespace {
/*! Compute the rotational energy of ground-state O2 at N and J
 * 
 * If the template argument evaluates true, the erot<false>(1, 0)
 * energy is removed from the output of erot<false>(N, J).
 * 
 * @param[in] N Main rotational number
 * @param[in] j Main rotational number plus spin (if j < 0 then J=N)
 * @return Rotational energy in Joule
 */
template <bool rescale_pure_rotational = true>
constexpr Numeric erot(const Rational N, const Rational j = -1) try {
  const Rational J = j < 0 ? N : j;

  if constexpr (rescale_pure_rotational) {
    return erot<false>(N, J) - erot<false>(1, 0);
  } else {
    using Conversion::mhz2joule;
    using Math::pow2;
    using Math::pow3;

    constexpr Numeric B0  = 43100.4425e0;
    constexpr Numeric D0  = .145123e0;
    constexpr Numeric H0  = 3.8e-08;
    constexpr Numeric xl0 = 59501.3435e0;
    constexpr Numeric xg0 = -252.58633e0;
    constexpr Numeric xl1 = 0.058369e0;
    constexpr Numeric xl2 = 2.899e-07;
    constexpr Numeric xg1 = -2.4344e-04;
    constexpr Numeric xg2 = -1.45e-09;

    const Numeric XN      = static_cast<Numeric>(N);
    const Numeric XX      = XN * (XN + 1);
    const Numeric xlambda = xl0 + xl1 * XX + xl2 * pow2(XX);
    const Numeric xgama   = xg0 + xg1 * XX + xg2 * pow2(XX);
    const Numeric C1      = B0 * XX - D0 * pow2(XX) + H0 * pow3(XX);

    if (J < N) {
      if (N == 1)  // erot<false>(1, 0)
        return mhz2joule(C1 - (xlambda + B0 * (2. * XN - 1.) + xgama * XN));
      return mhz2joule(C1 - (xlambda + B0 * (2. * XN - 1.) + xgama * XN) +
                       std::sqrt(pow2(B0 * (2. * XN - 1.)) + pow2(xlambda) -
                                 2. * B0 * xlambda));
    }
    if (J > N)
      return mhz2joule(C1 -
                       (xlambda - B0 * (2. * XN + 3.) - xgama * (XN + 1.)) -
                       std::sqrt(pow2(B0 * (2. * XN + 3.)) + pow2(xlambda) -
                                 2. * B0 * xlambda));
    return mhz2joule(C1);
  }
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void relaxation_matrix_offdiagonal(MatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesEnum broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm) try {
  using Conversion::kelvin2joule;

  ARTS_USER_ERROR_IF(bnd_qid.Isotopologue() != "O2-66"_isot,
                     "Bad isotopologue: {}",
                     bnd_qid.Isotopologue())

  if (bnd.size() == 0) return;

  const auto bk = [](const Rational& r) -> Numeric { return sqrt(2 * r + 1); };

  const auto n = bnd.size();

  auto& S           = bnd_qid.val[QuantumNumberType::S];
  const Rational Si = S.upp();
  const Rational Sf = S.low();

  const auto maxL = temp_init_size(bnd.max(QuantumNumberType::J),
                                   bnd.max(QuantumNumberType::N));

  const auto Om = [&]() {
    Vector out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Omega(atm.temperature,
                                bnd.front().ls.T0,
                                broadening_species == SpeciesEnum::Bath
                                    ? atm.mean_mass()
                                    : atm.mean_mass(broadening_species),
                                bnd_qid.Isotopologue().mass,
                                erot(i),
                                erot(i - 2));
    return out;
  }();

  const auto Q = [&]() {
    Vector out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Q(i, atm.temperature, bnd.front().ls.T0, erot(i));
    return out;
  }();

  arts_wigner_thread_init(maxL);
  for (Size i = 0; i < n; i++) {
    auto& J = bnd.lines[sorting[i]].qn.val[QuantumNumberType::J];
    auto& N = bnd.lines[sorting[i]].qn.val[QuantumNumberType::N];

    const Rational Ji = J.upp();
    const Rational Jf = J.low();
    const Rational Ni = N.upp();
    const Rational Nf = N.low();

    for (Size j = 0; j < n; j++) {
      if (i == j) continue;

      auto& J_p = bnd.lines[sorting[j]].qn.val[QuantumNumberType::J];
      auto& N_p = bnd.lines[sorting[j]].qn.val[QuantumNumberType::N];

      const Rational Ji_p = J_p.upp();
      const Rational Jf_p = J_p.low();
      const Rational Ni_p = N_p.upp();
      const Rational Nf_p = N_p.low();

      if (Jf_p > Jf) continue;

      // Tran etal 2006 symbol with modifications:
      //    1) [Ji] * [Ji_p] instead of [Ji_p] ^ 2 in partial accordance with Makarov etal 2013
      Numeric sum       = 0;
      const Numeric scl = (iseven(Ji_p + Ji + 1) ? 1 : -1) * bk(Ni) * bk(Nf) *
                          bk(Nf_p) * bk(Ni_p) * bk(Jf) * bk(Jf_p) * bk(Ji) *
                          bk(Ji_p);
      const auto [L0, L1] =
          wigner_limits(wigner3j_limits<3>(Ni_p, Ni),
                        {Rational(2), std::numeric_limits<Index>::max()});
      for (Rational L = L0; L <= L1; L += 2) {
        const Numeric a = wig3(Ni_p, Ni, L, 0, 0, 0);
        const Numeric b = wig3(Nf_p, Nf, L, 0, 0, 0);
        const Numeric c = wig6(L, Ji, Ji_p, Si, Ni_p, Ni);
        const Numeric d = wig6(L, Jf, Jf_p, Sf, Nf_p, Nf);
        const Numeric e = wig6(L, Ji, Ji_p, 1, Jf_p, Jf);
        sum += a * b * c * d * e * Numeric(2 * L + 1) * Q[L.toIndex()] /
               Om[L.toIndex()];
      }
      sum *= scl * Om[Ni.toIndex()];

      // Add to W and rescale to upwards element by the populations
      W[i, j] = sum;
      W[j, i] =
          sum * std::exp((bnd.lines[sorting[j]].e0 - bnd.lines[sorting[i]].e0) /
                         kelvin2joule(atm.temperature));
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

    for (Size j = i + 1; j < n; j++) {
      if (sumlw == 0) {
        W[j, i] = 0.0;
        W[i, j] = 0.0;
      } else {
        W[j, i] *= -sumup / sumlw;
        W[i, j] =
            W[j, i] *
            std::exp((bnd.lines[sorting[i]].e0 - bnd.lines[sorting[j]].e0) /
                     kelvin2joule(atm.temperature));
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs::makarov