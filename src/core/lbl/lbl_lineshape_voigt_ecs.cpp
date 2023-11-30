#include "lbl_lineshape_voigt_ecs.h"

#include <new_jacobian.h>
#include <physics_funcs.h>

#include <Faddeeva.hh>
#include <algorithm>

#include "configtypes.h"
#include "debug.h"
#include "lbl_lineshape_linemixing.h"
#include "partfun.h"
#include "sorting.h"
#include "wigner_functions.h"

#undef WIGNER3
#undef WIGNER6

namespace lbl::voigt::ecs {
#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

Numeric wig3(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) noexcept {
  return WIGNER3(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) noexcept {
  return WIGNER6(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

ComputeData::ComputeData(const ExhaustiveConstVectorView& f_grid,
                         const AtmPoint& atm,
                         const Vector2& los,
                         const zeeman::pol pol)
    : scl(f_grid.size()), shape(f_grid.size()) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 scl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -N * f * std::expm1(-r) * c;
                 });

  update_zeeman(los, atm.mag, pol);
}

void ComputeData::update_zeeman(const Vector2& los,
                                const Vector3& mag,
                                const zeeman::pol pol) {
  npm = zeeman::norm_view(pol, mag, los);
}

void ComputeData::core_calc_eqv() {
  /* FIXME:  (Added 2021-01-19; Richard Larsson)
    * 
    * This function cannot easily be used with partial derivatives
    * 
    * Doing so would allow the entire ECS approach to be analytical
    * in its derivatives.
    * 
    * The problem in short:
    *    There is an Eigenvalue decomposition happening (W = V diag(e) V^-1)
    * 
    *    We use V and e in a strange way.  We do not know how to get the
    *    partial derivatives of these two variables
    * 
    * The question:
    *    How do we compute dV and de?  We can get dW somewhat easily.
    */

  const auto n = pop.size();
  const auto m = Ws.size();

  for (Size k = 0; k < m; k++) {
    auto& V = Vs[k];
    auto& W = Ws[k];
    auto& eqv_str = eqv_strs[k];
    auto& eqv_val = eqv_vals[k];

    diagonalize(V, eqv_val, W);
    eqv_val += freq_offset;

    // Do the matrix forward multiplication
    for (Index i = 0; i < n; i++) {
      for (Index j = 0; j < n; j++) {
        eqv_str[i] += dip[j] * V(j, i);
      }
    }

    // Do the matrix backward multiplication
    inv(V, V);
    for (Index i = 0; i < n; i++) {
      Complex z(0, 0);
      for (Index j = 0; j < n; j++) {
        z += pop[j] * dip[j] * V(i, j);
      }
      eqv_str[i] *= z;
    }
  }
}

void ComputeData::core_calc(const ExhaustiveConstVectorView& f_grid) {
  core_calc_eqv();

  const auto n = f_grid.size();
  shape = 0;

  for (Size k = 0; k < Ws.size(); k++) {
    for (Index i = 0; i < n; i++) {
      const Numeric gamd = gd_fac * eqv_vals[k][i].real();
      const Numeric cte = Constant::sqrt_ln_2 / gamd;
      for (Index iv = 0; iv < f_grid.nelem(); iv++) {
        const Complex z = (eqv_vals[k][i] - f_grid[iv]) * cte;
        const Complex w = Faddeeva::w(z);
        shape[i] += vmrs[k] * eqv_strs[k][i] * w / gamd;
      }
    }
  }

  // for (Index i = 0; i < n; i++) {
  //   if (shape[i].real() < 0) shape[i] = 0;
  // }
}

/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N) {
  return (iseven(Jl + N) ? 1 : -1) * sqrt(6 * (2 * Jl + 1) * (2 * Ju + 1)) *
         wigner6j(1, 1, 1, Jl, Ju, N);
};

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
constexpr Numeric erot(const Rational N, const Rational j = -1) {
  const Rational J = j < 0 ? N : j;

  if constexpr (rescale_pure_rotational) {
    return erot<false>(N, J) - erot<false>(1, 0);
  } else {
    using Conversion::mhz2joule;
    using Math::pow2;
    using Math::pow3;

    constexpr Numeric B0 = 43100.4425e0;
    constexpr Numeric D0 = .145123e0;
    constexpr Numeric H0 = 3.8e-08;
    constexpr Numeric xl0 = 59501.3435e0;
    constexpr Numeric xg0 = -252.58633e0;
    constexpr Numeric xl1 = 0.058369e0;
    constexpr Numeric xl2 = 2.899e-07;
    constexpr Numeric xg1 = -2.4344e-04;
    constexpr Numeric xg2 = -1.45e-09;

    const auto XN = Numeric(N);
    const Numeric XX = XN * (XN + 1);
    const Numeric xlambda = xl0 + xl1 * XX + xl2 * pow2(XX);
    const Numeric xgama = xg0 + xg1 * XX + xg2 * pow2(XX);
    const Numeric C1 = B0 * XX - D0 * pow2(XX) + H0 * pow3(XX);

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

void relaxation_matrix_offdiagonal(Matrix& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const linemixing::species_data& rovib_data,
                                   const Numeric T) {
  using Conversion::kelvin2joule;

  const auto bk = [](const Rational& r) -> Numeric { return sqrt(2 * r + 1); };

  const auto n = bnd.size();

  auto& S = bnd_qid.val[QuantumNumberType::S];
  const Rational Si = S.upp();
  const Rational Sf = S.low();

  Vector dipr(n);
  for (Size i = 0; i < n; i++) {
    auto& J = bnd.lines[sorting[i]].qn.val[QuantumNumberType::J];
    auto& N = bnd.lines[sorting[i]].qn.val[QuantumNumberType::N];

    dipr[i] = reduced_dipole(J.upp(), J.low(), N.upp());
  }

  const auto maxL = temp_init_size(bnd.max(QuantumNumberType::J),
                                   bnd.max(QuantumNumberType::N));

  const auto Om = [&]() {
    std::vector<Numeric> out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Omega(T,
                                bnd.front().ls.T0,
                                bnd_qid.Isotopologue().mass,
                                erot(i),
                                erot(i - 2));
    return out;
  }();

  const auto Q = [&]() {
    std::vector<Numeric> out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Q(i, T, bnd.front().ls.T0, erot(i));
    return out;
  }();

  wig_thread_temp_init(maxL);
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
      Numeric sum = 0;
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
      W(i, j) = sum;
      W(j, i) = sum * std::exp((erot(Nf_p) - erot(Nf)) / kelvin2joule(T));
    }
  }
  wig_temp_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")

  // Sum rule correction
  for (Size i = 0; i < n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;

    for (Size j = 0; j < n; j++) {
      if (j > i) {
        sumlw += dipr[j] * W(j, i);
      } else {
        sumup += dipr[j] * W(j, i);
      }
    }

    const Rational Ni =
        bnd.lines[sorting[i]].qn.val[QuantumNumberType::N].low();
    for (Size j = i + 1; j < n; j++) {
      const Rational Nj =
          bnd.lines[sorting[j]].qn.val[QuantumNumberType::N].low();
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= -sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(Ni) - erot(Nj)) /
                                     kelvin2joule(T));  // This gives LTE
      }
    }
  }
}

void ComputeData::adapt(const QuantumIdentifier& bnd_qid,
                        const band_data& bnd,
                        const linemixing::species_data_map& rovib_data,
                        const AtmPoint& atm) {
  const auto n = bnd.size();

  pop.resize(n);
  dip.resize(n);
  sort.resize(n);
  Wimag.resize(n, n);

  vmrs.resize(1);
  eqv_strs.resize(1);
  eqv_strs[0].resize(n);
  eqv_vals.resize(1);
  eqv_vals[0].resize(n);
  Ws.resize(1);
  Ws[0].resize(n, n);
  Vs.resize(1);
  Vs[0].resize(n, n);

  Ws[0] = 0;
  vmrs[0] = 1;
  eqv_strs[0] = 0;
  eqv_vals[0] = 0;

  const Numeric T = atm.temperature;
  const Numeric QT = PartitionFunctions::Q(T, bnd_qid.Isotopologue());

  for (Size i = 0; i < n; i++) {
    const auto& line = bnd.lines[i];
    const Numeric s = line.s(T, QT);

    //! NOTE: Missing factor of (c^2 / 8pi) in pop[i], important??
    pop[i] = Math::pow3(line.f0) * s / line.a;
    dip[i] = std::sqrt(s / pop[i]);

    //! FIXME: SIGNS????
  }

  //! Must remember the sorting for the quantum numbers
  std::iota(sort.begin(), sort.end(), 0);
  bubble_sort_by(
      [&](Size i, Size j) {
        const auto a = bnd.lines[i].f0 * pop[i] * dip[i] * dip[i];
        const auto b = bnd.lines[j].f0 * pop[j] * dip[j] * dip[j];
        return b > a;
      },
      sort,
      pop,
      dip);

  const auto& ls = bnd.front().ls.single_models;
  Numeric vmr = 0;
  for (Size i = 0; i < bnd.front().ls.single_models.size(); i++) {
    const auto& rovib_data_it = rovib_data.find(ls[i].species);
    ARTS_USER_ERROR_IF(rovib_data_it == rovib_data.end(),
                       "No rovib data for species " , ls[i].species)
    const Numeric this_vmr =
        ls[i].species == Species::Species::Bath ? 1 - vmr : atm[ls[i].species];
    vmr += this_vmr;
    relaxation_matrix_offdiagonal(
        Wimag, bnd_qid, bnd, sort, rovib_data_it->second, atm.temperature);
    for (Size ir = 0; ir < n; ir++) {
      for (Size ic = 0; ic < n; ic++) {
        Ws[0](ir, ic) += 1i * this_vmr * Wimag(ir, ic);
      }
    }
  }
}

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true>,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const linemixing::species_data_map& rovib_data,
               const AtmPoint& atm,
               const zeeman::pol pol) {
  if (pol != zeeman::pol::no) {
    for (auto& line : bnd) {
      ARTS_USER_ERROR_IF(
          line.z.active(),
          "Zeeman effect and ECS in combination is not yet possible.")
    }
  }

  ARTS_USER_ERROR_IF(jacobian_targets.target_count() > 0,
                     "No Jacobian support.")

  ARTS_USER_ERROR_IF(bnd_qid.Isotopologue() != "O2-66", "Bad isotopologue.")

  if (bnd.size() == 0) return;

  const bool one_by_one = bnd.front().ls.one_by_one;

  if (one_by_one) {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            bnd, [](const auto& line) { return not line.ls.one_by_one; }),
        "One-by-one lineshape not supported for ECS.");

    ARTS_USER_ERROR("NOT IMPLEMENTED")
  } else {
    com_data.adapt(bnd_qid, bnd, rovib_data, atm);
  }

  for (Index i = 0; i < f_grid.size(); ++i) {
    pm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.shape[i]);
  }
}
}  // namespace lbl::voigt::ecs