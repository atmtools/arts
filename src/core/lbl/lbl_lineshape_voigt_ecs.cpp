#include "lbl_lineshape_voigt_ecs.h"

#include <jacobian.h>
#include <physics_funcs.h>

#include <Faddeeva.hh>
#include <algorithm>
#include <cmath>

#include "arts_omp.h"
#include "atm.h"
#include "configtypes.h"
#include "debug.h"
#include "isotopologues.h"
#include "lbl_lineshape_linemixing.h"
#include "partfun.h"
#include "sorting.h"
#include "species.h"

#undef WIGNER3
#undef WIGNER6

namespace lbl::voigt::ecs {
namespace makarov {
/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N);

void relaxation_matrix_offdiagonal(ExhaustiveMatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const Species::Species broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm);
}  // namespace makarov

namespace hartmann {
Numeric reduced_dipole(const Rational Jf,
                       const Rational Ji,
                       const Rational lf,
                       const Rational li,
                       const Rational k = 1);

void relaxation_matrix_offdiagonal(ExhaustiveMatrixView& W,
                                   const QuantumIdentifier& bnd_qid,
                                   const band_data& bnd,
                                   const ArrayOfIndex& sorting,
                                   const Species::Species broadening_species,
                                   const linemixing::species_data& rovib_data,
                                   const Vector& dipr,
                                   const AtmPoint& atm);
}  // namespace hartmann

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
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -N * f * std::expm1(-r);
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
  const auto m = vmrs.size();

  for (Index k = 0; k < m; k++) {
    auto V = Vs[k];
    auto W = Ws[k];
    auto eqv_str = eqv_strs[k];
    auto eqv_val = eqv_vals[k];

    diagonalize(V, eqv_val, W);

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

  const auto m = vmrs.size();
  const auto n = f_grid.size();
  shape = 0;

  for (Index k = 0; k < m; k++) {
    for (Index i = 0; i < eqv_strs[k].size(); i++) {
      const Numeric gamd = gd_fac * eqv_vals[k][i].real();
      const Numeric cte = Constant::sqrt_ln_2 / gamd;
      for (Index iv = 0; iv < n; iv++) {
        const Complex z = (eqv_vals[k][i] - f_grid[iv]) * cte;
        const Complex w = Faddeeva::w(z);
        shape[iv] += vmrs[k] * eqv_strs[k][i] * w / gamd;
      }
    }
  }
}

void ComputeData::adapt_multi(const QuantumIdentifier& bnd_qid,
                              const band_data& bnd,
                              const linemixing::species_data_map& rovib_data,
                              const AtmPoint& atm,
                              const bool presorted) {
  const auto n = bnd.size();
  const auto m = bnd.front().ls.single_models.size();

  ARTS_USER_ERROR_IF(m == 0, "No broadening species in the band")

  pop.resize(n);
  dip.resize(n);
  dipr.resize(n);
  sort.resize(n);
  Wimag.resize(n, n);

  vmrs.resize(m);
  eqv_strs.resize(m, n);
  eqv_vals.resize(m, n);
  Ws.resize(m, n, n);
  Vs.resize(m, n, n);

  Ws = 0;
  vmrs = 0;
  eqv_strs = 0;
  eqv_vals = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared *
                     atm.temperature / bnd_qid.Isotopologue().mass);

  const Numeric T = atm.temperature;
  const Numeric QT = PartitionFunctions::Q(T, bnd_qid.Isotopologue());

  for (Size i = 0; i < n; i++) {
    const auto& line = bnd.lines[i];
    pop[i] = line.gu * exp(-line.e0 / (Constant::k * T)) / QT;
    dip[i] = 0.5 * Constant::c *
             std::sqrt(line.a / (Math::pow3(line.f0) * Constant::two_pi));

    if (bnd.lineshape == Lineshape::VP_ECS_MAKAROV) {
      auto& J = bnd.lines[i].qn.val[QuantumNumberType::J];
      auto& N = bnd.lines[i].qn.val[QuantumNumberType::N];
      dipr[i] = makarov::reduced_dipole(J.upp(), J.low(), N.upp());
      dip[i] *= std::signbit(dipr[i]) ? -1 : 1;
    } else if (bnd.lineshape == Lineshape::VP_ECS_HARTMANN) {
      auto& J = bnd.lines[i].qn.val[QuantumNumberType::J];
      auto& l2 = bnd_qid.val[QuantumNumberType::l2];
      dipr[i] = hartmann::reduced_dipole(J.upp(), J.low(), l2.upp(), l2.low());
      dip[i] *= std::signbit(dipr[i]) ? -1 : 1;
    }
  }

  //! Must remember the sorting for the quantum numbers
  if (not presorted) {
    std::iota(sort.begin(), sort.end(), 0);
    bubble_sort_by(
        [&](Size i, Size j) {
          const auto a = bnd.lines[sort[i]].f0 * pop[i] * dip[i] * dip[i];
          const auto b = bnd.lines[sort[j]].f0 * pop[j] * dip[j] * dip[j];
          return b > a;
        },
        sort,
        pop,
        dip,
        dipr);
  } else {
    const auto presorter = [this](const auto& vec) {
      auto out = vec;
      for (Size i : sort) {
        out[i] = vec[sort[i]];
      }
      return out;
    };

    pop = presorter(pop);
    dip = presorter(dip);
    dipr = presorter(dipr);
  }

  for (auto& line : bnd) {
    const bool lines_have_same_species =
        std::transform_reduce(line.ls.single_models.begin(),
                              line.ls.single_models.end(),
                              bnd.lines[0].ls.single_models.begin(),
                              true,
                              std::logical_or<>(),

                              [](const auto& lsl, const auto& lsr) {
                                return lsl.species == lsr.species;
                              });

    ARTS_USER_ERROR_IF(not lines_have_same_species, "Bad species combination")
  }

  for (Size i = 0; i < m; i++) {
    const auto spec = bnd.front().ls.single_models[i].species;

    const auto& rovib_data_it = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(
        rovib_data_it == rovib_data.end(), "No rovib data for species ", spec)

    vmrs[i] = spec == Species::Species::Bath ? 1 - sum(vmrs) : atm[spec];

    Wimag = 0.0;
    for (Size k = 0; k < n; k++) {
      Wimag(k, k) = bnd.lines[sort[k]].ls.single_models[i].G0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
      Ws[i](k, k) = bnd.lines[sort[k]].ls.single_models[i].D0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
    }

    if (bnd.lineshape == Lineshape::VP_ECS_MAKAROV) {
      makarov::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else if (bnd.lineshape == Lineshape::VP_ECS_HARTMANN) {
      hartmann::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else {
      ARTS_USER_ERROR("UNKNOWN ECS LINE SHAPE ", bnd.lineshape)
    }

    Ws[i].imag() = Wimag;
  }

  for (Size i = 0; i < n; i++) {
    Ws(joker, i, i) += bnd.lines[sort[i]].f0;
  }
}

void ComputeData::adapt_single(const QuantumIdentifier& bnd_qid,
                               const band_data& bnd,
                               const linemixing::species_data_map& rovib_data,
                               const AtmPoint& atm,
                               const bool presorted) {
  const auto n = bnd.size();

  pop.resize(n);
  dip.resize(n);
  dipr.resize(n);
  sort.resize(n);
  Wimag.resize(n, n);

  vmrs.resize(1);
  eqv_strs.resize(1, n);
  eqv_vals.resize(1, n);
  Ws.resize(1, n, n);
  Vs.resize(1, n, n);

  Ws = 0;
  vmrs = 1;
  eqv_strs = 0;
  eqv_vals = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared *
                     atm.temperature / bnd_qid.Isotopologue().mass);

  const Numeric T = atm.temperature;
  const Numeric QT = PartitionFunctions::Q(T, bnd_qid.Isotopologue());

  for (Size i = 0; i < n; i++) {
    const auto& line = bnd.lines[i];
    pop[i] = line.gu * exp(-line.e0 / (Constant::k * T)) / QT;
    dip[i] = 0.5 * Constant::c *
             std::sqrt(line.a / (Math::pow3(line.f0) * Constant::two_pi));

    if (bnd.lineshape == Lineshape::VP_ECS_MAKAROV) {
      auto& J = bnd.lines[i].qn.val[QuantumNumberType::J];
      auto& N = bnd.lines[i].qn.val[QuantumNumberType::N];
      dipr[i] = makarov::reduced_dipole(J.upp(), J.low(), N.upp());
      dip[i] *= std::signbit(dipr[i]) ? -1 : 1;
    } else if (bnd.lineshape == Lineshape::VP_ECS_HARTMANN) {
      auto& J = bnd.lines[i].qn.val[QuantumNumberType::J];
      auto& l2 = bnd_qid.val[QuantumNumberType::l2];
      dipr[i] = hartmann::reduced_dipole(J.upp(), J.low(), l2.upp(), l2.low());
      dip[i] *= std::signbit(dipr[i]) ? -1 : 1;
    }
  }

  //! Must remember the sorting for the quantum numbers
  if (not presorted) {
    std::iota(sort.begin(), sort.end(), 0);
    bubble_sort_by(
        [&](Size i, Size j) {
          const auto a = bnd.lines[sort[i]].f0 * pop[i] * dip[i] * dip[i];
          const auto b = bnd.lines[sort[j]].f0 * pop[j] * dip[j] * dip[j];
          return b > a;
        },
        sort,
        pop,
        dip,
        dipr);
  } else {
    const auto presorter = [this](const auto& vec) {
      auto out = vec;
      for (Size i : sort) {
        out[i] = vec[sort[i]];
      }
      return out;
    };

    pop = presorter(pop);
    dip = presorter(dip);
    dipr = presorter(dipr);
  }

  for (auto& line : bnd) {
    const bool lines_have_same_species =
        std::transform_reduce(line.ls.single_models.begin(),
                              line.ls.single_models.end(),
                              bnd.lines[0].ls.single_models.begin(),
                              true,
                              std::logical_or<>(),

                              [](const auto& lsl, const auto& lsr) {
                                return lsl.species == lsr.species;
                              });

    ARTS_USER_ERROR_IF(not lines_have_same_species, "Bad species combination")
  }

  Numeric vmr = 0;
  for (Size i = 0; i < bnd.front().ls.single_models.size(); i++) {
    const auto spec = bnd.front().ls.single_models[i].species;

    const auto& rovib_data_it = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(
        rovib_data_it == rovib_data.end(), "No rovib data for species ", spec)

    const Numeric this_vmr =
        spec == Species::Species::Bath ? 1 - vmr : atm[spec];

    Wimag = 0.0;
    for (Size k = 0; k < n; k++) {
      Wimag(k, k) = bnd.lines[sort[k]].ls.single_models[i].G0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
    }

    if (bnd.lineshape == Lineshape::VP_ECS_MAKAROV) {
      makarov::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else if (bnd.lineshape == Lineshape::VP_ECS_HARTMANN) {
      hartmann::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else {
      ARTS_USER_ERROR("UNKNOWN ECS LINE SHAPE ", bnd.lineshape)
    }

    for (Size ir = 0; ir < n; ir++) {
      for (Size ic = 0; ic < n; ic++) {
        imag_val(Ws[0](ir, ic)) += this_vmr * Wimag(ir, ic);
      }
    }

    vmr += this_vmr;
  }

  if (vmr != 1.0) Ws[0] /= vmr;

  for (Size i = 0; i < n; i++) {
    real_val(Ws[0](i, i)) =
        bnd.lines[sort[i]].f0 + bnd.lines[sort[i]].ls.D0(atm);
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
               const zeeman::pol pol,
               const bool no_negative_absorption) {
  if (pol != zeeman::pol::no) {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            bnd, [](auto& zee) { return zee.on; }, &line::z),
        "Zeeman effect and ECS in combination is not yet possible.")
    return;
  }

  ARTS_USER_ERROR_IF(jacobian_targets.target_count() > 0,
                     "No Jacobian support.")

  if (bnd.size() == 0) return;

  const bool one_by_one = bnd.front().ls.one_by_one;

  if (one_by_one) {
    com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm);
  } else {
    com_data.adapt_single(bnd_qid, bnd, rovib_data, atm);
  }

  com_data.core_calc(f_grid);

  for (Index i = 0; i < f_grid.size(); ++i) {
    const auto F = Constant::sqrt_ln_2 / Constant::sqrt_pi *
                   atm[bnd_qid.Species()] * atm[bnd_qid.Isotopologue()] *
                   com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
  }
}

void equivalent_values(ExhaustiveComplexTensor3View eqv_str,
                       ExhaustiveComplexTensor3View eqv_val,
                       ComputeData& com_data,
                       const QuantumIdentifier& bnd_qid,
                       const band_data& bnd,
                       const linemixing::species_data_map& rovib_data,
                       const AtmPoint& atm,
                       const Vector& T) {
  const auto k = eqv_str.npages();
  const auto m = eqv_str.ncols();

  ARTS_USER_ERROR_IF(eqv_val.shape() != eqv_val.shape(),
                     "eqv_str and eqv_val must have the same shape.")
  ARTS_USER_ERROR_IF(T.size() != k,
                     "T must have the same size as eqv_str pages.")
  ARTS_USER_ERROR_IF(bnd.size() != static_cast<Size>(m),
                     "bnd must have the same size as eqv_str cols.")

  if (bnd.size() == 0) return;

  const bool one_by_one = bnd.front().ls.one_by_one;

  if (one_by_one) {
    com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm, false);
  } else {
    com_data.adapt_single(bnd_qid, bnd, rovib_data, atm, false);
  }

  if (not arts_omp_in_parallel() and arts_omp_get_max_threads() > 1) {
#pragma omp parallel for firstprivate(com_data)
    for (Index i = 0; i < k; ++i) {
      AtmPoint atm_copy = atm;
      atm_copy.temperature = T[i];
      if (one_by_one) {
        com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm_copy, true);
      } else {
        com_data.adapt_single(bnd_qid, bnd, rovib_data, atm_copy, true);
      }
      com_data.core_calc_eqv();
      eqv_str[i] = com_data.eqv_strs;
      eqv_val[i] = com_data.eqv_vals;
    }
  } else {
    for (Index i = 0; i < k; ++i) {
      AtmPoint atm_copy = atm;
      atm_copy.temperature = T[i];
      if (one_by_one) {
        com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm_copy, true);
      } else {
        com_data.adapt_single(bnd_qid, bnd, rovib_data, atm_copy, true);
      }
      com_data.core_calc_eqv();
      eqv_str[i] = com_data.eqv_strs;
      eqv_val[i] = com_data.eqv_vals;
    }
  }
}
}  // namespace lbl::voigt::ecs