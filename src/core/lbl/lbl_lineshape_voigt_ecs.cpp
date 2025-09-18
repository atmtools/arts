#include "lbl_lineshape_voigt_ecs.h"

#include <arts_omp.h>
#include <atm.h>
#include <configtypes.h>
#include <debug.h>
#include <isotopologues.h>
#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <exception>
#include <ranges>
#include <stdexcept>

#include "lbl_lineshape_linemixing.h"
#include "lbl_lineshape_model.h"
#include "lbl_lineshape_voigt_ecs_hartmann.h"
#include "lbl_lineshape_voigt_ecs_makarov.h"

#undef WIGNER3
#undef WIGNER6

namespace lbl::voigt::ecs {
ComputeData::ComputeData(const ConstVectorView& f_grid,
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

  for (Size k = 0; k < m; k++) {
    auto V = Vs[k];
    auto W = Ws[k];
    inplace_transpose(W);
    auto eqv_str = eqv_strs[k];
    auto eqv_val = eqv_vals[k];

    diagonalize(V, eqv_val, W);

    // Do the matrix forward multiplication
    for (Size i = 0; i < n; i++) {
      for (Size j = 0; j < n; j++) {
        eqv_str[i] += dip[j] * V[j, i];
      }
    }

    // Do the matrix backward multiplication
    inv(V, V);
    for (Size i = 0; i < n; i++) {
      Complex z(0, 0);
      for (Size j = 0; j < n; j++) {
        z += pop[j] * dip[j] * V[i, j];
      }
      eqv_str[i] *= z;
    }
  }
}

void ComputeData::core_calc(const ConstVectorView& f_grid) try {
  core_calc_eqv();

  const auto m = vmrs.size();
  const auto n = f_grid.size();
  shape        = 0;

  for (Size k = 0; k < m; k++) {
    for (Size i = 0; i < eqv_strs[k].size(); i++) {
      const Numeric gamd = gd_fac * eqv_vals[k][i].real();
      const Numeric cte  = Constant::sqrt_ln_2 / gamd;
      for (Size iv = 0; iv < n; iv++) {
        const Complex z  = (eqv_vals[k][i] - f_grid[iv]) * cte;
        const Complex w  = Faddeeva::w(z);
        shape[iv]       += vmrs[k] * eqv_strs[k][i] * w / gamd;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
void get_vmrs(VectorView vmrs,
              const line_shape::model::map_t& mod,
              const AtmPoint& atm) {
  std::transform(mod.begin(), mod.end(), vmrs.begin(), [&atm](auto& m) {
    return m.first == SpeciesEnum::Bath ? 0.0 : atm[m.first];
  });

  const Size bath_spec =
      std::distance(mod.begin(), mod.find(SpeciesEnum::Bath));

  if (bath_spec != mod.size()) {
    vmrs[bath_spec] = 1 - sum(vmrs);
  } else {
    vmrs /= sum(vmrs);
  }
}
}  // namespace

void ComputeData::adapt_multi(const QuantumIdentifier& bnd_qid,
                              const band_data& bnd,
                              const LinemixingSpeciesEcsData& rovib_data,
                              const AtmPoint& atm,
                              const bool presorted) try {
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

  Ws       = 0;
  vmrs     = 0;
  eqv_strs = 0;
  eqv_vals = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared *
                     atm.temperature / bnd_qid.isot.mass);

  const Numeric T  = atm.temperature;
  const Numeric QT = PartitionFunctions::Q(T, bnd_qid.isot);

  for (Size i = 0; i < n; i++) {
    const auto& line = bnd.lines[i];
    pop[i]           = line.gu * exp(-line.e0 / (Constant::k * T)) / QT;
    dip[i]           = 0.5 * Constant::c *
             std::sqrt(line.a / (Math::pow3(line.f0) * Constant::two_pi));

    if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) {
      auto& J  = bnd.lines[i].qn.at(QuantumNumberType::J);
      auto& N  = bnd.lines[i].qn.at(QuantumNumberType::N);
      dipr[i]  = makarov::reduced_dipole(J.upper.get<Rational>(),
                                        J.lower.get<Rational>(),
                                        N.upper.get<Rational>());
      dip[i]  *= std::signbit(dipr[i]) ? -1 : 1;
    } else if (bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN) {
      auto& J   = bnd.lines[i].qn.at(QuantumNumberType::J);
      auto& l2  = bnd_qid.state.at(QuantumNumberType::l2);
      dipr[i]   = hartmann::reduced_dipole(J.upper.get<Rational>(),
                                         J.lower.get<Rational>(),
                                         l2.upper.get<Rational>(),
                                         l2.lower.get<Rational>());
      dip[i]   *= std::signbit(dipr[i]) ? -1 : 1;
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

    pop  = presorter(pop);
    dip  = presorter(dip);
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
                                return lsl.first == lsr.first;
                              });

    ARTS_USER_ERROR_IF(not lines_have_same_species, "Bad species combination")
  }

  get_vmrs(vmrs, bnd.front().ls.single_models, atm);

  Size i = 0;
  for (auto& spec : bnd.front().ls.single_models | stdv::keys) {
    const auto& rovib_data_it = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(
        rovib_data_it == rovib_data.end(), "No rovib data for species {}", spec)

    Wimag = 0.0;
    for (Size k = 0; k < n; k++) {
      Wimag[k, k] = bnd.lines[sort[k]].ls.single_models.at(spec).G0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
      Ws[i][k, k] = bnd.lines[sort[k]].ls.single_models.at(spec).D0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
    }

    if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) {
      makarov::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else if (bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN) {
      hartmann::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else {
      ARTS_USER_ERROR("UNKNOWN ECS LINE SHAPE {}", bnd.lineshape)
    }

    Ws[i].imag() = Wimag;
    i++;
  }

  for (Size i = 0; i < n; i++) {
    Ws[joker, i, i] += bnd.lines[sort[i]].f0;
  }
}
ARTS_METHOD_ERROR_CATCH

void ComputeData::adapt_single(const QuantumIdentifier& bnd_qid,
                               const band_data& bnd,
                               const LinemixingSpeciesEcsData& rovib_data,
                               const AtmPoint& atm,
                               const bool presorted) try {
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

  Ws       = 0;
  vmrs     = 1;
  eqv_strs = 0;
  eqv_vals = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared *
                     atm.temperature / bnd_qid.isot.mass);

  const Numeric T  = atm.temperature;
  const Numeric QT = PartitionFunctions::Q(T, bnd_qid.isot);

  for (Size i = 0; i < n; i++) {
    const auto& line = bnd.lines[i];
    pop[i]           = line.gu * exp(-line.e0 / (Constant::k * T)) / QT;
    dip[i]           = 0.5 * Constant::c *
             std::sqrt(line.a / (Math::pow3(line.f0) * Constant::two_pi));

    if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) {
      auto& J  = bnd.lines[i].qn.at(QuantumNumberType::J);
      auto& N  = bnd.lines[i].qn.at(QuantumNumberType::N);
      dipr[i]  = makarov::reduced_dipole(J.upper.get<Rational>(),
                                        J.lower.get<Rational>(),
                                        N.upper.get<Rational>());
      dip[i]  *= std::signbit(dipr[i]) ? -1 : 1;
    } else if (bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN) {
      auto& J   = bnd.lines[i].qn.at(QuantumNumberType::J);
      auto& l2  = bnd_qid.state.at(QuantumNumberType::l2);
      dipr[i]   = hartmann::reduced_dipole(J.upper.get<Rational>(),
                                         J.lower.get<Rational>(),
                                         l2.upper.get<Rational>(),
                                         l2.lower.get<Rational>());
      dip[i]   *= std::signbit(dipr[i]) ? -1 : 1;
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

    pop  = presorter(pop);
    dip  = presorter(dip);
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
                                return lsl.first == lsr.first;
                              });

    ARTS_USER_ERROR_IF(not lines_have_same_species, "Bad species combination")
  }

  get_vmrs(vmrs, bnd.front().ls.single_models, atm);

  Size i = 0;
  for (auto& spec : bnd.front().ls.single_models | stdv::keys) {
    const auto& rovib_data_it = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(
        rovib_data_it == rovib_data.end(), "No rovib data for species {}", spec)

    const Numeric this_vmr = vmrs[i];

    Wimag = 0.0;
    for (Size k = 0; k < n; k++) {
      Wimag[k, k] = bnd.lines[sort[k]].ls.single_models.at(spec).G0(
          bnd.lines[sort[k]].ls.T0, atm.temperature, atm.pressure);
    }

    if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) {
      makarov::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else if (bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN) {
      hartmann::relaxation_matrix_offdiagonal(
          Wimag, bnd_qid, bnd, sort, spec, rovib_data_it->second, dipr, atm);
    } else {
      ARTS_USER_ERROR("UNKNOWN ECS LINE SHAPE {}", bnd.lineshape)
    }

    for (Size ir = 0; ir < n; ir++) {
      for (Size ic = 0; ic < n; ic++) {
        imag_val(Ws[0][ir, ic]) += this_vmr * Wimag[ir, ic];
      }
    }

    i++;
  }

  for (Size i = 0; i < n; i++) {
    real_val(Ws[0][i, i]) =
        bnd.lines[sort[i]].f0 + bnd.lines[sort[i]].ls.D0(atm);
  }
}
ARTS_METHOD_ERROR_CATCH

void calculate(PropmatVectorView pm_,
               PropmatMatrixView,
               ComputeData& com_data,
               const ConstVectorView f_grid_,
               const Range& f_range,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const LinemixingSpeciesEcsData& rovib_data,
               const AtmPoint& atm,
               const zeeman::pol pol,
               const bool no_negative_absorption) try {
  if (pol != zeeman::pol::no) {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            bnd, [](auto& zee) { return zee.on; }, &line::z),
        "Zeeman effect and ECS in combination is not yet possible.")
    return;
  }

  PropmatVectorView pm         = pm_[f_range];
  const ConstVectorView f_grid = f_grid_[f_range];

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

  for (Size i = 0; i < f_grid.size(); ++i) {
    const auto F = Constant::sqrt_ln_2 / Constant::sqrt_pi *
                   atm[bnd_qid.isot.spec] * atm[bnd_qid.isot] *
                   com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
  }
}
ARTS_METHOD_ERROR_CATCH

void equivalent_values(ComplexTensor3View eqv_str,
                       ComplexTensor3View eqv_val,
                       ComputeData& com_data,
                       const QuantumIdentifier& bnd_qid,
                       const band_data& bnd,
                       const LinemixingSpeciesEcsData& rovib_data,
                       const AtmPoint& atm,
                       const Vector& T) try {
  const auto k = eqv_str.npages();
  const auto m = eqv_str.ncols();

  ARTS_USER_ERROR_IF(eqv_val.shape() != eqv_val.shape(),
                     "eqv_str and eqv_val must have the same shape.")
  ARTS_USER_ERROR_IF(T.size() != static_cast<Size>(k),
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

  std::string err{};
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(com_data)
  for (Index i = 0; i < k; ++i) {
    try {
      AtmPoint atm_copy    = atm;
      atm_copy.temperature = T[i];
      if (one_by_one) {
        com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm_copy, true);
      } else {
        com_data.adapt_single(bnd_qid, bnd, rovib_data, atm_copy, true);
      }
      com_data.core_calc_eqv();
      eqv_str[i] = com_data.eqv_strs;
      eqv_val[i] = com_data.eqv_vals;
    } catch (std::exception& e) {
#pragma omp critical
      err += std::format("{}\n", e.what());
    }
  }

  if (not err.empty()) throw std::runtime_error(err);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs