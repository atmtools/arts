/*!
 * @file   predefined_absorption_models.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */
#include "predefined_absorption_models.h"

#include <predef.h>

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <stdexcept>
#include <variant>

#include "atm.h"
#include "debug.h"
#include "isotopologues.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "predefined/predef_data.h"
#include "species.h"

namespace Absorption::PredefinedModel {
VMRS::VMRS(const AtmPoint& atm)
    : CO2(atm[species[0]]),
      O2(atm[species[1]]),
      N2(atm[species[2]]),
      H2O(atm[species[3]]),
      LWC(atm[species[4]]) {}

/** Compute the selected model and returns if it can be computed
 *
 * Remember to use the "if constexpr (not check_exist)" statement
 * as this allows us to not keep any copy of the list of available
 * models
 * 
 * @tparam check_exist Perform no computations if false
 * @param[inout] pm A local propagation matrix
 * @param[in] model A single isotope record
 * @param[in] f A local frequency grid
 * @param[in] p A local pressure
 * @param[in] t A local temperature
 * @param[in] vmr A VMRS object defined from the WSVs abs_species and rtp_vmr
 * @return true When there are computations that can be or have been performed
 * @return false When there are no computations that can be or have been performed
 */
template <bool check_exist>
bool compute_selection(
    PropmatVector& pm [[maybe_unused]],
    const SpeciesIsotope& model,
    const Vector& f [[maybe_unused]],
    const Numeric& p [[maybe_unused]],
    const Numeric& t [[maybe_unused]],
    const VMRS& vmr [[maybe_unused]],
    const Absorption::PredefinedModel::ModelVariant& predefined_model_data
    [[maybe_unused]]) try {
  switch (Species::find_species_index(model)) {
    case "H2O-ForeignContCKDMT400"_isot_index:
      if constexpr (not check_exist)
        MT_CKD400::compute_foreign_h2o(
            pm,
            f,
            p,
            t,
            vmr.H2O,
            std::get<MT_CKD400::WaterData>(predefined_model_data));
      return true;
    case "H2O-SelfContCKDMT400"_isot_index:
      if constexpr (not check_exist)
        MT_CKD400::compute_self_h2o(
            pm,
            f,
            p,
            t,
            vmr.H2O,
            std::get<MT_CKD400::WaterData>(predefined_model_data));
      return true;
    case "O2-MPM2020"_isot_index:
      if constexpr (not check_exist) MPM2020::compute(pm, f, p, t, vmr.O2);
      return true;
    case "O2-PWR2021"_isot_index:
      if constexpr (not check_exist)
        PWR20xx::compute_o2_2021(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "H2O-PWR2021"_isot_index:
      if constexpr (not check_exist)
        PWR20xx::compute_h2o_2021(pm, f, p, t, vmr.H2O);
      return true;
    case "N2-SelfContPWR2021"_isot_index:
      if constexpr (not check_exist)
        PWR20xx::compute_n2(pm, f, p, t, vmr.N2, vmr.H2O);
      return true;
    case "O2-PWR2022"_isot_index:
      if constexpr (not check_exist)
        PWR20xx::compute_o2_2022(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "H2O-PWR2022"_isot_index:
      if constexpr (not check_exist)
        PWR20xx::compute_h2o_2022(pm, f, p, t, vmr.H2O);
      return true;
    case "O2-PWR98"_isot_index:
      if constexpr (not check_exist)
        PWR98::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "O2-TRE05"_isot_index:
      if constexpr (not check_exist)
        TRE05::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "H2O-PWR98"_isot_index:
      if constexpr (not check_exist) PWR98::water(pm, f, p, t, vmr.H2O);
      return true;
    case "O2-MPM89"_isot_index:
      if constexpr (not check_exist)
        MPM89::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "H2O-MPM89"_isot_index:
      if constexpr (not check_exist) MPM89::water(pm, f, p, t, vmr.H2O);
      return true;
    case "N2-SelfContMPM93"_isot_index:
      if constexpr (not check_exist)
        MPM93::nitrogen(pm, f, p, t, vmr.N2, vmr.H2O);
      return true;
    case "H2O-ForeignContCKDMT350"_isot_index:
      if constexpr (not check_exist)
        CKDMT350::compute_foreign_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case "H2O-SelfContCKDMT350"_isot_index:
      if constexpr (not check_exist)
        CKDMT350::compute_self_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case "H2O-ForeignContCKDMT320"_isot_index:
      if constexpr (not check_exist)
        CKDMT320::compute_foreign_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case "H2O-SelfContCKDMT320"_isot_index:
      if constexpr (not check_exist)
        CKDMT320::compute_self_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case "H2O-ForeignContStandardType"_isot_index:
      if constexpr (not check_exist)
        Standard::water_foreign(pm, f, p, t, vmr.H2O);
      return true;
    case "H2O-SelfContStandardType"_isot_index:
      if constexpr (not check_exist) Standard::water_self(pm, f, p, t, vmr.H2O);
      return true;
    case "O2-SelfContStandardType"_isot_index:
      if constexpr (not check_exist)
        Standard::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case "N2-SelfContStandardType"_isot_index:
      if constexpr (not check_exist) Standard::nitrogen(pm, f, p, t, vmr.N2);
      return true;
    case "CO2-CKDMT252"_isot_index:
      if constexpr (not check_exist)
        MT_CKD252::carbon_dioxide(pm, f, p, t, vmr.CO2);
      return true;
    case "O2-visCKDMT252"_isot_index:
      if constexpr (not check_exist) MT_CKD252::oxygen_vis(pm, f, p, t, vmr.O2);
      return true;
    case "N2-CIAfunCKDMT252"_isot_index:
      if constexpr (not check_exist)
        MT_CKD252::nitrogen_fun(pm, f, p, t, vmr.N2, vmr.H2O, vmr.O2);
      return true;
    case "N2-CIArotCKDMT252"_isot_index:
      if constexpr (not check_exist)
        MT_CKD252::nitrogen_rot(pm, f, p, t, vmr.N2, vmr.H2O, vmr.O2);
      return true;
    case "O2-CIAfunCKDMT100"_isot_index:
      if constexpr (not check_exist) MT_CKD100::oxygen_cia(pm, f, p, t, vmr.O2);
      return true;
    case "O2-v0v0CKDMT100"_isot_index:
      if constexpr (not check_exist)
        MT_CKD100::oxygen_v0v0(pm, f, p, t, vmr.O2, vmr.N2);
      return true;
    case "O2-v1v0CKDMT100"_isot_index:
      if constexpr (not check_exist)
        MT_CKD100::oxygen_v0v1(pm, f, p, t, vmr.O2);
      return true;
    case "liquidcloud-ELL07"_isot_index:
      if constexpr (not check_exist) ELL07::compute(pm, f, t, vmr.LWC);
      return true;
  }
  return false;
} catch (std::bad_variant_access& e) {
  throw std::runtime_error(
      var_string("Data for model ", model, " is not of correct type"));
} catch (std::exception&) {
  throw;
}

bool can_compute(const SpeciesIsotope& model) {
  PropmatVector pm;
  return compute_selection<true>(pm, model, {}, {}, {}, {}, {});
}

/** Compute the partial VMR derivative
 *
 * Sets dpm to the expected value (does not add, but really sets)
 *
 * Note that there is an extra special case when the template argument is
 * true to handle the 0-vmr case
 * 
 * @tparam special Whether or not this is called for special derivatives
 * @param[inout] dpm The 
 * @param[in] pm A local propagation matrix
 * @param[in] model A single isotope record
 * @param[in] f A local frequency grid
 * @param[in] p A local pressure
 * @param[in] t A local temperature
 * @param[in] vmr A VMRS object defined from the WSVs abs_species and rtp_vmr
 * @param[in] spec The species whose derivative is computed
 */
void compute_vmr_deriv(
    PropmatVector& dpm,
    const PropmatVector& pm,
    const SpeciesIsotope& model,
    const Vector& f,
    const Numeric& p,
    const Numeric& t,
    VMRS vmr,
    const Numeric dvmr,
    const SpeciesEnum spec,
    const Absorption::PredefinedModel::ModelVariant& predefined_model_data
    [[maybe_unused]]) {
  switch (spec) {
    case SpeciesEnum::Oxygen:
      vmr.O2 += dvmr;
      break;
    case SpeciesEnum::Water:
      vmr.H2O += dvmr;
      break;
    case SpeciesEnum::Nitrogen:
      vmr.N2 += dvmr;
      break;
    case SpeciesEnum::CarbonDioxide:
      vmr.CO2 += dvmr;
      break;
    case SpeciesEnum::liquidcloud:
      vmr.LWC += dvmr;
      break;
    default:
      ARTS_ASSERT(false);
  }
  static_assert(
      sizeof(VMRS) / sizeof(Numeric) == 5,
      "It seems you have changed VMRS.  Please check that the derivatives are up-to-date above this assert");

  dpm = 0;
  compute_selection<false>(dpm, model, f, p, t, vmr, predefined_model_data);
  dpm -= pm;
  dpm /= dvmr;
}

void compute(
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotope& model,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const VMRS& vmr,
    const JacobianTargets& jacobian_targets,
    const Absorption::PredefinedModel::ModelVariant& predefined_model_data) {
  if (not compute_selection<true>(propmat_clearsky,
                                  model,
                                  f_grid,
                                  rtp_pressure,
                                  rtp_temperature,
                                  vmr,
                                  predefined_model_data))
    return;

  using enum SpeciesEnum;
  const auto freq_jac = jacobian_targets.find_all<Jacobian::AtmTarget>(
      AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w);
  const auto temp_jac = jacobian_targets.find<Jacobian::AtmTarget>(AtmKey::t);
  const auto vmrs_jac = jacobian_targets.find_all<Jacobian::AtmTarget>(
      CarbonDioxide, Oxygen, Nitrogen, Water, liquidcloud);
  const bool do_freq_jac =
      std::ranges::any_of(freq_jac, [](auto& x) { return x.first; });
  const bool do_temp_jac = temp_jac.first;
  const bool do_vmrs_jac =
      std::ranges::any_of(vmrs_jac, [](auto& x) { return x.first; });

  if (do_freq_jac or do_temp_jac or do_vmrs_jac) {
    PropmatVector pm(f_grid.nelem());
    PropmatVector dpm(f_grid.nelem());
    compute_selection<false>(pm,
                             model,
                             f_grid,
                             rtp_pressure,
                             rtp_temperature,
                             vmr,
                             predefined_model_data);

    // Add absorption to the forward parameter
    propmat_clearsky += pm;

    if (do_temp_jac) {
      const Numeric d = temp_jac.second->d;
      const auto iq   = temp_jac.second->target_pos;
      ARTS_ASSERT(d not_eq 0)

      dpm = 0;
      compute_selection<false>(dpm,
                               model,
                               f_grid,
                               rtp_pressure,
                               rtp_temperature + d,
                               vmr,
                               predefined_model_data);
      dpm                      -= pm;
      dpm                      /= d;
      dpropmat_clearsky_dx[iq] += dpm;
    }

    for (auto& j : freq_jac) {
      if (j.first) {
        const Numeric d = j.second->d;
        const auto iq   = j.second->target_pos;
        ARTS_ASSERT(d not_eq 0)

        Vector f_grid_d{f_grid};
        f_grid_d += d;

        dpm = 0.0;
        compute_selection<false>(dpm,
                                 model,
                                 f_grid_d,
                                 rtp_pressure,
                                 rtp_temperature,
                                 vmr,
                                 predefined_model_data);
        dpm                      -= pm;
        dpm                      /= d;
        dpropmat_clearsky_dx[iq] += dpm;
      }
    }

    for (auto& j : vmrs_jac) {
      if (j.first) {
        const Numeric d = j.second->d;
        const auto iq   = j.second->target_pos;
        compute_vmr_deriv(dpm,
                          pm,
                          model,
                          f_grid,
                          rtp_pressure,
                          rtp_temperature,
                          vmr,
                          d,
                          *std::get_if<SpeciesEnum>(&j.second->type),
                          predefined_model_data);
        dpropmat_clearsky_dx[iq] += dpm;
      }
    }
  } else {
    compute_selection<false>(propmat_clearsky,
                             model,
                             f_grid,
                             rtp_pressure,
                             rtp_temperature,
                             vmr,
                             predefined_model_data);
  }
}

std::ostream& operator<<(std::ostream& os, const VMRS& vmrs) {
  return os << "O2: " << vmrs.O2 << '\n'
            << "N2: " << vmrs.N2 << '\n'
            << "H2O: " << vmrs.H2O << '\n'
            << "LWC: " << vmrs.LWC << '\n'
            << "CO2: " << vmrs.CO2 << '\n';
}
}  // namespace Absorption::PredefinedModel