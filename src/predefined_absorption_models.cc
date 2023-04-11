/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * @file   predefined_absorption_models.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */
#include "predefined_absorption_models.h"

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <iomanip>
#include <predef.h>

#include "debug.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "linescaling.h"
#include "matpack_data.h"
#include "predefined/predef_data.h"
#include "propagationmatrix.h"
#include "quantum_numbers.h"
#include "species.h"

namespace Absorption::PredefinedModel {
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
bool compute_selection(PropagationMatrix& pm [[maybe_unused]],
                       const SpeciesIsotopeRecord& model,
                       const Vector& f [[maybe_unused]],
                       const Numeric& p [[maybe_unused]],
                       const Numeric& t [[maybe_unused]],
                       const VMRS& vmr [[maybe_unused]],
                       const PredefinedModelData& predefined_model_data
                       [[maybe_unused]]) {
  switch (Species::find_species_index(model)) {
    case find_species_index(Species::Species::Water, "ForeignContCKDMT400"):
      if constexpr (not check_exist) MT_CKD400::compute_foreign_h2o(pm, f, p, t, vmr.H2O, predefined_model_data.get<MT_CKD400::WaterData>());
      return true;
    case find_species_index(Species::Species::Water, "SelfContCKDMT400"):
      if constexpr (not check_exist) MT_CKD400::compute_self_h2o(pm, f, p, t, vmr.H2O, predefined_model_data.get<MT_CKD400::WaterData>());
      return true;
    case find_species_index(Species::Species::Oxygen, "MPM2020"):
      if constexpr (not check_exist) MPM2020::compute(pm, f, p, t, vmr.O2);
      return true;
    case find_species_index(Species::Species::Oxygen, "PWR2021"):
      if constexpr (not check_exist) PWR20xx::compute_o2_2021(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "PWR2021"):
      if constexpr (not check_exist) PWR20xx::compute_h2o_2021(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Nitrogen, "SelfContPWR2021"):
      if constexpr (not check_exist) PWR20xx::compute_n2(pm, f, p, t, vmr.N2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Oxygen, "PWR2022"):
      if constexpr (not check_exist) PWR20xx::compute_o2_2022(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "PWR2022"):
      if constexpr (not check_exist) PWR20xx::compute_h2o_2022(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Oxygen, "PWR98"):
      if constexpr (not check_exist) PWR98::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Oxygen, "TRE05"):
      if constexpr (not check_exist) TRE05::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "PWR98"):
      if constexpr (not check_exist) PWR98::water(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Oxygen, "MPM89"):
      if constexpr (not check_exist) MPM89::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "MPM89"):
      if constexpr (not check_exist) MPM89::water(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "ForeignContCKDMT350"):
      if constexpr (not check_exist) CKDMT350::compute_foreign_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "SelfContCKDMT350"):
      if constexpr (not check_exist) CKDMT350::compute_self_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "ForeignContStandardType"):
      if constexpr (not check_exist) Standard::water_foreign(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "SelfContStandardType"):
      if constexpr (not check_exist) Standard::water_self(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Oxygen, "SelfContStandardType"):
      if constexpr (not check_exist) Standard::oxygen(pm, f, p, t, vmr.O2, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Nitrogen, "SelfContStandardType"):
      if constexpr (not check_exist) Standard::nitrogen(pm, f, p, t, vmr.N2);
      return true;
    case find_species_index(Species::Species::CarbonDioxide, "CKDMT252"):
      if constexpr (not check_exist) MT_CKD252::carbon_dioxide(pm, f, p, t, vmr.CO2);
      return true;
    case find_species_index(Species::Species::Oxygen, "visCKDMT252"):
      if constexpr (not check_exist) MT_CKD252::oxygen_vis(pm, f, p, t, vmr.O2);
      return true;
    case find_species_index(Species::Species::Nitrogen, "CIAfunCKDMT252"):
      if constexpr (not check_exist) MT_CKD252::nitrogen_fun(pm, f, p, t, vmr.N2, vmr.H2O, vmr.O2);
      return true;
    case find_species_index(Species::Species::Nitrogen, "CIArotCKDMT252"):
      if constexpr (not check_exist) MT_CKD252::nitrogen_rot(pm, f, p, t, vmr.N2, vmr.H2O, vmr.O2);
      return true;
    case find_species_index(Species::Species::Oxygen, "CIAfunCKDMT100"):
      if constexpr (not check_exist) MT_CKD100::oxygen_cia(pm, f, p, t, vmr.O2);
      return true;
    case find_species_index(Species::Species::Oxygen, "v0v0CKDMT100"):
      if constexpr (not check_exist) MT_CKD100::oxygen_v0v0(pm, f, p, t, vmr.O2, vmr.N2);
      return true;
    case find_species_index(Species::Species::Oxygen, "v1v0CKDMT100"):
      if constexpr (not check_exist) MT_CKD100::oxygen_v0v1(pm, f, p, t, vmr.O2);
      return true;
    case find_species_index(Species::Species::liquidcloud, "ELL07"):
      if constexpr (not check_exist) ELL07::compute(pm, f, t, vmr.LWC);
      return true;
    case find_species_index(Species::Species::FINAL, "Not A Model"): break;
  }
  return false;
}

bool can_compute(const SpeciesIsotopeRecord& model) {
  PropagationMatrix pm;
  return compute_selection<true>(pm, model, {}, {}, {}, {}, {});
}

/** Sets a VMR perturbation
 * 
 * @tparam special Whether or not this is called for special derivatives
 * @param[in] x Original VMR
 * @return constexpr Numeric [1e-6 if x < 1e-10 else x * 1e-6] if not special else x
 */
template <bool special>
constexpr Numeric dvmr_calc(Numeric x) noexcept {
  if constexpr (special) {
    return x;
  } else {
    constexpr Numeric d = 1e-6;
    constexpr Numeric l = d * 1e-4;
    return x < l ? d : x * d;
  }
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
template <bool special>
bool compute_vmr_deriv(PropagationMatrix& dpm,
                       const PropagationMatrix& pm,
                       const SpeciesIsotopeRecord& model,
                       const Vector& f,
                       const Numeric& p,
                       const Numeric& t,
                       VMRS vmr,
                       const Species::Species spec,
                       const PredefinedModelData& predefined_model_data
                       [[maybe_unused]]) {
  Numeric dvmr = 0;

  switch (spec) {
    case Species::Species::Oxygen:
      dvmr = dvmr_calc<special>(vmr.O2);
      if constexpr (not special) vmr.O2 += dvmr;
      break;
    case Species::Species::Water:
      dvmr = dvmr_calc<special>(vmr.H2O);
      if constexpr (not special) vmr.H2O += dvmr;
      break;
    case Species::Species::Nitrogen:
      dvmr = dvmr_calc<special>(vmr.N2);
      if constexpr (not special) vmr.N2 += dvmr;
      break;
    case Species::Species::CarbonDioxide:
      dvmr = dvmr_calc<special>(vmr.CO2);
      if constexpr (not special) vmr.CO2 += dvmr;
      break;
    case Species::Species::liquidcloud:
      dvmr = dvmr_calc<special>(vmr.LWC);
      if constexpr (not special) vmr.LWC += dvmr;
      break;
    default:
      return false;  // Escape mechanism when nothing should be done
  }
  static_assert(
      sizeof(VMRS) / sizeof(Numeric) == 5,
      "It seems you have changed VMRS.  Please check that the derivatives are up-to-date above this assert");

  if constexpr (not special) {
    dpm.SetZero();
    compute_selection<false>(dpm, model, f, p, t, vmr, predefined_model_data);
    dpm -= pm;
    dpm /= dvmr;
  } else {
    if (dvmr not_eq 0) {
      dpm = pm;
      dpm /= dvmr;
    } else {
      compute_vmr_deriv<false>(
          dpm, pm, model, f, p, t, vmr, spec, predefined_model_data);
    }
  }
  return true;
}

void compute(PropagationMatrix& propmat_clearsky,
             ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
             const SpeciesIsotopeRecord& model,
             const Vector& f_grid,
             const Numeric& rtp_pressure,
             const Numeric& rtp_temperature,
             const VMRS& vmr,
             const ArrayOfRetrievalQuantity& jacobian_quantities,
             const PredefinedModelData& predefined_model_data) {
  if (not compute_selection<true>(propmat_clearsky,
                                  model,
                                  f_grid,
                                  rtp_pressure,
                                  rtp_temperature,
                                  vmr,
                                  predefined_model_data))
    return;

  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  const bool do_vmrs_jac =
      std::any_of(jacobian_quantities.begin(),
                  jacobian_quantities.end(),
                  [](auto& deriv) { return deriv == Jacobian::Line::VMR; }) or
      std::any_of(jacobian_quantities.begin(),
                  jacobian_quantities.end(),
                  [model](auto& deriv) {
                    return deriv == Jacobian::Special::ArrayOfSpeciesTagVMR and
                           std::any_of(deriv.Target().species_array_id.begin(),
                                       deriv.Target().species_array_id.end(),
                                       [model](auto& tag) {
                                         return tag.Isotopologue() == model;
                                       });
                  });

  if (do_freq_jac or do_temp_jac or do_vmrs_jac) {
    //! Set simple propagation matrices (NOTE: stokes_dim == 1, if any model cahnges that, fix this!)
    PropagationMatrix pm(f_grid.nelem());
    PropagationMatrix dpm(f_grid.nelem());
    compute_selection<false>(pm,
                             model,
                             f_grid,
                             rtp_pressure,
                             rtp_temperature,
                             vmr,
                             predefined_model_data);

    // Add absorption to the forward parameter
    propmat_clearsky.Kjj() += pm.Kjj();

    if (do_temp_jac) {
      const Numeric d = temperature_perturbation(jacobian_quantities);
      ARTS_ASSERT(d not_eq 0)

      dpm.SetZero();
      compute_selection<false>(dpm,
                               model,
                               f_grid,
                               rtp_pressure,
                               rtp_temperature + d,
                               vmr,
                               predefined_model_data);
      dpm -= pm;
      dpm /= d;
      for (Index iq = 0; iq < dpropmat_clearsky_dx.nelem(); iq++) {
        if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
        }
      }
    }

    if (do_freq_jac) {
      const Numeric d = frequency_perturbation(jacobian_quantities);
      ARTS_ASSERT(d not_eq 0)

      Vector f_grid_d{f_grid};
      f_grid_d += d;

      dpm.SetZero();
      compute_selection<false>(dpm,
                               model,
                               f_grid_d,
                               rtp_pressure,
                               rtp_temperature,
                               vmr,
                               predefined_model_data);
      dpm -= pm;
      dpm /= d;
      for (Index iq = 0; iq < dpropmat_clearsky_dx.nelem(); iq++) {
        if (is_frequency_parameter(jacobian_quantities[iq])) {
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
        }
      }
    }

    for (Index iq = 0; iq < dpropmat_clearsky_dx.nelem(); iq++) {
      auto& deriv = jacobian_quantities[iq];
      if (deriv == Jacobian::Line::VMR) {
        if (compute_vmr_deriv<false>(dpm,
                                     pm,
                                     model,
                                     f_grid,
                                     rtp_pressure,
                                     rtp_temperature,
                                     vmr,
                                     deriv.QuantumIdentity().Species(),
                                     predefined_model_data))
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
      } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR and
                 std::any_of(deriv.Target().species_array_id.begin(),
                             deriv.Target().species_array_id.end(),
                             [model](auto& tag) {
                               return tag.Isotopologue() == model;
                             })) {
        if (compute_vmr_deriv<true>(
                dpm,
                pm,
                model,
                f_grid,
                rtp_pressure,
                rtp_temperature,
                vmr,
                deriv.Target().species_array_id.front().Spec(),
                predefined_model_data))
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
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
}  // namespace Absorption::PredefinedModel