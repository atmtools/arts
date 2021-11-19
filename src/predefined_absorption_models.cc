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

#include "debug.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "linescaling.h"
#include "matpack.h"
#include "matpackI.h"
#include "predefined/predef.h"
#include "propagationmatrix.h"
#include "quantum.h"
#include "species.h"

namespace Absorption::PredefinedModel {
/** Compute the selected model
*
* Remember to use the "if constexpr (not check_exist)" statement
* as this allows us to not keep any copy of the list of available
* models
*
*/
template <bool check_exist = false>
bool compute_selection(PropagationMatrix& pm [[maybe_unused]],
                       const SpeciesIsotopeRecord& model,
                       const Vector& f [[maybe_unused]],
                       const Numeric& p [[maybe_unused]],
                       const Numeric& t [[maybe_unused]],
                       const VMRS& vmr [[maybe_unused]]) {
  switch (Species::find_species_index(model)) {
    case find_species_index(Species::Species::Oxygen, "MPM2020"):
      if constexpr (not check_exist) MPM2020::compute(pm, f, p, t, vmr.O2);
      return true;
    case find_species_index(Species::Species::Water, "ForeignContCKDMT350"):
      if constexpr (not check_exist)
        CKDMT350::compute_foreign_h2o(pm, f, p, t, vmr.H2O);
      return true;
    case find_species_index(Species::Species::Water, "SelfContCKDMT350"):
      if constexpr (not check_exist)
        CKDMT350::compute_self_h2o(pm, f, p, t, vmr.H2O);
      return true;
  }
  return false;
}

void compute(PropagationMatrix& propmat_clearsky,
             ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
             const SpeciesIsotopeRecord& model,
             const Vector& f_grid,
             const Numeric& rtp_pressure,
             const Numeric& rtp_temperature,
             const VMRS& vmr,
             const ArrayOfRetrievalQuantity& jacobian_quantities) {
  if (not compute_selection<true>(
          propmat_clearsky, model, f_grid, rtp_pressure, rtp_temperature, vmr))
    return;

  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);

  if (do_freq_jac or do_temp_jac) {
    //! Set simple propagation matrices
    PropagationMatrix pm(
        f_grid.nelem());  // NOTE: Change if ever stokes_dim not_eq 1
    PropagationMatrix dpm(
        f_grid.nelem());  // NOTE: Change if ever stokes_dim not_eq 1
    compute_selection(pm, model, f_grid, rtp_pressure, rtp_temperature, vmr);

    // Add absorption to the forward parameter
    propmat_clearsky.Kjj() += pm.Kjj();

    if (do_temp_jac) {
      const Numeric d = temperature_perturbation(jacobian_quantities);
      ARTS_ASSERT(d not_eq 0)

      compute_selection(
          dpm, model, f_grid, rtp_pressure, rtp_temperature + d, vmr);
      dpm -= pm;
      dpm /= d;
      for (Index iq = 0; iq < dpropmat_clearsky_dx.nelem(); iq++) {
        if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
        }
      }

      dpm.SetZero();
    }

    if (do_freq_jac) {
      const Numeric d = temperature_perturbation(jacobian_quantities);
      ARTS_ASSERT(d not_eq 0)

      Vector f_grid_d{f_grid};
      f_grid_d += d;

      compute_selection(
          dpm, model, f_grid_d, rtp_pressure, rtp_temperature, vmr);
      dpm -= pm;
      dpm /= d;
      for (Index iq = 0; iq < dpropmat_clearsky_dx.nelem(); iq++) {
        if (is_frequency_parameter(jacobian_quantities[iq])) {
          dpropmat_clearsky_dx[iq].Kjj() += dpm.Kjj();
        }
      }

      // ADD IF MORE DERIVATIVES: dpm.SetZero();
    }

    // Done!
  } else {
    compute_selection(
        propmat_clearsky, model, f_grid, rtp_pressure, rtp_temperature, vmr);
  }
}
}  // namespace Absorption::PredefinedModel