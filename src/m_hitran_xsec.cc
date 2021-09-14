/* Copyright (C) 2017 Oliver Lemke <oliver.lemke@uni-hamburg.de> and
                      Stefan Buehler <stefan.buehler@uni-hamburg.de>.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   m_hitran_xsec.cc
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2021-02-23

  \brief  Workspace methods for HITRAN absorption cross section data.

*/

#include "arts.h"
#include "global_data.h"
#include "hitran_xsec.h"
#include "jacobian.h"
#include "m_xml.h"
#include "messages.h"
#include "physics_funcs.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadXsecData(ArrayOfXsecRecord& hitran_xsec_data,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const String& basename,
                  const Verbosity& verbosity) {
  // Build a set of species indices. Duplicates are ignored.
  std::set<Species::Species> unique_species;
  for (auto& asp : abs_species) {
    for (auto& sp : asp) {
      if (sp.Type() == Species::TagType::HitranXsec) {
        unique_species.insert(sp.Spec());
      }
    }
  }

  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }

  // Read xsec data for all active species and collect them in hitran_xsec_data
  hitran_xsec_data.clear();
  for (auto& species_name : unique_species) {
    XsecRecord xsec_coeffs;
    const String filename{tmpbasename +
                          String(Species::toShortName(species_name)) + ".xml"};

    try {
      ReadXML(xsec_coeffs, "", filename, "", verbosity);

      hitran_xsec_data.push_back(xsec_coeffs);
    } catch (const std::exception& e) {
      ARTS_USER_ERROR(
          "Error reading coefficients file:\n", filename, "\n", e.what());
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddHitranXsec(  // WS Output:
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const Vector& rtp_vmr,
    const ArrayOfXsecRecord& hitran_xsec_data,
    const Numeric& force_p,
    const Numeric& force_t,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUTS;

  // Forward simulations and their error handling
  ARTS_USER_ERROR_IF(rtp_vmr.nelem() not_eq abs_species.nelem(),
                     "Mismatch dimensions on species and VMR inputs");
  ARTS_USER_ERROR_IF(
      propmat_clearsky.NumberOfFrequencies() not_eq f_grid.nelem(),
      "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  if (dpropmat_clearsky_dx.nelem()) {
    ARTS_USER_ERROR_IF(
        dpropmat_clearsky_dx.nelem() not_eq jacobian_quantities.nelem(),
        "Mismatch dimensions on xsec derivatives and Jacobian grids");
    ARTS_USER_ERROR_IF(
        std::any_of(dpropmat_clearsky_dx.cbegin(),
                    dpropmat_clearsky_dx.cend(),
                    [&f_grid](auto& x) {
                      return x.NumberOfFrequencies() not_eq f_grid.nelem();
                    }),
        "Mismatch dimensions on internal matrices of xsec derivatives and frequency");
  }

  // Jacobian overhead START
  /* NOTE:  The calculations below are inefficient and could
              be made much better by using interp in Extract to
              return the derivatives as well. */
  // Jacobian vectors START
  Vector dxsec_temp_dT;
  Vector dxsec_temp_dF;
  Vector dfreq;
  // Jacobian vectors END
  const bool do_jac = supports_hitran_xsec(jacobian_quantities);
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);
  if (do_freq_jac) {
    dfreq.resize(f_grid.nelem());
    dfreq = f_grid;
    dfreq += df;
    dxsec_temp_dF.resize(f_grid.nelem());
  }
  if (do_temp_jac) {
    dxsec_temp_dT.resize(f_grid.nelem());
  }
  // Jacobian overhead END

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  // Allocate a vector with dimension frequencies for constructing our
  // cross-sections before adding them (more efficient to allocate this here
  // outside of the loops)
  Vector xsec_temp(f_grid.nelem(), 0.);

  ArrayOfString fail_msg;
  bool do_abort = false;
  // Loop over Xsec data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (Index i = 0; i < abs_species.nelem(); i++) {
    for (Index s = 0; s < abs_species[i].nelem(); s++) {
      const SpeciesTag& this_species = abs_species[i][s];

      // Check if this is a HITRAN cross section tag
      if (this_species.Type() != Species::TagType::HitranXsec) continue;

      Index this_xdata_index =
          hitran_xsec_get_index(hitran_xsec_data, this_species.Spec());
      ARTS_USER_ERROR_IF(this_xdata_index < 0,
                         "Cross-section species ",
                         this_species.Name(),
                         " not found in *hitran_xsec_data*.")
      const XsecRecord& this_xdata = hitran_xsec_data[this_xdata_index];
      // ArrayOfMatrix& this_dxsec = do_jac ? dpropmat_clearsky_dx[i] : empty;

      if (do_abort) continue;
      const Numeric current_p = force_p < 0 ? rtp_pressure : force_p;
      const Numeric current_t = force_t < 0 ? rtp_temperature : force_t;

      // Get the absorption cross sections from the HITRAN data:
      this_xdata.Extract(xsec_temp,
                         f_grid,
                         current_p,
                         current_t,
                         verbosity);
      if (do_freq_jac)
        this_xdata.Extract(dxsec_temp_dF,
                           dfreq,
                           current_p,
                           current_t,
                           verbosity);
      if (do_temp_jac)
        this_xdata.Extract(dxsec_temp_dT,
                           f_grid,
                           current_p,
                           current_t + dt,
                           verbosity);
    }

    // Add to result variable:
    Numeric nd = number_density(rtp_pressure, rtp_temperature);
    if (!do_jac) {
      xsec_temp *= nd * rtp_vmr[i];
      propmat_clearsky.Kjj() += xsec_temp;
    } else {
      for (Index f = 0; f < f_grid.nelem(); f++) {
        propmat_clearsky.Kjj()[f] += xsec_temp[f] * nd * rtp_vmr[i];
        for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
          const auto& deriv = jacobian_quantities[iq];

          if (is_frequency_parameter(deriv)) {
            dpropmat_clearsky_dx[iq].Kjj()[f] +=
                (dxsec_temp_dF[f] - xsec_temp[f]) * nd * rtp_vmr[i] / df;
          }

          if (!deriv.propmattype()) continue;

          if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR ||
              deriv == Jacobian::Line::VMR) {
            if (species_match(deriv, abs_species[i])) {
              dpropmat_clearsky_dx[iq].Kjj()[f] +=
                  xsec_temp[f] * nd * rtp_vmr[i];
            }
          } else if (deriv == Jacobian::Atm::Temperature) {
            dpropmat_clearsky_dx[iq].Kjj()[f] +=
                (dxsec_temp_dT[f] - xsec_temp[f]) * nd * rtp_vmr[i] / dt;
          }
        }
      }
    }
  }
}
