/*!
  \file   m_hitran_xsec.cc
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2021-02-23

  \brief  Workspace methods for HITRAN absorption cross section data.

*/

#include <workspace.h>

#include "jacobian.h"
#include "physics_funcs.h"
#include "species.h"
#include "xml_io.h"
#include "xsec_fit.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_xsec_fit_dataReadSpeciesSplitCatalog(
    ArrayOfXsecRecord& absorption_xsec_fit_data,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& basename,
    const Index& ignore_missing_) try {
  const bool ignore_missing = static_cast<bool>(ignore_missing_);

  // Build a set of species indices. Duplicates are ignored.
  std::set<SpeciesEnum> unique_species;
  for (auto& asp : abs_species) {
    for (auto& sp : asp) {
      if (sp.Type() == SpeciesTagType::XsecFit) {
        unique_species.insert(sp.Spec());
      }
    }
  }

  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }

  // Read xsec data for all active species and collect them in absorption_xsec_fit_data
  absorption_xsec_fit_data.clear();
  for (auto& species_name : unique_species) {
    XsecRecord xsec_coeffs;
    String filename{std::format("{}{}-XFIT.xml", tmpbasename, species_name)};

    if (not find_xml_file_existence(filename)) {
      if (ignore_missing) continue;
      ARTS_USER_ERROR("File {} not found", filename);
    }

    xml_read_from_file_base(filename, xsec_coeffs);

    absorption_xsec_fit_data.push_back(std::move(xsec_coeffs));
  }
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixAddXsecFit(  // WS Output:
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    // WS Input:
    const SpeciesEnum& select_species,
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& f_grid,
    const AtmPoint& atm_point,
    const ArrayOfXsecRecord& absorption_xsec_fit_data,
    const Numeric& force_p,
    const Numeric& force_t) {
  // Forward simulations and their error handling
  ARTS_USER_ERROR_IF(
      propagation_matrix.size() not_eq f_grid.size(),
      "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  ARTS_USER_ERROR_IF(
      static_cast<Size>(propagation_matrix_jacobian.nrows()) not_eq
              jacobian_targets.target_count() or
          propagation_matrix_jacobian.ncols() not_eq f_grid.size(),
      "Mismatch dimensions on internal matrices of xsec derivatives and frequency");

  // Jacobian overhead START
  /* NOTE:  The calculations below are inefficient and could
              be made much better by using interp in Extract to
              return the derivatives as well. */
  // Jacobian vectors START
  Vector dxsec_temp_dT;
  Vector dxsec_temp_dF;
  Vector dfreq;
  // Jacobian vectors END
  const auto freq_jac = jacobian_targets.find_all<Jacobian::AtmTarget>(
      AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w);
  const auto temp_jac = jacobian_targets.find<Jacobian::AtmTarget>(AtmKey::t);
  const bool do_freq_jac =
      std::ranges::any_of(freq_jac, [](auto& x) { return x.first; });
  const bool do_temp_jac = temp_jac.first;
  const Numeric df       = field_perturbation(freq_jac);
  const Numeric dt       = temp_jac.first;
  if (do_freq_jac) {
    dfreq.resize(f_grid.size());
    dfreq  = f_grid;
    dfreq += df;
    dxsec_temp_dF.resize(f_grid.size());
  }
  if (do_temp_jac) {
    dxsec_temp_dT.resize(f_grid.size());
  }
  // Jacobian overhead END

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  // Allocate a vector with dimension frequencies for constructing our
  // cross-sections before adding them (more efficient to allocate this here
  // outside of the loops)
  Vector xsec_temp(f_grid.size(), 0.);

  ArrayOfString fail_msg;
  bool do_abort = false;
  // Loop over Xsec data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (auto& this_xdata : absorption_xsec_fit_data) {
    if (select_species != SpeciesEnum::Bath and
        this_xdata.Species() != select_species)
      continue;
    if (do_abort) continue;

    const Numeric vmr       = atm_point[this_xdata.Species()];
    const Numeric current_p = force_p < 0 ? atm_point.pressure : force_p;
    const Numeric current_t = force_t < 0 ? atm_point.temperature : force_t;

    // Get the absorption cross sections from the HITRAN data:
    this_xdata.Extract(xsec_temp, f_grid, current_p, current_t);
    if (do_freq_jac) {
      this_xdata.Extract(dxsec_temp_dF, dfreq, current_p, current_t);
    }
    if (do_temp_jac) {
      this_xdata.Extract(dxsec_temp_dT, f_grid, current_p, current_t + dt);
    }

    // Add to result variable:
    Numeric nd = number_density(atm_point.pressure, atm_point.temperature);

    Numeric dnd_dt =
        dnumber_density_dt(atm_point.pressure, atm_point.temperature);
    for (Index f = 0; f < f_grid.size(); f++) {
      propagation_matrix[f].A() += xsec_temp[f] * nd * vmr;

      if (temp_jac.first) {
        const auto iq = temp_jac.second->target_pos;
        propagation_matrix_jacobian[iq, f].A() +=
            ((dxsec_temp_dT[f] - xsec_temp[f]) / dt * nd +
             xsec_temp[f] * dnd_dt) *
            vmr;
      }

      for (auto& j : freq_jac) {
        if (j.first) {
          const auto iq = j.second->target_pos;
          propagation_matrix_jacobian[iq, f].A() +=
              (dxsec_temp_dF[f] - xsec_temp[f]) * nd * vmr / df;
        }
      }

      if (const auto j =
              jacobian_targets.find<Jacobian::AtmTarget>(this_xdata.Species());
          j.first) {
        const auto iq                           = j.second->target_pos;
        propagation_matrix_jacobian[iq, f].A() += xsec_temp[f] * nd * vmr;
      }
    }
  }
}
