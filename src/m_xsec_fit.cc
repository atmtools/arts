/*!
  \file   m_hitran_xsec.cc
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2021-02-23

  \brief  Workspace methods for HITRAN absorption cross section data.

*/

#include <workspace.h>
#include "xsec_fit.h"
#include "jacobian.h"
#include "m_xml.h"
#include "physics_funcs.h"
#include "arts_omp.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadXsecData(ArrayOfXsecRecord& xsec_fit_data,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const String& basename) {
  // Build a set of species indices. Duplicates are ignored.
  std::set<Species::Species> unique_species;
  for (auto& asp : abs_species) {
    for (auto& sp : asp) {
      if (sp.Type() == Species::TagType::XsecFit) {
        unique_species.insert(sp.Spec());
      }
    }
  }

  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }

  // Read xsec data for all active species and collect them in xsec_fit_data
  xsec_fit_data.clear();
  for (auto& species_name : unique_species) {
    XsecRecord xsec_coeffs;
    const String filename{tmpbasename +
                          String(Species::toShortName(species_name)) + ".xml"};

    try {
      xml_read_from_file(filename, xsec_coeffs);

      xsec_fit_data.push_back(xsec_coeffs);
    } catch (const std::exception& e) {
      ARTS_USER_ERROR(
          "Error reading coefficients file:\n", filename, "\n", e.what());
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddXsecFit(  // WS Output:
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const ArrayOfXsecRecord& xsec_fit_data,
    const Numeric& force_p,
    const Numeric& force_t) {
  // Forward simulations and their error handling
  ARTS_USER_ERROR_IF(
      propmat_clearsky.nelem() not_eq f_grid.nelem(),
      "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  ARTS_USER_ERROR_IF(
    dpropmat_clearsky_dx.nrows() not_eq jacobian_quantities.nelem()
    or
    dpropmat_clearsky_dx.ncols() not_eq f_grid.nelem(),
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
    if (select_abs_species.nelem() and abs_species[i] not_eq select_abs_species)
      continue;

    const Numeric vmr = atm_point[abs_species[i]];

    for (Index s = 0; s < abs_species[i].nelem(); s++) {
      const SpeciesTag& this_species = abs_species[i][s];

      // Check if this is a HITRAN cross section tag
      if (this_species.Type() != Species::TagType::XsecFit) continue;

      Index this_xdata_index =
          hitran_xsec_get_index(xsec_fit_data, this_species.Spec());
      ARTS_USER_ERROR_IF(this_xdata_index < 0,
                         "Cross-section species ",
                         this_species.Name(),
                         " not found in *xsec_fit_data*.")
      const XsecRecord& this_xdata = xsec_fit_data[this_xdata_index];
      // ArrayOfMatrix& this_dxsec = do_jac ? dpropmat_clearsky_dx[i] : empty;

      if (do_abort) continue;
      const Numeric current_p = force_p < 0 ? atm_point.pressure : force_p;
      const Numeric current_t = force_t < 0 ? atm_point.temperature : force_t;

      // Get the absorption cross sections from the HITRAN data:
      this_xdata.Extract(xsec_temp, f_grid, current_p, current_t);
      if (do_freq_jac) {
        this_xdata.Extract(
            dxsec_temp_dF, dfreq, current_p, current_t);
      }
      if (do_temp_jac) {
        this_xdata.Extract(
            dxsec_temp_dT, f_grid, current_p, current_t + dt);
      }
    }

    // Add to result variable:
    Numeric nd = number_density(atm_point.pressure, atm_point.temperature);
    if (!do_jac) {
      xsec_temp *= nd * vmr;
      for (Index iv=0; iv<f_grid.nelem(); iv++) propmat_clearsky[iv].A() += xsec_temp[iv];
    } else {
      Numeric dnd_dt = dnumber_density_dt(atm_point.pressure, atm_point.temperature);
      for (Index f = 0; f < f_grid.nelem(); f++) {
        propmat_clearsky[f].A() += xsec_temp[f] * nd * vmr;
        for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
          const auto& deriv = jacobian_quantities[iq];

          if (!deriv.propmattype()) continue;

          if (is_frequency_parameter(deriv)) {
            dpropmat_clearsky_dx(iq, f).A() +=
                (dxsec_temp_dF[f] - xsec_temp[f]) * nd * vmr / df;
          } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR ||
                     deriv == Jacobian::Line::VMR) {
            if (species_match(deriv, abs_species[i])) {
              dpropmat_clearsky_dx(iq, f).A() +=
                  xsec_temp[f] * nd * vmr;
            }
          } else if (deriv == Jacobian::Atm::Temperature) {
            dpropmat_clearsky_dx(iq, f).A() +=
                ((dxsec_temp_dT[f] - xsec_temp[f]) / dt * nd +
                 xsec_temp[f] * dnd_dt) *
                vmr;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddXsecFit(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfArrayOfMatrix& /* dabs_xsec_per_species_dx */,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& /* jacobian_quantities */,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const ArrayOfXsecRecord& xsec_fit_data,
    const Numeric& force_p,
    const Numeric& force_t) {
  {
    // Check that all parameters that should have the number of tag
    // groups as a dimension are consistent:
    const Index n_tgs = abs_species.nelem();
    const Index n_xsec = abs_xsec_per_species.nelem();

    ARTS_USER_ERROR_IF(
        n_tgs != n_xsec,
        "The following variables must all have the same dimension:\n"
        "abs_species:          ",
        n_tgs,
        "\n"
        "abs_xsec_per_species: ",
        n_xsec)
  }

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  {
    // Check that all parameters that should have the the dimension of p_grid
    // are consistent:
    const Index n_p = abs_p.nelem();
    const Index n_t = abs_t.nelem();

    ARTS_USER_ERROR_IF(
        n_p != n_t,
        "The following variables must all have the same dimension:\n"
        "abs_p:          ",
        n_p,
        "\n"
        "abs_t:          ",
        n_t)
  }

  // Allocate a vector with dimension frequencies for constructing our
  // cross-sections before adding them (more efficient to allocate this here
  // outside of the loops)
  Vector xsec_temp(f_grid.nelem(), 0.);

  ArrayOfString fail_msg;
  bool do_abort = false;

  // Loop over Xsec data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (Index ii = 0; ii < abs_species_active.nelem(); ii++) {
    const Index i = abs_species_active[ii];

    for (Index s = 0; s < abs_species[i].nelem(); s++) {
      const SpeciesTag& this_species = abs_species[i][s];

      // Check if this is a HITRAN cross section tag
      if (this_species.Type() != Species::TagType::XsecFit) continue;

      Index this_xdata_index =
          hitran_xsec_get_index(xsec_fit_data, this_species.Spec());
      ARTS_USER_ERROR_IF(this_xdata_index < 0,
                         "Cross-section species ",
                         this_species.Name(),
                         " not found in *xsec_fit_data*.")
      const XsecRecord& this_xdata = xsec_fit_data[this_xdata_index];
      Matrix& this_xsec = abs_xsec_per_species[i];

      // Loop over pressure:
#pragma omp parallel for if (!arts_omp_in_parallel() && abs_p.nelem() >= 1) \
    firstprivate(xsec_temp)
      for (Index ip = 0; ip < abs_p.nelem(); ip++) {
        if (do_abort) continue;
        const Numeric current_p = force_p < 0 ? abs_p[ip] : force_p;
        const Numeric current_t = force_t < 0 ? abs_t[ip] : force_t;

        // Get the absorption cross sections from the HITRAN data:
        try {
          this_xdata.Extract(xsec_temp,
                             f_grid,
                             current_p,
                             current_t);
        } catch (std::runtime_error& e) {
          std::ostringstream os;
          os << "Problem with HITRAN cross section species "
             << this_species.Name() << " at pressure level " << ip << " ("
             << abs_p[ip] / 100. << " hPa):\n"
             << e.what() << "\n";
#pragma omp critical(abs_xsec_per_speciesAddXsecFit)
          {
            do_abort = true;
            fail_msg.push_back(os.str());
          }
        }
        this_xsec(joker, ip) += xsec_temp;
      }
    }
  }

  if (do_abort) {
    ARTS_USER_ERROR(fail_msg);
  }
}
