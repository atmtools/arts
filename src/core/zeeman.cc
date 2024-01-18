/**
 * @file   zeeman.cc
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2014-10-14
 * 
 * @brief Implementations of Zeeman propagation matrix calculations
 * 
 * This file implements Zeeman propagation matrix calculations while
 * also computing the derivatives that might be interesting for 
 * later Jacobian deductions.
 */

#include "zeeman.h"

#include "arts_conversions.h"
#include "lineshape.h"
#include "jacobian.h"
#include "nlte.h"

#include <math_funcs.h>

void zeeman_on_the_fly(
    PropmatVector& propmat_clearsky,
    StokvecVector& nlte_source,
    PropmatMatrix& dpropmat_clearsky_dx,
    StokvecMatrix& dnlte_source_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const JacobianTargets& jacobian_targets,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const VibrationalEnergyLevels& nlte_vib_energies,
    const Vector& rtp_los,
    const Index& nlte_do,
    const Index& manual_tag,
    const Numeric& H0,
    const Numeric& theta0,
    const Numeric& eta0) {
  // Size of problem
  const Index nf = f_grid.size();
  const Index nq = jacobian_targets.target_count();
  const Index ns = abs_species.size();

  // Possible things that can go wrong in this code (excluding line parameters)
  check_abs_species(abs_species);
  ARTS_USER_ERROR_IF(propmat_clearsky.size() not_eq nf,
                     "*f_grid* must match *propmat_clearsky*")
  ARTS_USER_ERROR_IF(nlte_source.size() not_eq nf,
                     "*f_grid* must match *nlte_source*")
  ARTS_USER_ERROR_IF(
      nq not_eq dpropmat_clearsky_dx.nrows() or
          nf not_eq dpropmat_clearsky_dx.ncols(),
      "*dpropmat_clearsky_dx* must match derived form of *jacobian_quantities* times the lenght of *f_grid*")
  ARTS_USER_ERROR_IF(
      nlte_do and (nq not_eq dnlte_source_dx.nrows() or
                   nf not_eq dnlte_source_dx.ncols()),
      "*dnlte_source_dx* must match derived form of *jacobian_quantities* when non-LTE is on")
  ARTS_USER_ERROR_IF(any_negative(f_grid),
                     "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(atm_point.temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(atm_point.pressure <= 0, "Non-positive pressure")
  ARTS_USER_ERROR_IF(manual_tag and H0 < 0,
                     "Negative manual magnetic field strength")

  // Magnetic field internals and derivatives...
  const auto X =
      manual_tag
          ? Zeeman::FromPreDerived(
                H0, Conversion::deg2rad(theta0), Conversion::deg2rad(eta0))
          : Zeeman::FromGrids(atm_point.mag[0],
                              atm_point.mag[1],
                              atm_point.mag[2],
                              Conversion::deg2rad(rtp_los[0]),
                              Conversion::deg2rad(rtp_los[1]));

  // Polarization
  const auto polarization_scale_data = Zeeman::AllPolarization(X.theta, X.eta);
  const auto polarization_scale_dtheta_data =
      Zeeman::AllPolarization_dtheta(X.theta, X.eta);
  const auto polarization_scale_deta_data =
      Zeeman::AllPolarization_deta(X.theta, X.eta);

  // Deal with sparse computational grid
  const Vector f_grid_sparse(0);
  const Numeric sparse_limit = 0;

  for (auto polar : {Zeeman::Polarization::SigmaMinus,
                     Zeeman::Polarization::Pi,
                     Zeeman::Polarization::SigmaPlus}) {
    // Calculations data
    LineShape::ComputeData com(f_grid, jacobian_targets, nlte_do);
    LineShape::ComputeData sparse_com(f_grid_sparse, jacobian_targets, nlte_do);

    auto& pol = Zeeman::SelectPolarization(polarization_scale_data, polar);
    auto& dpol_dtheta =
        Zeeman::SelectPolarization(polarization_scale_dtheta_data, polar);
    auto& dpol_deta =
        Zeeman::SelectPolarization(polarization_scale_deta_data, polar);

    for (Index ispecies = 0; ispecies < ns; ispecies++) {
      // Skip it if there are no species or there is no Zeeman
      if (not abs_species[ispecies].size() or
          not abs_species[ispecies].Zeeman() or
          not abs_lines_per_species[ispecies].size())
        continue;
      if (select_abs_species.size() and
          select_abs_species not_eq abs_species[ispecies])
        continue;

      for (auto& band : abs_lines_per_species[ispecies]) {
        LineShape::compute(com,
                           sparse_com,
                           band,
                           jacobian_targets,
                           atm_point.is_lte()
                               ? std::pair{0., 0.}
                               : atm_point.levels(band.quantumidentity),
                           nlte_vib_energies,
                           band.BroadeningSpeciesVMR(atm_point),
                           abs_species[ispecies],
                           atm_point[band.Species()],
                           atm_point[band.Isotopologue()],
                           atm_point.pressure,
                           atm_point.temperature,
                           X.H,
                           sparse_limit,
                           polar,
                           Options::LblSpeedup::None,
                           false);
      }
    }

    // Sum up the propagation matrix
    Zeeman::sum_propmat(propmat_clearsky, com.F, pol);

    // Sum up the Jacobian
    for (auto& atm : jacobian_targets.atm()) {
      const auto j = atm.target_pos;

      if (atm.type == Atm::Key::mag_u) {
        Zeeman::dsum_propmat(dpropmat_clearsky_dx[j],
                             com.F,
                             com.dF(joker, j),
                             pol,
                             dpol_dtheta,
                             dpol_deta,
                             X.dH_du,
                             X.dtheta_du,
                             X.deta_du);
      } else if (atm.type == Atm::Key::mag_v) {
        Zeeman::dsum_propmat(dpropmat_clearsky_dx[j],
                             com.F,
                             com.dF(joker, j),
                             pol,
                             dpol_dtheta,
                             dpol_deta,
                             X.dH_dv,
                             X.dtheta_dv,
                             X.deta_dv);
      } else if (atm.type == Atm::Key::mag_w) {
        Zeeman::dsum_propmat(dpropmat_clearsky_dx[j],
                             com.F,
                             com.dF(joker, j),
                             pol,
                             dpol_dtheta,
                             dpol_deta,
                             X.dH_dw,
                             X.dtheta_dw,
                             X.deta_dw);
      } else {
        Zeeman::sum_propmat(dpropmat_clearsky_dx[j], com.dF(joker, j), pol);
      }
    }

    if (nlte_do) {
      // Sum up the source vector
      Zeeman::sum_stokvec(nlte_source, com.N, pol);

      // Sum up the Jacobian
      for (auto& atm : jacobian_targets.atm()) {
        const auto j = atm.target_pos;

        if (atm.type == Atm::Key::mag_u) {
          Zeeman::dsum_stokvec(dnlte_source_dx[j],
                               com.N,
                               com.dN(joker, j),
                               pol,
                               dpol_dtheta,
                               dpol_deta,
                               X.dH_du,
                               X.dtheta_du,
                               X.deta_du);
        } else if (atm.type == Atm::Key::mag_v) {
          Zeeman::dsum_stokvec(dnlte_source_dx[j],
                               com.N,
                               com.dN(joker, j),
                               pol,
                               dpol_dtheta,
                               dpol_deta,
                               X.dH_dv,
                               X.dtheta_dv,
                               X.deta_dv);
        } else if (atm.type == Atm::Key::mag_w) {
          Zeeman::dsum_stokvec(dnlte_source_dx[j],
                               com.N,
                               com.dN(joker, j),
                               pol,
                               dpol_dtheta,
                               dpol_deta,
                               X.dH_dw,
                               X.dtheta_dw,
                               X.deta_dw);
        } else {
          Zeeman::sum_stokvec(dnlte_source_dx[j], com.dN(joker, j), pol);
        }
      }

      // FIXME: Add Jacobian for line parameters
    }
  }
}
