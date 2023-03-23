/* Copyright (C) 2014
   Richard Larsson <ric.larsson@gmail.com>

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
#include "linescaling.h"
#include "lineshape.h"
#include "nlte.h"
#include "species_info.h"

void zeeman_on_the_fly(
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const VibrationalEnergyLevels& nlte_vib_levels,
    const Vector& rtp_los,
    const Index& nlte_do,
    const Index& manual_tag,
    const Numeric& H0,
    const Numeric& theta0,
    const Numeric& eta0) {

  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();
  const Index ns = abs_species.nelem();

  // Possible things that can go wrong in this code (excluding line parameters)
  check_abs_species(abs_species);
  ARTS_USER_ERROR_IF(propmat_clearsky.NumberOfFrequencies() not_eq nf,
    "*f_grid* must match *propmat_clearsky*")
  ARTS_USER_ERROR_IF(propmat_clearsky.StokesDimensions() not_eq 4,
    "*propmat_clearsky* must have *stokes_dim* 4")
  ARTS_USER_ERROR_IF(nlte_source.NumberOfFrequencies() not_eq nf,
    "*f_grid* must match *nlte_source*")
  ARTS_USER_ERROR_IF(nlte_source.StokesDimensions() not_eq 4,
    "*nlte_source* must have *stokes_dim* 4")
  ARTS_USER_ERROR_IF(not nq and (nq not_eq dpropmat_clearsky_dx.nelem()),
    "*dpropmat_clearsky_dx* must match derived form of *jacobian_quantities*")
  ARTS_USER_ERROR_IF(not nq and bad_propmat(dpropmat_clearsky_dx, f_grid, 4),
    "*dpropmat_clearsky_dx* must have Stokes dim 4 and frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(nlte_do and (nq not_eq dnlte_source_dx.nelem()),
    "*dnlte_source_dx* must match derived form of *jacobian_quantities* when non-LTE is on")
  ARTS_USER_ERROR_IF(nlte_do and bad_propmat(dnlte_source_dx, f_grid, 4),
    "*dnlte_source_dx* must have Stokes dim 4 and frequency dim same as *f_grid* when non-LTE is on")
  ARTS_USER_ERROR_IF(any_negative(f_grid), "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(atm_point.temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(atm_point.pressure <= 0, "Non-positive pressure")
  ARTS_USER_ERROR_IF(manual_tag and H0 < 0, "Negative manual magnetic field strength")

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
    LineShape::ComputeData com(f_grid, jacobian_quantities, nlte_do);
    LineShape::ComputeData sparse_com(f_grid_sparse, jacobian_quantities, nlte_do);
    
    auto& pol = Zeeman::SelectPolarization(polarization_scale_data, polar);
    auto& dpol_dtheta =
        Zeeman::SelectPolarization(polarization_scale_dtheta_data, polar);
    auto& dpol_deta =
        Zeeman::SelectPolarization(polarization_scale_deta_data, polar);

    for (Index ispecies = 0; ispecies < ns; ispecies++) {
      // Skip it if there are no species or there is no Zeeman
      if (not abs_species[ispecies].nelem() or
          not abs_species[ispecies].Zeeman() or
          not abs_lines_per_species[ispecies].nelem())
        continue;
      if (select_abs_species.nelem() and
          select_abs_species not_eq abs_species[ispecies])
        continue;

      for (auto& band : abs_lines_per_species[ispecies]) {
        LineShape::compute(
            com, sparse_com, band, jacobian_quantities,
            atm_point.is_lte() ? std::pair{0., 0.}
                               : atm_point.levels(band.quantumidentity),
            nlte_vib_levels, band.BroadeningSpeciesVMR(atm_point),
            abs_species[ispecies], atm_point[abs_species[ispecies]],
            isotopologue_ratios[band.Isotopologue()], atm_point.pressure,
            atm_point.temperature, X.H, sparse_limit, polar,
            Options::LblSpeedup::None, false);
      }
    }
      
    // Sum up the propagation matrix
    Zeeman::sum(propmat_clearsky, com.F, pol);
    
    // Sum up the Jacobian
    for (Index j=0; j<nq; j++) {
      auto& deriv = jacobian_quantities[j];
      
      if (not deriv.propmattype()) continue;
      
      if (deriv == Jacobian::Atm::MagneticU) {
        Zeeman::dsum(dpropmat_clearsky_dx[j], com.F, com.dF(joker, j),
                      pol, dpol_dtheta, dpol_deta,
                      X.dH_du, X.dtheta_du, X.deta_du);
      } else if (deriv == Jacobian::Atm::MagneticV) {
        Zeeman::dsum(dpropmat_clearsky_dx[j], com.F, com.dF(joker, j),
                      pol, dpol_dtheta, dpol_deta,
                      X.dH_dv, X.dtheta_dv, X.deta_dv);
      } else if (deriv == Jacobian::Atm::MagneticW) {
        Zeeman::dsum(dpropmat_clearsky_dx[j], com.F, com.dF(joker, j),
                      pol, dpol_dtheta, dpol_deta,
                      X.dH_dw, X.dtheta_dw, X.deta_dw);
      } else {
        Zeeman::sum(dpropmat_clearsky_dx[j], com.dF(joker, j), pol);
      }
    }
    
    if (nlte_do) {
      // Sum up the source vector
      Zeeman::sum(nlte_source, com.N, pol, false);
      
      // Sum up the Jacobian
      for (Index j=0; j<nq; j++) {
        auto& deriv = jacobian_quantities[j];
        
        if (not deriv.propmattype()) continue;
        
        if (deriv == Jacobian::Atm::MagneticU) {
          Zeeman::dsum(dnlte_source_dx[j], com.N, com.dN(joker, j),
                        pol, dpol_dtheta, dpol_deta,
                        X.dH_du, X.dtheta_du, X.deta_du, false);
        } else if (deriv == Jacobian::Atm::MagneticV) {
          Zeeman::dsum(dnlte_source_dx[j], com.N, com.dN(joker, j),
                        pol, dpol_dtheta, dpol_deta,
                        X.dH_dv, X.dtheta_dv, X.deta_dv, false);
        } else if (deriv == Jacobian::Atm::MagneticW) {
          Zeeman::dsum(dnlte_source_dx[j], com.N, com.dN(joker, j),
                        pol, dpol_dtheta, dpol_deta,
                        X.dH_dw, X.dtheta_dw, X.deta_dw, false);
        } else {
          Zeeman::sum(dnlte_source_dx[j], com.dN(joker, j), pol, false);
        }
      }
    }
  }
}
