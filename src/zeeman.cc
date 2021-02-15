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
#include "constants.h"
#include "linefunctions.h"
#include "linescaling.h"
#include "species_info.h"

/** Checks if a Propagation Matrix or something similar has good grids */
template <class T>
bool bad_propmat(const Array<T>& main,
                 const Vector& f_grid,
                 const Index sd = 4) noexcept {
  const Index nf = f_grid.nelem();
  for (auto& var : main) {
    const bool bad_stokes = sd not_eq var.StokesDimensions();
    const bool bad_freq = nf not_eq var.NumberOfFrequencies();
    if (bad_freq or bad_stokes) return true;
  }
  return false;
}

/** Checks if abs_species is acceptable */
bool bad_abs_species(const ArrayOfArrayOfSpeciesTag& abs_species) noexcept {
  for (auto& species : abs_species) {
    if (species.nelem()) {
      for (auto& spec : species)
        if (species[0].Species() not_eq spec.Species() or
            species[0].Isotopologue() not_eq spec.Isotopologue() or
            species[0].Type() not_eq spec.Type())
          return true;
    } else
      return true;
  }
  return false;
}

/** Checks for negative values */
template <typename MatpackType> constexpr
bool any_negative(const MatpackType& var) noexcept {
  if (var.empty())
    return false;
  else if (min(var) < 0)
    return true;
  else
    return false;
}

void zeeman_on_the_fly(
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Vector& f_grid,
    const Vector& rtp_vmr,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_mag,
    const Vector& rtp_los,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
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
  ARTS_USER_ERROR_IF(bad_abs_species(abs_species),
    "*abs_species* sub-arrays must have the same species, isotopologue, and type as first sub-array.")
  ARTS_USER_ERROR_IF((rtp_mag.nelem() not_eq 3) and (not manual_tag),
    "Only for 3D *rtp_mag* or a manual magnetic field")
  ARTS_USER_ERROR_IF(rtp_vmr.nelem() not_eq abs_species.nelem(),
    "*rtp_vmr* must match *abs_species*")
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
  ARTS_USER_ERROR_IF(not nq and bad_propmat(dpropmat_clearsky_dx, f_grid),
    "*dpropmat_clearsky_dx* must have Stokes dim 4 and frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(nlte_do and (nq not_eq dnlte_source_dx.nelem()),
    "*dnlte_source_dx* must match derived form of *jacobian_quantities* when non-LTE is on")
  ARTS_USER_ERROR_IF(nlte_do and bad_propmat(dnlte_source_dx, f_grid),
    "*dnlte_source_dx* must have Stokes dim 4 and frequency dim same as *f_grid* when non-LTE is on")
  ARTS_USER_ERROR_IF(any_negative(f_grid), "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(any_negative(rtp_vmr), "Negative VMR (at least one value).")
  ARTS_USER_ERROR_IF(any_negative(rtp_nlte.Data()), "Negative NLTE (at least one value).")
  ARTS_USER_ERROR_IF(rtp_temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(rtp_pressure <= 0, "Non-positive pressure")
  ARTS_USER_ERROR_IF(manual_tag and H0 < 0, "Negative manual magnetic field strength")

  // Pressure information
  const Numeric dnumdens_dmvr = number_density(rtp_pressure, rtp_temperature);
  const Numeric dnumdens_dt_dmvr =
      dnumber_density_dt(rtp_pressure, rtp_temperature);

  // Main compute vectors
  Linefunctions::InternalData scratch(nf, nq), sum(nf, nq);

  // Magnetic field internals and derivatives...
  const auto X =
      manual_tag
          ? Zeeman::FromPreDerived(
                H0, Conversion::deg2rad(theta0), Conversion::deg2rad(eta0))
          : Zeeman::FromGrids(rtp_mag[0],
                              rtp_mag[1],
                              rtp_mag[2],
                              Conversion::deg2rad(rtp_los[0]),
                              Conversion::deg2rad(rtp_los[1]));

  // Polarization
  const auto polarization_scale_data = Zeeman::AllPolarization(X.theta, X.eta);
  const auto polarization_scale_dtheta_data =
      Zeeman::AllPolarization_dtheta(X.theta, X.eta);
  const auto polarization_scale_deta_data =
      Zeeman::AllPolarization_deta(X.theta, X.eta);
      
  // Non-LTE
  const Vector B = nlte_do ? planck(f_grid, rtp_temperature) : Vector(0);
  const Vector dBdT = nlte_do ? dplanck_dt(f_grid, rtp_temperature) : Vector(0);
  const Vector dBdf = nlte_do ? dplanck_df(f_grid, rtp_temperature) : Vector(0);
  const auto eB = MapToEigen(B);
  const auto edBdT = MapToEigen(dBdT);
  const auto edBdf = MapToEigen(dBdf);

  for (auto polar : {Zeeman::Polarization::SigmaMinus,
                     Zeeman::Polarization::Pi,
                     Zeeman::Polarization::SigmaPlus}) {
    auto& pol = Zeeman::SelectPolarization(polarization_scale_data, polar);
    auto& dpol_dtheta =
        Zeeman::SelectPolarization(polarization_scale_dtheta_data, polar);
    auto& dpol_deta =
        Zeeman::SelectPolarization(polarization_scale_deta_data, polar);

    for (Index ispecies = 0; ispecies < ns; ispecies++) {
      
      // Skip it if there are no species or there is no Zeeman
      if (not abs_species[ispecies].nelem() or not is_zeeman(abs_species[ispecies]) or not abs_lines_per_species[ispecies].nelem())
        continue;
      
      for (auto& band : abs_lines_per_species[ispecies]) {
        // Constants for these lines
        const Numeric QT0 = single_partition_function(band.T0(),
                                                      partition_functions.getParamType(band.QuantumIdentity()),
                                                      partition_functions.getParam(band.QuantumIdentity()));
        const Numeric QT = single_partition_function(rtp_temperature,
                                                     partition_functions.getParamType(band.QuantumIdentity()),
                                                     partition_functions.getParam(band.QuantumIdentity()));
        const Numeric dQTdT = dsingle_partition_function_dT(rtp_temperature,
                                                            partition_functions.getParamType(band.QuantumIdentity()),
                                                            partition_functions.getParam(band.QuantumIdentity()));
        const Numeric DC = Linefunctions::DopplerConstant(rtp_temperature, band.SpeciesMass());
        const Numeric dDCdT = Linefunctions::dDopplerConstant_dT(rtp_temperature, DC);
        const Vector line_shape_vmr = band.BroadeningSpeciesVMR(rtp_vmr, abs_species);
        const Numeric numdens = rtp_vmr[ispecies] * dnumdens_dmvr;
        const Numeric dnumdens_dT = rtp_vmr[ispecies] * dnumdens_dt_dmvr;
        const Numeric isotop_ratio = isotopologue_ratios.getIsotopologueRatio(band.QuantumIdentity());
          
        Linefunctions::set_cross_section_of_band(
          scratch,
          sum,
          f_grid,
          band,
          jacobian_quantities,
          line_shape_vmr,
          rtp_nlte,  // This must be turned into a map of some kind...
          rtp_pressure,
          rtp_temperature,
          isotop_ratio,
          X.H,
          DC,
          dDCdT,
          QT,
          dQTdT,
          QT0,
          false,
          true,
          polar);
        
        auto pol_real = pol.attenuation();
        auto pol_imag = pol.dispersion();
        auto abs = propmat_clearsky.Data()(0, 0, joker, joker);

        // Propagation matrix calculations
        MapToEigen(abs).leftCols<4>().noalias() += numdens * sum.F.real() * pol_real;
        MapToEigen(abs).rightCols<3>().noalias() += numdens * sum.F.imag() * pol_imag;

        if (nq) {
          for (Index j = 0; j < nq; j++) {
            if (not propmattype_index(jacobian_quantities, j)) continue;
            const auto& deriv = jacobian_quantities[j];
            Eigen::Map<
                Eigen::Matrix<Numeric, Eigen::Dynamic, 7, Eigen::RowMajor>>
                dabs(dpropmat_clearsky_dx[j].Data().get_c_array(),
                    f_grid.nelem(), 7);

            if (deriv == Jacobian::Atm::Temperature) {
              dabs.leftCols<4>().noalias() +=
                  numdens * sum.dF.col(j).real() * pol_real +
                  dnumdens_dT * sum.F.real() * pol_real;
              dabs.rightCols<3>().noalias() +=
                  numdens * sum.dF.col(j).imag() * pol_imag +
                  dnumdens_dT * sum.F.imag() * pol_imag;
            } else if (deriv == Jacobian::Atm::MagneticU) {
              dabs.leftCols<4>().noalias() +=
                  numdens * X.dH_du * sum.dF.col(j).real() * pol_real +
                  numdens * X.deta_du * sum.F.real() *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_du * sum.F.real() *
                      dpol_dtheta.attenuation();
              dabs.rightCols<3>().noalias() +=
                  numdens * X.dH_du * sum.dF.col(j).imag() * pol_imag +
                  numdens * X.deta_du * sum.F.imag() *
                      dpol_deta.dispersion() +
                  numdens * X.dtheta_du * sum.F.imag() *
                      dpol_dtheta.dispersion();
            } else if (deriv == Jacobian::Atm::MagneticV) {
              dabs.leftCols<4>().noalias() +=
                  numdens * X.dH_dv * sum.dF.col(j).real() * pol_real +
                  numdens * X.deta_dv * sum.F.real() *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_dv * sum.F.real() *
                      dpol_dtheta.attenuation();
              dabs.rightCols<3>().noalias() +=
                  numdens * X.dH_dv * sum.dF.col(j).imag() * pol_imag +
                  numdens * X.deta_dv * sum.F.imag() *
                      dpol_deta.dispersion() +
                  numdens * X.dtheta_dv * sum.F.imag() *
                      dpol_dtheta.dispersion();
            } else if (deriv == Jacobian::Atm::MagneticW) {
              dabs.leftCols<4>().noalias() +=
                  numdens * X.dH_dw * sum.dF.col(j).real() * pol_real +
                  numdens * X.deta_dw * sum.F.real() *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_dw * sum.F.real() *
                      dpol_dtheta.attenuation();
              dabs.rightCols<3>().noalias() +=
                  numdens * X.dH_dw * sum.dF.col(j).imag() * pol_imag +
                  numdens * X.deta_dw * sum.F.imag() *
                      dpol_deta.dispersion() +
                  numdens * X.dtheta_dw * sum.F.imag() *
                      dpol_dtheta.dispersion();
            } else if (deriv == Jacobian::Line::VMR and
                      deriv.QuantumIdentity().In(band.QuantumIdentity())) {
              dabs.leftCols<4>().noalias() +=
                  numdens * sum.dF.col(j).real() * pol_real +
                  dnumdens_dmvr * sum.F.real() * pol_real;
              dabs.rightCols<3>().noalias() +=
                  numdens * sum.dF.col(j).imag() * pol_imag +
                  dnumdens_dmvr * sum.F.imag() * pol_imag;
            } else if (deriv == abs_species[ispecies]) {
              dabs.leftCols<4>().noalias() += numdens * sum.F.real() * pol_real;
              dabs.rightCols<3>().noalias() += numdens * sum.F.imag() * pol_imag;
              
            } else {
              dabs.leftCols<4>().noalias() +=
                  numdens * sum.dF.col(j).real() * pol_real;
              dabs.rightCols<3>().noalias() +=
                  numdens * sum.dF.col(j).imag() * pol_imag;
            }
          }
        }

        // Source vector calculations
        if (nlte_do) {
          auto nlte_src = nlte_source.Data()(0, 0, joker, joker);

          MapToEigen(nlte_src)
              .leftCols<4>()
              .noalias() += numdens * eB.cwiseProduct(sum.N.real()) * pol_real;

          for (Index j = 0; j < nq; j++) {
            if (not propmattype_index(jacobian_quantities, j)) continue;
            const auto& deriv = jacobian_quantities[j];

            Eigen::Map<
                Eigen::Matrix<Numeric, Eigen::Dynamic, 4, Eigen::RowMajor>>
                dnlte_dx_src(dnlte_source_dx[j].Data().get_c_array(),
                            f_grid.nelem(), 4),
                nlte_dsrc_dx(dnlte_source_dx[j].Data().get_c_array(),
                            f_grid.nelem(), 4);

            if (deriv == Jacobian::Atm::Temperature) {
              dnlte_dx_src.noalias() +=
                  dnumdens_dT * eB.cwiseProduct(sum.N.real()) * pol_real +
                  numdens * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real;

              nlte_dsrc_dx.noalias() +=
                  numdens * edBdT.cwiseProduct(sum.N.real()) * pol_real;
            } else if (deriv.Target().isWind()) {
              dnlte_dx_src.noalias() +=
              dnumdens_dT * eB.cwiseProduct(sum.N.real()) * pol_real +
              numdens * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real;
              
              nlte_dsrc_dx.noalias() +=
              numdens * edBdf.cwiseProduct(sum.N.real()) * pol_real;
            } else if (deriv == Jacobian::Atm::MagneticU) {
              dnlte_dx_src.noalias() +=
                  numdens * X.dH_du * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real +
                  numdens * X.deta_du * eB.cwiseProduct(sum.N.real()) *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_du * eB.cwiseProduct(sum.N.real()) *
                      dpol_dtheta.attenuation();
            } else if (deriv == Jacobian::Atm::MagneticV) {
              dnlte_dx_src.noalias() +=
                  numdens * X.dH_dv * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real +
                  numdens * X.deta_dv * eB.cwiseProduct(sum.N.real()) *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_dv * eB.cwiseProduct(sum.N.real()) *
                      dpol_dtheta.attenuation();
            } else if (deriv == Jacobian::Atm::MagneticW) {
              dnlte_dx_src.noalias() +=
                  numdens * X.dH_dw * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real +
                  numdens * X.deta_dw * eB.cwiseProduct(sum.N.real()) *
                      dpol_deta.attenuation() +
                  numdens * X.dtheta_dw * eB.cwiseProduct(sum.N.real()) *
                      dpol_dtheta.attenuation();
            } else if (deriv == Jacobian::Line::VMR and
                     deriv.QuantumIdentity().In(band.QuantumIdentity())) {
              dnlte_dx_src.noalias() +=
                  dnumdens_dmvr * eB.cwiseProduct(sum.N.real()) * pol_real +
                  numdens * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real;
            } else if (deriv == abs_species[ispecies]) {
              dnlte_dx_src.noalias() += numdens * eB.cwiseProduct(sum.N.real()) * pol_real;
            } else {
              dnlte_dx_src.noalias() +=
                  numdens * eB.cwiseProduct(sum.dN.col(j).real()) * pol_real;
            }
          }
        }
      }
    }
  }
}
