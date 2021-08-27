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
 * @file   m_fullmodel.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */


#include "predefined_absorption_models.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddPredefinedO2MPM2020(PropagationMatrix& propmat_clearsky,
                                            ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                            const ArrayOfArrayOfSpeciesTag& abs_species,
                                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                                            const Vector& f_grid,
                                            const Numeric& rtp_pressure,
                                            const Numeric& rtp_temperature,
                                            const Vector& rtp_vmr,
                                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Forward simulations and their error handling
  ARTS_USER_ERROR_IF (rtp_vmr.nelem() not_eq abs_species.nelem(),
    "Mismatch dimensions on species and VMR inputs");
  ARTS_USER_ERROR_IF (propmat_clearsky.NumberOfFrequencies() not_eq f_grid.nelem(),
    "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  if (dpropmat_clearsky_dx.nelem()) {
    ARTS_USER_ERROR_IF (dpropmat_clearsky_dx.nelem() not_eq jacobian_quantities.nelem(),
      "Mismatch dimensions on xsec derivatives and Jacobian grids");
    ARTS_USER_ERROR_IF (std::any_of(dpropmat_clearsky_dx.cbegin(), dpropmat_clearsky_dx.cend(),
      [&f_grid](auto& x){return x.NumberOfFrequencies() not_eq f_grid.nelem();}),
      "Mismatch dimensions on internal matrices of xsec derivatives and frequency");
  }
  
  // Perform calculations if there is any oxygen
  if (find_first_species_tag(abs_species, SpeciesTag("O2-MPM2020")) >= 0) {
    const Index h2o = find_first_species(abs_species, Species::fromShortName("H2O"));
    Absorption::PredefinedModel::makarov2020_o2_lines_mpm(propmat_clearsky, dpropmat_clearsky_dx,
                                                          f_grid, rtp_pressure, rtp_temperature, h2o == -1 ? 0.0 : rtp_vmr[h2o], jacobian_quantities);
  }
}
