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
void abs_xsec_per_speciesAddPredefinedO2MPM2020(ArrayOfMatrix& abs_xsec_per_species,
                                                ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
                                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                                const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                const Vector& f_grid,
                                                const Vector& abs_p,
                                                const Vector& abs_t,
                                                const Matrix& abs_vmrs,
                                                const Verbosity&)
{
  // Forward simulations and their error handling
  ARTS_USER_ERROR_IF (abs_vmrs.ncols() not_eq abs_p.nelem(),
    "Mismatch dimensions on pressure and VMR inputs")
  ARTS_USER_ERROR_IF (abs_t.nelem() not_eq abs_p.nelem(),
    "Mismatch dimensions on pressure and temperature inputs");
  ARTS_USER_ERROR_IF (abs_vmrs.nrows() not_eq abs_species.nelem(),
    "Mismatch dimensions on species and VMR inputs");
  ARTS_USER_ERROR_IF (abs_xsec_per_species.nelem() not_eq abs_species.nelem(),
    "Mismatch dimensions on xsec and VMR inputs");
  ARTS_USER_ERROR_IF (std::any_of(abs_xsec_per_species.cbegin(), abs_xsec_per_species.cend(), 
    [&abs_p](auto x){return x.ncols() not_eq abs_p.nelem();}),
    "Mismatch dimensions on internal matrices of xsec and pressure");
  ARTS_USER_ERROR_IF (std::any_of(abs_xsec_per_species.cbegin(), abs_xsec_per_species.cend(),
    [&f_grid](auto x){return x.nrows() not_eq f_grid.nelem();}),
    "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  if (dabs_xsec_per_species_dx.nelem()) {
    ARTS_USER_ERROR_IF(dabs_xsec_per_species_dx.nelem() not_eq abs_species.nelem(),
      "Mismatch dimensions on species inputs and xsec derivatives");
    ARTS_USER_ERROR_IF (std::any_of(dabs_xsec_per_species_dx.cbegin(), dabs_xsec_per_species_dx.cend(), 
      [&jacobian_quantities](auto x){return x.nelem() not_eq jacobian_quantities.nelem();}),
      "Mismatch dimensions on xsec derivatives and Jacobian grids");
    ARTS_USER_ERROR_IF (std::any_of(dabs_xsec_per_species_dx.cbegin(), dabs_xsec_per_species_dx.cend(),
      [&abs_p](auto x1){return std::any_of(x1.cbegin(), x1.cend(), 
        [&abs_p](auto x2){return x2.ncols() not_eq abs_p.nelem();});}),
      "Mismatch dimensions on internal matrices of xsec derivatives and pressure");
    ARTS_USER_ERROR_IF (std::any_of(dabs_xsec_per_species_dx.cbegin(), dabs_xsec_per_species_dx.cend(),
      [&f_grid](auto x1){return std::any_of(x1.cbegin(), x1.cend(), 
        [&f_grid](auto x2){return x2.nrows() not_eq f_grid.nelem();});}),
      "Mismatch dimensions on internal matrices of xsec derivatives and frequency");
  }
  
  // Positions of important species and VMR of water
  auto o2_mpm2020 =  find_first_species_tag(abs_species, SpeciesTag("O2-MPM2020"));
  auto h2o = find_first_species(abs_species, Species::fromShortName("H2O"));
  auto h2o_vmr = h2o == -1 ? Vector(abs_p.nelem(), 0) : abs_vmrs(h2o, joker);
  ArrayOfMatrix empty(0);
  
  // Perform calculations if there is any oxygen
  if (o2_mpm2020 >= 0 and o2_mpm2020 < abs_xsec_per_species.nelem()) {
    Absorption::PredefinedModel::makarov2020_o2_lines_mpm(abs_xsec_per_species[o2_mpm2020], 
                                                          dabs_xsec_per_species_dx.nelem() ? dabs_xsec_per_species_dx[o2_mpm2020] : empty,
                                                          f_grid, abs_p, abs_t, h2o_vmr, jacobian_quantities);
  }
}
