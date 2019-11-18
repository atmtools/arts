/* Copyright (C) 2012
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
 * @date   2012-08-14
 * 
 * @brief Public methods of ARTS to compute Zeeman effects
 * 
 * Several methods to change and alter and in other way set up
 * Zeeman effect calculations are implemented in this file
 */

#include "global_data.h"
#include "propagationmatrix.h"
#include "zeeman.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfStokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_dx_source,
    ArrayOfStokesVector& nlte_dsource_dx,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Vector& rtp_mag,
    const Vector& ppath_los,
    const Index& atmosphere_dim,
    const Index& lbl_checked,
    const Index& manual_zeeman_tag,
    const Numeric& manual_zeeman_magnetic_field_strength,
    const Numeric& manual_zeeman_theta,
    const Numeric& manual_zeeman_eta,
    const Verbosity&) try {
  if (abs_lines_per_species.nelem() == 0) return;

  if ((atmosphere_dim not_eq 3) and (not manual_zeeman_tag))
    throw "Only for 3D *atmosphere_dim* or a manual magnetic field";
  
  if ((ppath_los.nelem() not_eq 2) and (not manual_zeeman_tag))
    throw "Only for 2D *ppath_los* or a manual magnetic field";
  
  if (not lbl_checked)
    throw "Please set lbl_checked true to use this function";

  // Change to LOS by radiation
  Vector rtp_los;
  if (not manual_zeeman_tag) mirror_los(rtp_los, ppath_los, atmosphere_dim);

  // Main computations
  zeeman_on_the_fly(propmat_clearsky,
                    nlte_source,
                    dpropmat_clearsky_dx,
                    dnlte_dx_source,
                    nlte_dsource_dx,
                    abs_species,
                    jacobian_quantities,
                    abs_lines_per_species,
                    isotopologue_ratios,
                    partition_functions,
                    f_grid,
                    rtp_vmr,
                    rtp_nlte,
                    rtp_mag,
                    rtp_los,
                    rtp_pressure,
                    rtp_temperature,
                    manual_zeeman_tag,
                    manual_zeeman_magnetic_field_strength,
                    manual_zeeman_theta,
                    manual_zeeman_eta);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *propmat_clearskyAddZeeman*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *propmat_clearskyAddZeeman*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}
