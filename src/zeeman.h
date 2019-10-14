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
 * @brief Header of Zeeman propagation matrix calculations
 */

#include "abs_species_tags.h"
#include "global_data.h"
#include "linerecord.h"
#include "physics_funcs.h"
#include "quantum.h"
#include "rte.h"

/** Creates a Zeeman ArrayOfArrayOfLineRecord
 * 
 * Sets the Zeeman LineRecord(s) for all the
 * lines.  Sets them by best computations if
 * this is possible or by simplified Hund cases
 * elsewise.  Can be forced to set zeroes instead
 * of real values (e.g., to have them replaced later)
 * 
 * Tests that the computations are sane by 
 * summing up relative strengths.
 * 
 * @param[out] aoaol List of list of lines with set Zeeman effects
 * @param[in]  abs_species as WSV
 * @param[in]  abs_lines_per_species as WSV
 * @param[in]  zero_values Sets Zeeman splitting coefficients to 0
 */
void create_Zeeman_linerecordarrays(
    ArrayOfArrayOfLineRecord& aoaol,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const bool zero_values);

/** Main and only way to compute Zeeman effect
 * 
 * Computes the effect and the derivatives.
 * 
 * Should work in NLTE settings but this is
 * not well-tested
 * 
 * @param[in,out] propmat_clearsky as WSV
 * @param[in,out] nlte_source as WSV
 * @param[in,out] dpropmat_clearsky_dx as WSV
 * @param[in,out] dnlte_dx_source as WSV
 * @param[in,out] nlte_dsource_dx as WSV
 * @param[in]  abs_species as WSV
 * @param[in]  jacobian_quantities as WSV
 * @param[in]  zeeman_linerecord_precalc as WSV
 * @param[in]  isotopologue_ratios as WSV
 * @param[in]  partition_functions as WSV
 * @param[in]  f_grid as WSV
 * @param[in]  rtp_vmr as WSV
 * @param[in]  rtp_nlte as WSV
 * @param[in]  rtp_mag as WSV
 * @param[in]  rtp_los as WSV
 * @param[in]  rtp_pressure as WSV
 * @param[in]  rtp_temperature as WSV
 * @param[in]  manual_zeeman_tag Sets whether the the magnetic field is input manually
 * @param[in]  manual_zeeman_magnetic_field_strength Magnetic field strength
 * @param[in]  manual_zeeman_theta Magnetic field theta angle
 * @param[in]  manual_zeeman_eta Magnetic field eta angle
 */
void zeeman_on_the_fly(
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfStokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_dx_source,
    ArrayOfStokesVector& nlte_dsource_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Vector& f_grid,
    const Vector& rtp_vmr,
    const Vector& rtp_nlte,
    const Vector& rtp_mag,
    const Vector& rtp_los,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const Index& manual_zeeman_tag,
    const Numeric& manual_zeeman_magnetic_field_strength,
    const Numeric& manual_zeeman_theta,
    const Numeric& manual_zeeman_eta);


/** Main and only way to compute Zeeman effect
 * 
 * Computes the effect and the derivatives.
 * 
 * Should work in NLTE settings but this is
 * not well-tested
 * 
 * @param[in,out] propmat_clearsky as WSV
 * @param[in,out] nlte_source as WSV
 * @param[in,out] dpropmat_clearsky_dx as WSV
 * @param[in,out] dnlte_dx_source as WSV
 * @param[in,out] nlte_dsource_dx as WSV
 * @param[in]  abs_species as WSV
 * @param[in]  jacobian_quantities as WSV
 * @param[in]  abs_lines_per_species as WSV
 * @param[in]  isotopologue_ratios as WSV
 * @param[in]  partition_functions as WSV
 * @param[in]  f_grid as WSV
 * @param[in]  rtp_vmr as WSV
 * @param[in]  rtp_nlte as WSV
 * @param[in]  rtp_mag as WSV
 * @param[in]  rtp_los as WSV
 * @param[in]  rtp_pressure as WSV
 * @param[in]  rtp_temperature as WSV
 * @param[in]  manual_zeeman_tag Sets whether the the magnetic field is input manually
 * @param[in]  manual_zeeman_magnetic_field_strength Magnetic field strength
 * @param[in]  manual_zeeman_theta Magnetic field theta angle
 * @param[in]  manual_zeeman_eta Magnetic field eta angle
 */
void zeeman_on_the_fly2(
  ArrayOfPropagationMatrix& propmat_clearsky,
  ArrayOfStokesVector& nlte_source,
  ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
  ArrayOfStokesVector& dnlte_dx_source,
  ArrayOfStokesVector& nlte_dsource_dx,
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
  const Index& manual_tag,
  const Numeric& H0,
  const Numeric& theta0,
  const Numeric& eta0);
