/**
 * @file   zeeman.cc
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2014-10-14
 * 
 * @brief Header of Zeeman propagation matrix calculations
 */

#include "physics_funcs.h"
#include "quantum_numbers.h"
#include "rte.h"
#include "species_tags.h"

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
 * @param[in,out] dnlte_source_dx as WSV
 * @param[in]  abs_species as WSV
 * @param[in]  select_abs_species as WSV
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
 * @param[in]  nlte_do as WSV
 * @param[in]  manual_zeeman_tag Sets whether the the magnetic field is input manually
 * @param[in]  manual_zeeman_magnetic_field_strength Magnetic field strength
 * @param[in]  manual_zeeman_theta Magnetic field theta angle
 * @param[in]  manual_zeeman_eta Magnetic field eta angle
 */
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
    const Numeric& eta0);
