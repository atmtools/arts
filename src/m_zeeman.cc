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


#include "new_jacobian.h"
#include "rte.h"
#include "zeeman.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(
    PropmatVector& propmat_clearsky,
    StokvecVector& nlte_source,
    PropmatMatrix& dpropmat_clearsky_dx,
    StokvecMatrix& dnlte_source_dx,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
<<<<<<< HEAD
    const JacobianTargets& jacobian_targets,
=======
<<<<<<< Updated upstream
    const ArrayOfRetrievalQuantity& jacobian_quantities,
>>>>>>> d240b0157 (???)
    const SpeciesIsotopologueRatios& isotopologue_ratios,
=======
    const JacobianTargets& jacobian_targets,
>>>>>>> Stashed changes
    const AtmPoint& atm_point,
    const VibrationalEnergyLevels& nlte_vib_levels,
    const Vector& ppath_los,
    const Index& nlte_do,
    const Index& lbl_checked,
    const Index& manual_zeeman_tag,
    const Numeric& manual_zeeman_magnetic_field_strength,
    const Numeric& manual_zeeman_theta,
    const Numeric& manual_zeeman_eta) {
  if (abs_lines_per_species.size() == 0) return;
  
  ARTS_USER_ERROR_IF((ppath_los.size() not_eq 2) and (not manual_zeeman_tag),
    "Only for 2D *ppath_los* or a manual magnetic field");
  
  ARTS_USER_ERROR_IF(not lbl_checked,
    "Please set lbl_checked true to use this function")

  // Change to LOS by radiation
  Vector rtp_los;
  if (not manual_zeeman_tag) mirror_los(rtp_los, ppath_los);

  // Main computations
  zeeman_on_the_fly(propmat_clearsky,
                    nlte_source,
                    dpropmat_clearsky_dx,
                    dnlte_source_dx,
                    abs_species,
                    select_abs_species,
                    jacobian_targets,
                    abs_lines_per_species,
                    f_grid,
                    atm_point,
                    nlte_vib_levels,
                    rtp_los,
                    nlte_do,
                    manual_zeeman_tag,
                    manual_zeeman_magnetic_field_strength,
                    manual_zeeman_theta,
                    manual_zeeman_eta);
}

void abs_linesZeemanCoefficients(ArrayOfAbsorptionLines& abs_lines,
                                 const ArrayOfQuantumIdentifier& qid,
                                 const Vector& gs) {
  ARTS_USER_ERROR_IF (qid.size() not_eq static_cast<Size>(gs.size()), "Inputs not matching in size");
  for (Size i=0; i<qid.size(); i++) {
    const QuantumIdentifier& id = qid[i];
    const Numeric g = gs[i];
    
    for (AbsorptionLines& band: abs_lines) {
      if (id.isotopologue_index not_eq band.quantumidentity.isotopologue_index) continue;

      for (auto& line: band.lines) {
        auto test = id.part_of(band.quantumidentity, line.localquanta);

        if (test.upp) line.zeeman.gu(g);
        if (test.low) line.zeeman.gl(g);
      }
    }
  }
}

void abs_lines_per_speciesZeemanCoefficients(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                const ArrayOfQuantumIdentifier& qid,
                                                const Vector& gs) {
  for (auto& abs_lines: abs_lines_per_species) {
    for (Size i=0; i<qid.size(); i++) {
      abs_linesZeemanCoefficients(abs_lines, qid, gs);
    }
  }
}
