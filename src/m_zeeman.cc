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

#include "debug.h"
#include "global_data.h"
#include "messages.h"
#include "ppath_struct.h"
#include "propagationmatrix.h"
#include "zeeman.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Vector& rtp_mag,
    const Vector& ppath_los,
    const Index& atmosphere_dim,
    const Index& nlte_do,
    const Index& lbl_checked,
    const Index& manual_zeeman_tag,
    const Numeric& manual_zeeman_magnetic_field_strength,
    const Numeric& manual_zeeman_theta,
    const Numeric& manual_zeeman_eta,
    const Verbosity&) {
  if (abs_lines_per_species.nelem() == 0) return;

  ARTS_USER_ERROR_IF((atmosphere_dim not_eq 3) and (not manual_zeeman_tag),
                     "Only for 3D *atmosphere_dim* or a manual magnetic field")

  ARTS_USER_ERROR_IF((ppath_los.nelem() not_eq 2) and (not manual_zeeman_tag),
                     "Only for 2D *ppath_los* or a manual magnetic field");

  ARTS_USER_ERROR_IF(not lbl_checked,
                     "Please set lbl_checked true to use this function")

  // Change to LOS by radiation
  Vector rtp_los;
  if (not manual_zeeman_tag) mirror_los(rtp_los, ppath_los, atmosphere_dim);

  // Main computations
  zeeman_on_the_fly(propmat_clearsky,
                    nlte_source,
                    dpropmat_clearsky_dx,
                    dnlte_source_dx,
                    abs_species,
                    select_abs_species,
                    jacobian_quantities,
                    abs_lines_per_species,
                    isotopologue_ratios,
                    f_grid,
                    rtp_vmr,
                    rtp_nlte,
                    rtp_mag,
                    rtp_los,
                    rtp_pressure,
                    rtp_temperature,
                    nlte_do,
                    manual_zeeman_tag,
                    manual_zeeman_magnetic_field_strength,
                    manual_zeeman_theta,
                    manual_zeeman_eta);
}

void abs_linesZeemanCoefficients(ArrayOfAbsorptionLines& abs_lines,
                                 const ArrayOfQuantumIdentifier& qid,
                                 const Vector& gs,
                                 const Verbosity&) {
  ARTS_USER_ERROR_IF(qid.nelem() not_eq gs.nelem(),
                     "Inputs not matching in size");
  for (Index i = 0; i < qid.nelem(); i++) {
    const QuantumIdentifier& id = qid[i];
    const Numeric g = gs[i];

    for (AbsorptionLines& band : abs_lines) {
      if (id.isotopologue_index not_eq band.quantumidentity.isotopologue_index)
        continue;

      for (auto& line : band.lines) {
        auto test = id.part_of(band.quantumidentity, line.localquanta);

        if (test.upp) line.zeeman.gu(g);
        if (test.low) line.zeeman.gl(g);
      }
    }
  }
}

void abs_lines_per_speciesZeemanCoefficients(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& qid,
    const Vector& gs,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (Index i = 0; i < qid.nelem(); i++) {
      abs_linesZeemanCoefficients(abs_lines, qid, gs, verbosity);
    }
  }
}

void zeeman_magnetic_fieldFromPath(Matrix& B,
                                   const Ppath& ppath,
                                   const Matrix& ppvar_mag,
                                   const Verbosity&) {
  const Index np = ppath.np;
  ARTS_USER_ERROR_IF(np != ppath.los.nrows(),
                     "Inconsistent number of points in ppath, np: ",
                     np,
                     " ppath.los.nrows(): ",
                     ppath.los.nrows());
  ARTS_USER_ERROR_IF(ppath.los.ncols() != 2,
                     "ppath.los must have 2 columns, have ",
                     ppath.los.ncols());
  ARTS_USER_ERROR_IF(np != ppvar_mag.ncols(),
                     "Inconsistent number of points in ppvar_mag, np: ",
                     np,
                     " ppvar_mag.ncols(): ",
                     ppvar_mag.ncols());
  ARTS_USER_ERROR_IF(ppvar_mag.nrows() != 3,
                     "ppvar_mag must have 3 columns, is ",
                     ppvar_mag.nrows());
  ARTS_USER_ERROR_IF(ppath.dim != 3, "ppath.dim must be 3, is ", ppath.dim);

  B.resize(np, 3);

  for (Index ip = 0; ip < np; ip++) {
    const auto X = Zeeman::FromGrids(ppvar_mag(0, ip),
                                     ppvar_mag(1, ip),
                                     ppvar_mag(2, ip),
                                     Conversion::deg2rad(ppath.los(ip, 0)),
                                     Conversion::deg2rad(ppath.los(ip, 1)));
    B(ip, 0) = X.H;
    B(ip, 1) = Conversion::rad2deg(X.theta);
    B(ip, 2) = Conversion::rad2deg(X.eta);
  }
}
