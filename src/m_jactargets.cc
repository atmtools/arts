#include <enums.h>
#include <jacobian.h>

#include <iterator>

#include "configtypes.h"
#include "debug.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "quantum_numbers.h"

void jacobian_targetsInit(JacobianTargets& jacobian_targets) {
  jacobian_targets.clear();
}

void jacobian_targetsFinalize(JacobianTargets& jacobian_targets,
                              const AtmField& atmospheric_field,
                              const SurfaceField& surface_field,
                              const ArrayOfAbsorptionBand& absorption_bands) {
  jacobian_targets.finalize(atmospheric_field, surface_field, absorption_bands);
}

void jacobian_targetsAddTemperature(JacobianTargets& jacobian_targets,
                                    const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      AtmKey::t, d, jacobian_targets.target_count());
}

void jacobian_targetsAddPressure(JacobianTargets& jacobian_targets,
                                 const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      AtmKey::p, d, jacobian_targets.target_count());
}

void jacobian_targetsAddMagneticField(JacobianTargets& jacobian_targets,
                                      const String& component,
                                      const Numeric& d) {
  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::mag_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::mag_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::mag_w, d, jacobian_targets.target_count());
      break;
  }
}

void jacobian_targetsAddWindField(JacobianTargets& jacobian_targets,
                                  const String& component,
                                  const Numeric& d) {
  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::wind_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::wind_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          AtmKey::wind_w, d, jacobian_targets.target_count());
      break;
  }
}

void jacobian_targetsAddSpeciesVMR(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& species,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      species, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSpeciesVMR(JacobianTargets& jacobian_targets,
                                   const String& species,
                                   const Numeric& d) {
  jacobian_targetsAddSpeciesVMR(jacobian_targets, to<SpeciesEnum>(species), d);
}
#include <iostream>
void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const Numeric& d) {
      std::cout <<"SpeciesIsotope\n";
  ARTS_USER_ERROR_IF(
      std::ranges::none_of(Species::Isotopologues, Cmp::eq(species)),
      "Unknown isotopologue: \"",
      species.FullName(),
      '"');

  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      species, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const String& species,
    const Numeric& d) {
      std::cout <<"String\n";
  jacobian_targetsAddSpeciesIsotopologueRatio(
      jacobian_targets, SpeciesIsotope{species}, d);
}

void jacobian_targetsAddLineParameter(
    JacobianTargets& jacobian_targets,
    const ArrayOfAbsorptionBand& absorption_bands,
    const QuantumIdentifier& qid,
    const Index& line_index,
    const String& parameter,
    const String& species,
    const String& coefficient,
    const Numeric& d) {
  lbl::line_key key{qid};
  key.line = static_cast<Size>(line_index);

  const AbsorptionBand& band = [&]() {
    auto ptr = std::ranges::find(absorption_bands, qid, &lbl::band::key);
    ARTS_USER_ERROR_IF(ptr == absorption_bands.end(),
                       "No band with quantum identifier: ",
                       qid);
    return *ptr;
  }();

  ARTS_USER_ERROR_IF(key.line >= band.data.lines.size() or line_index < 0,
                     "Line index out of range: ",
                     line_index,
                     " band has ",
                     band.data.lines.size(),
                     " absorption lines.");

  if (coefficient.empty()) {
    key.var = to<LineByLineVariable>(parameter);
  } else {
    key.ls_var   = to<LineShapeModelVariable>(parameter);
    key.ls_coeff = to<LineShapeModelCoefficient>(coefficient);
    key.spec     = [&]() {
      auto ptr = std::ranges::find(
          band.data.lines[line_index].ls.single_models,
          [](const String& spec) { return to<SpeciesEnum>(spec); }(species),
          &lbl::line_shape::species_model::species);
      ARTS_USER_ERROR_IF(
          ptr == band.data.lines[line_index].ls.single_models.end(),
          "No species model for species: \"",
          species,
          "\" in line: ",
          band.data.lines[line_index],
          " for quantum identifier: ",
          qid);
      return std::distance(band.data.lines[line_index].ls.single_models.begin(),
                           ptr);
    }();
  }

  jacobian_targets.target<Jacobian::LineTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}
