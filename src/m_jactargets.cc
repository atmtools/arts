#include <enums.h>
#include <new_jacobian.h>

#include <iterator>

#include "configtypes.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "quantum_numbers.h"

ENUMCLASS(FieldComponent, char, u, v, w)

void jacobian_targetsInit(JacobianTargets& jacobian_targets) {
  jacobian_targets.clear();
}

void jacobian_targetsFinalize(JacobianTargets& jacobian_targets,
                              const AtmField& atm_field) {
  jacobian_targets.finalize(atm_field);
}

void jacobian_targetsAddTemperature(JacobianTargets& jacobian_targets,
                                    const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      Atm::Key::t, d, jacobian_targets.target_count());
}

void jacobian_targetsAddPressure(JacobianTargets& jacobian_targets,
                                 const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      Atm::Key::p, d, jacobian_targets.target_count());
}

void jacobian_targetsAddMagneticField(JacobianTargets& jacobian_targets,
                                      const String& component,
                                      const Numeric& d) {
  using enum FieldComponent;
  switch (toFieldComponentOrThrow(component)) {
    case u:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::mag_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::mag_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::mag_w, d, jacobian_targets.target_count());
      break;
    case FINAL:;
  }
}

void jacobian_targetsAddWindField(JacobianTargets& jacobian_targets,
                                  const String& component,
                                  const Numeric& d) {
  using enum FieldComponent;
  switch (toFieldComponentOrThrow(component)) {
    case u:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::wind_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::wind_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
          Atm::Key::wind_w, d, jacobian_targets.target_count());
      break;
    case FINAL:;
  }
}

void jacobian_targetsAddSpeciesVMR(JacobianTargets& jacobian_targets,
                                   const String& species,
                                   const Numeric& d) {
  Species::Species s = Species::fromShortName(species);
  if (not good_enum(s)) s = Species::toSpeciesOrThrow(species);

  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      s, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const String& species,
    const Numeric& d) {
  const Index i = Species::find_species_index(species);
  ARTS_USER_ERROR_IF(i < 0, "Unknown isotopologue: ", std::quoted(species));

  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      Species::Isotopologues[i], d, jacobian_targets.target_count());
}

void jacobian_targetsAddLineParameter(JacobianTargets& jacobian_targets,
                                      const AbsorptionBands& absorption_bands,
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

  ARTS_USER_ERROR_IF(
      key.line >= band.data.lines.size() or line_index < 0,
      "Line index out of range: ",
      line_index,
      " band has ",
      band.data.lines.size(),
      " absorption lines.");

  if (coefficient.empty()) {
    key.var = lbl::tovariableOrThrow(parameter);
  } else {
    key.ls_var = lbl::line_shape::tovariableOrThrow(parameter);
    key.ls_coeff = lbl::temperature::tocoefficientOrThrow(coefficient);
    key.spec = [&]() {
      auto ptr = std::ranges::find(
          band.data.lines[line_index].ls.single_models,
          [](const String& spec) {
            const auto s = Species::fromShortName(spec);
            return good_enum(s) ? s : Species::toSpeciesOrThrow(spec);
          }(species),
          &lbl::line_shape::species_model::species);
      ARTS_USER_ERROR_IF(
          ptr == band.data.lines[line_index].ls.single_models.end(),
          "No species model for species: ",
          std::quoted(species),
          " in line: ",
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
