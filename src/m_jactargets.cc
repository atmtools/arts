#include <enums.h>
#include <new_jacobian.h>

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
