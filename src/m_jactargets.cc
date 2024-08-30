#include <jacobian.h>
#include <enumsFieldComponent.h>

#include <iterator>
#include <limits>

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

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceKey& key,
                                const Numeric& d) {
  jacobian_targets.target<Jacobian::SurfaceTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceTypeTag& key,
                                const Numeric& d) {
  jacobian_targets.target<Jacobian::SurfaceTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfacePropertyTag& key,
                                const Numeric& d) {
  jacobian_targets.target<Jacobian::SurfaceTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const AtmKey& key,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& key,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesIsotope& key,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const QuantumIdentifier& key,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const ParticulatePropertyTag& key,
                                   const Numeric& d) {
  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      key, d, jacobian_targets.target_count());
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

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const Numeric& d) {
  ARTS_USER_ERROR_IF(
      std::ranges::none_of(Species::Isotopologues, Cmp::eq(species)),
      "Unknown isotopologue: \"",
      species.FullName(),
      '"');

  jacobian_targets.target<Jacobian::AtmTarget>().emplace_back(
      species, d, jacobian_targets.target_count());
}
