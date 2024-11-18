#include <enumsFieldComponent.h>
#include <jacobian.h>

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
                              const AbsorptionBands& absorption_bands) {
  jacobian_targets.finalize(
      atmospheric_field, surface_field, absorption_bands, {});
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceKey& key,
                                const Numeric& d) {
  jacobian_targets.surf().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceTypeTag& key,
                                const Numeric& d) {
  jacobian_targets.surf().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfacePropertyTag& key,
                                const Numeric& d) {
  jacobian_targets.surf().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const AtmKey& key,
                                   const Numeric& d) {
  jacobian_targets.atm().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& key,
                                   const Numeric& d) {
  jacobian_targets.atm().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesIsotope& key,
                                   const Numeric& d) {
  jacobian_targets.atm().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const QuantumIdentifier& key,
                                   const Numeric& d) {
  jacobian_targets.atm().emplace_back(key, d, jacobian_targets.target_count());
}

void jacobian_targetsAddTemperature(JacobianTargets& jacobian_targets,
                                    const Numeric& d) {
  jacobian_targets.atm().emplace_back(
      AtmKey::t, d, jacobian_targets.target_count());
}

void jacobian_targetsAddPressure(JacobianTargets& jacobian_targets,
                                 const Numeric& d) {
  jacobian_targets.atm().emplace_back(
      AtmKey::p, d, jacobian_targets.target_count());
}

void jacobian_targetsAddMagneticField(JacobianTargets& jacobian_targets,
                                      const String& component,
                                      const Numeric& d) {
  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u:
      jacobian_targets.atm().emplace_back(
          AtmKey::mag_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.atm().emplace_back(
          AtmKey::mag_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.atm().emplace_back(
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
      jacobian_targets.atm().emplace_back(
          AtmKey::wind_u, d, jacobian_targets.target_count());
      break;
    case v:
      jacobian_targets.atm().emplace_back(
          AtmKey::wind_v, d, jacobian_targets.target_count());
      break;
    case w:
      jacobian_targets.atm().emplace_back(
          AtmKey::wind_w, d, jacobian_targets.target_count());
      break;
  }
}

void jacobian_targetsAddSpeciesVMR(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& species,
                                   const Numeric& d) {
  jacobian_targets.atm().emplace_back(
      species, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const Numeric& d) {
  ARTS_USER_ERROR_IF(
      std::ranges::none_of(Species::Isotopologues, Cmp::eq(species)),
      "Unknown isotopologue: \"{}\"",
      species.FullName());

  jacobian_targets.atm().emplace_back(
      species, d, jacobian_targets.target_count());
}

void jacobian_targetsAddSensor(JacobianTargets& jacobian_targets,
                               const ArrayOfSensorObsel& measurement_sensor,
                               const SensorKeyType& key,
                               const Numeric& d,
                               const Index& elem) {
  ARTS_USER_ERROR_IF(measurement_sensor.size() <= static_cast<Size>(elem),
                     "Sensor element out of bounds: {}",
                     elem)

  jacobian_targets.sensor().push_back({
      .type       = {.type = key, .elem = elem, .model = SensorJacobianModelType::None},
      .d          = d,
      .target_pos = jacobian_targets.target_count(),
  });
}

void jacobian_targetsAddSensorPolyFit(
    JacobianTargets& jacobian_targets,
    const ArrayOfSensorObsel& measurement_sensor,
    const SensorKeyType& key,
    const Numeric& d,
    const Index& elem,
    const Index& polyorder) {
  ARTS_USER_ERROR_IF(measurement_sensor.size() <= static_cast<Size>(elem),
                     "Sensor element out of bounds: {}",
                     elem)
  ARTS_USER_ERROR_IF(
      polyorder < 0, "Polyorder must be non-negative: {}", polyorder)

  Jacobian::SensorTarget target{
      .type       = {.type          = key,
                     .elem          = elem,
                     .model         = SensorJacobianModelType::PolynomialOffset,
                     .polyorder     = polyorder},
                //     .original_grid = sensor_grid(measurement_sensor[elem], key)},
      .d          = d,
      .target_pos = jacobian_targets.target_count(),
  };
}
