#include <enumsFieldComponent.h>
#include <jacobian.h>

void jacobian_targetsInit(JacobianTargets& jacobian_targets) {
  jacobian_targets.clear();
}

void jacobian_targetsFinalize(JacobianTargets& jacobian_targets,
                              const AtmField& atmospheric_field,
                              const SurfaceField& surface_field,
                              const AbsorptionBands& absorption_bands,
                              const ArrayOfSensorObsel& measurement_sensor) {
  jacobian_targets.finalize(
      atmospheric_field, surface_field, absorption_bands, measurement_sensor);
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceKey& key,
                                const Numeric& d) {
  jacobian_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceTypeTag& key,
                                const Numeric& d) {
  jacobian_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfacePropertyTag& key,
                                const Numeric& d) {
  jacobian_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const AtmKey& key,
                                   const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& key,
                                   const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesIsotope& key,
                                   const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const QuantumIdentifier& key,
                                   const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddTemperature(JacobianTargets& jacobian_targets,
                                    const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{AtmKey::t}, d);
}

void jacobian_targetsAddPressure(JacobianTargets& jacobian_targets,
                                 const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{AtmKey::p}, d);
}

void jacobian_targetsAddMagneticField(JacobianTargets& jacobian_targets,
                                      const String& component,
                                      const Numeric& d) {
  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::mag_u}, d); break;
    case v: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::mag_v}, d); break;
    case w: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::mag_w}, d); break;
  }
}

void jacobian_targetsAddWindField(JacobianTargets& jacobian_targets,
                                  const String& component,
                                  const Numeric& d) {
  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::wind_u}, d); break;
    case v: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::wind_v}, d); break;
    case w: jacobian_targets.emplace_back(AtmKeyVal{AtmKey::wind_w}, d); break;
  }
}

void jacobian_targetsAddSpeciesVMR(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& species,
                                   const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{species}, d);
}

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const Numeric& d) {
  jacobian_targets.emplace_back(AtmKeyVal{species}, d);
}

void jacobian_targetsAddSensorFrequencyPolyFit(
    JacobianTargets& jacobian_targets,
    const ArrayOfSensorObsel& measurement_sensor,
    const Numeric& d,
    const Index& elem,
    const Index& polyorder) {
  ARTS_USER_ERROR_IF(measurement_sensor.size() <= static_cast<Size>(elem),
                     "Sensor element out of bounds: {}",
                     elem)
  ARTS_USER_ERROR_IF(
      polyorder < 0, "Polyorder must be non-negative: {}", polyorder)

  jacobian_targets.emplace_back(
      SensorKey{
          .type          = SensorKeyType::f,
          .elem          = elem,
          .model         = SensorJacobianModelType::PolynomialOffset,
          .polyorder     = polyorder,
          .original_grid = measurement_sensor[elem].flat(SensorKeyType::f)},
      d);
}
