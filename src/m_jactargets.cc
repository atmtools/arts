#include <workspace.h>

#include <ranges>

void jacobian_targetsOff(JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  jacobian_targets.clear();
}

void jacobian_targetsInit(JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  jacobian_targetsOff(jacobian_targets);
}

void jacobian_targetsFinalize(JacobianTargets& jacobian_targets,
                              const AtmField& atmospheric_field,
                              const SurfaceField& surface_field,
                              const AbsorptionBands& absorption_bands,
                              const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  jacobian_targets.finalize(
      atmospheric_field, surface_field, absorption_bands, measurement_sensor);
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfaceKey& key,
                                const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jacobian_targetsAddSurface(JacobianTargets& jacobian_targets,
                                const SurfacePropertyTag& key,
                                const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const AtmKey& key,
                                   const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesEnum& key,
                                   const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const SpeciesIsotope& key,
                                   const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddAtmosphere(JacobianTargets& jacobian_targets,
                                   const QuantumIdentifier& key,
                                   const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{key}, d);
}

void jacobian_targetsAddTemperature(JacobianTargets& jacobian_targets,
                                    const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{AtmKey::t}, d);
}

void jacobian_targetsAddPressure(JacobianTargets& jacobian_targets,
                                 const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{AtmKey::p}, d);
}

void jacobian_targetsAddMagneticField(JacobianTargets& jacobian_targets,
                                      const String& component,
                                      const Numeric& d) {
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{species}, d);
}

void jacobian_targetsAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    const SpeciesIsotope& species,
    const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(AtmKeyVal{species}, d);
}

void jacobian_targetsAddSensorFrequencyPolyFit(
    JacobianTargets& jacobian_targets,
    const ArrayOfSensorObsel& measurement_sensor,
    const Numeric& d,
    const Index& sensor_elem,
    const Index& polyorder) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      polyorder < 0, "Polyorder must be non-negative: {}", polyorder)

  std::vector<std::pair<Index, const void*>> sensor_grid_ptrs;
  for (Size i = 0; i < measurement_sensor.size(); i++) {
    const void* ptr =
        reinterpret_cast<const void*>(measurement_sensor[i].f_grid_ptr().get());
    if (std::ranges::none_of(sensor_grid_ptrs | std::views::values,
                             Cmp::eq(ptr))) {
      sensor_grid_ptrs.emplace_back(i, ptr);
    }
  }

  ARTS_USER_ERROR_IF(
      sensor_grid_ptrs.size() <= static_cast<Size>(sensor_elem),
      "There are only {0} independent sensor frequency grids.  "
      "Cannot select index: {1}, please choose an index less than {0}.",
      sensor_grid_ptrs.size(),
      sensor_elem)

  const Index measurement_elem = sensor_grid_ptrs[sensor_elem].first;

  jacobian_targets.emplace_back(
      SensorKey{.type             = SensorKeyType::f,
                .measurement_elem = measurement_elem,
                .model            = SensorJacobianModelType::PolynomialOffset,
                .polyorder        = polyorder,
                .original_grid    = measurement_sensor[measurement_elem].flat(
                    SensorKeyType::f)},
      d);
}
