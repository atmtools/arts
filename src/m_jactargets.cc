#include <jac_polyfit.h>
#include <workspace.h>

void jacobian_targetsOff(JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  jacobian_targets.clear();
}

void jacobian_targetsConditionalClear(JacobianTargets& jacobian_targets,
                                      const Index& do_jacobian) {
  ARTS_TIME_REPORT

  if (do_jacobian == 0) {
    jacobian_targetsOff(jacobian_targets);
  }
}

void jacobian_targetsInit(JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  jacobian_targetsOff(jacobian_targets);
}

void jacobian_targetsFinalize(JacobianTargets& jacobian_targets,
                              const AtmField& atmospheric_field,
                              const SurfaceField& surface_field,
                              const SubsurfaceField& subsurface_field,
                              const AbsorptionBands& absorption_bands,
                              const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  jacobian_targets.finalize(atmospheric_field,
                            surface_field,
                            subsurface_field,
                            absorption_bands,
                            measurement_sensor);
}

void jacobian_targetsAddSubsurface(JacobianTargets& jacobian_targets,
                                   const SubsurfaceKey& key,
                                   const Numeric& d) {
  ARTS_TIME_REPORT

  jacobian_targets.emplace_back(SubsurfaceKeyVal{key}, d);
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
                                   const QuantumLevelIdentifier& key,
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

void jacobian_targetsAddOverlappingMagneticField(
    JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  const auto f = stdr::find_if(jacobian_targets.atm, Jacobian::is_mag);

  ARTS_USER_ERROR_IF(f == jacobian_targets.atm.end(),
                     "Cannot find original magnetic field for component");

  const AtmKey component = std::get<AtmKey>(f->type);
  const Size N           = jacobian_targets.target_count();

  jacobian_targets.atm.emplace_back(*f);
  jacobian_targets.atm.back().overlap     = true;
  jacobian_targets.atm.back().overlap_key = component;
  jacobian_targets.atm.back().target_pos  = N;
  switch (component) {
    using enum AtmKey;
    case mag_u: jacobian_targets.atm.back().type = mag_v; break;
    case mag_v:
    case mag_w: jacobian_targets.atm.back().type = mag_u; break;
    default:    throw std::out_of_range("Invalid state");
  }

  jacobian_targets.atm.emplace_back(*f);
  jacobian_targets.atm.back().overlap     = true;
  jacobian_targets.atm.back().overlap_key = component;
  jacobian_targets.atm.back().target_pos  = N + 1;
  switch (component) {
    using enum AtmKey;
    case mag_u:
    case mag_v: jacobian_targets.atm.back().type = mag_w; break;
    case mag_w: jacobian_targets.atm.back().type = mag_v; break;
    default:    throw std::out_of_range("Invalid state");
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

void jacobian_targetsAddOverlappingWindField(
    JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  const auto f = stdr::find_if(jacobian_targets.atm, Jacobian::is_wind);

  ARTS_USER_ERROR_IF(f == jacobian_targets.atm.end(),
                     "Cannot find original wind field for component");

  const AtmKey component = std::get<AtmKey>(f->type);
  const Size N           = jacobian_targets.target_count();

  jacobian_targets.atm.emplace_back(*f);
  jacobian_targets.atm.back().overlap     = true;
  jacobian_targets.atm.back().overlap_key = component;
  jacobian_targets.atm.back().target_pos  = N;
  switch (component) {
    using enum AtmKey;
    case wind_u: jacobian_targets.atm.back().type = wind_v; break;
    case wind_v:
    case wind_w: jacobian_targets.atm.back().type = wind_u; break;
    default:     throw std::out_of_range("Invalid state");
  }

  jacobian_targets.atm.emplace_back(*f);
  jacobian_targets.atm.back().overlap     = true;
  jacobian_targets.atm.back().overlap_key = component;
  jacobian_targets.atm.back().target_pos  = N + 1;
  switch (component) {
    using enum AtmKey;
    case wind_u:
    case wind_v: jacobian_targets.atm.back().type = wind_w; break;
    case wind_w: jacobian_targets.atm.back().type = wind_v; break;
    default:     throw std::out_of_range("Invalid state");
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

void jacobian_targetsAddSensorFrequencyPolyOffset(
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

  make_polyoffset(jacobian_targets.emplace_back(
                      SensorKey{.type             = SensorKeyType::f,
                                .measurement_elem = measurement_elem},
                      d),
                  static_cast<Size>(polyorder),
                  measurement_sensor);
}
