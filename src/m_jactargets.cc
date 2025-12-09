#include <jac_polyfit.h>
#include <workspace.h>

void jac_targetsOff(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  jac_targets.clear();
}

void jac_targetsConditionalClear(JacobianTargets& jac_targets,
                                 const Index& do_jac) {
  ARTS_TIME_REPORT

  if (do_jac == 0) jac_targetsOff(jac_targets);
}

void jac_targetsInit(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  jac_targetsOff(jac_targets);
}

void jac_targetsFinalize(JacobianTargets& jac_targets,
                         const AtmField& atm_field,
                         const SurfaceField& surf_field,
                         const SubsurfaceField& subsurf_field,
                         const AbsorptionBands& abs_bands,
                         const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  jac_targets.finalize(
      atm_field, surf_field, subsurf_field, abs_bands, measurement_sensor);
}

void jac_targetsAddSubsurface(JacobianTargets& jac_targets,
                              const SubsurfaceKey& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(SubsurfaceKeyVal{key}, d);
}

void jac_targetsAddSubsurface(JacobianTargets& jac_targets,
                              const SubsurfacePropertyTag& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(SubsurfaceKeyVal{key}, d);
}

void jac_targetsAddSurface(JacobianTargets& jac_targets,
                           const SurfaceKey& key,
                           const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jac_targetsAddSurface(JacobianTargets& jac_targets,
                           const SurfacePropertyTag& key,
                           const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(SurfaceKeyVal{key}, d);
}

void jac_targetsAddAtmosphere(JacobianTargets& jac_targets,
                              const AtmKey& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{key}, d);
}

void jac_targetsAddAtmosphere(JacobianTargets& jac_targets,
                              const SpeciesEnum& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{key}, d);
}

void jac_targetsAddAtmosphere(JacobianTargets& jac_targets,
                              const SpeciesIsotope& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{key}, d);
}

void jac_targetsAddAtmosphere(JacobianTargets& jac_targets,
                              const QuantumLevelIdentifier& key,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{key}, d);
}

void jac_targetsAddTemperature(JacobianTargets& jac_targets, const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{AtmKey::t}, d);
}

void jac_targetsAddPressure(JacobianTargets& jac_targets, const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{AtmKey::p}, d);
}

void jac_targetsAddMagneticField(JacobianTargets& jac_targets,
                                 const String& component,
                                 const Numeric& d) {
  ARTS_TIME_REPORT

  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u: jac_targets.emplace_back(AtmKeyVal{AtmKey::mag_u}, d); break;
    case v: jac_targets.emplace_back(AtmKeyVal{AtmKey::mag_v}, d); break;
    case w: jac_targets.emplace_back(AtmKeyVal{AtmKey::mag_w}, d); break;
  }
}

void jac_targetsAddOverlappingMagneticField(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  const auto fptr = stdr::find_if(jac_targets.atm, Jacobian::is_mag);

  ARTS_USER_ERROR_IF(fptr == jac_targets.atm.end(),
                     "Cannot find original magnetic field for component");

  const AtmKey component = std::get<AtmKey>(fptr->type);
  const Size N           = jac_targets.target_count();

  // Must copy because emplace_back destroys the pointer
  auto copy        = *fptr;
  copy.overlap     = true;
  copy.overlap_key = component;
  copy.target_pos  = N;

  jac_targets.atm.emplace_back(copy);
  switch (component) {
    using enum AtmKey;
    case mag_u: jac_targets.atm.back().type = mag_v; break;
    case mag_v:
    case mag_w: jac_targets.atm.back().type = mag_u; break;
    default:    throw std::out_of_range("Invalid state");
  }

  copy.target_pos = N + 1;
  jac_targets.atm.emplace_back(std::move(copy));
  switch (component) {
    using enum AtmKey;
    case mag_u:
    case mag_v: jac_targets.atm.back().type = mag_w; break;
    case mag_w: jac_targets.atm.back().type = mag_v; break;
    default:    throw std::out_of_range("Invalid state");
  }
}

void jac_targetsAddWindField(JacobianTargets& jac_targets,
                             const String& component,
                             const Numeric& d) {
  ARTS_TIME_REPORT

  using enum FieldComponent;
  switch (to<FieldComponent>(component)) {
    case u: jac_targets.emplace_back(AtmKeyVal{AtmKey::wind_u}, d); break;
    case v: jac_targets.emplace_back(AtmKeyVal{AtmKey::wind_v}, d); break;
    case w: jac_targets.emplace_back(AtmKeyVal{AtmKey::wind_w}, d); break;
  }
}

void jac_targetsAddOverlappingWindField(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  const auto fptr = stdr::find_if(jac_targets.atm, Jacobian::is_wind);

  ARTS_USER_ERROR_IF(fptr == jac_targets.atm.end(),
                     "Cannot find original wind field for component");

  const AtmKey component = std::get<AtmKey>(fptr->type);
  const Size N           = jac_targets.target_count();

  // Must copy because emplace_back destroys the pointer
  auto copy        = *fptr;
  copy.overlap     = true;
  copy.overlap_key = component;
  copy.target_pos  = N;

  jac_targets.atm.emplace_back(std::move(copy));
  switch (component) {
    using enum AtmKey;
    case wind_u: jac_targets.atm.back().type = wind_v; break;
    case wind_v:
    case wind_w: jac_targets.atm.back().type = wind_u; break;
    default:     throw std::out_of_range("Invalid state");
  }

  copy.target_pos = N + 1;
  jac_targets.atm.emplace_back(std::move(copy));
  switch (component) {
    using enum AtmKey;
    case wind_u:
    case wind_v: jac_targets.atm.back().type = wind_w; break;
    case wind_w: jac_targets.atm.back().type = wind_v; break;
    default:     throw std::out_of_range("Invalid state");
  }
}

void jac_targetsAddSpeciesVMR(JacobianTargets& jac_targets,
                              const SpeciesEnum& species,
                              const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{species}, d);
}

void jac_targetsAddSpeciesIsotopologueRatio(JacobianTargets& jac_targets,
                                            const SpeciesIsotope& species,
                                            const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{species}, d);
}

void jac_targetsAddSensorFrequencyPolyOffset(
    JacobianTargets& jac_targets,
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
    if (stdr::none_of(sensor_grid_ptrs | stdv::values, Cmp::eq(ptr))) {
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

  make_polyoffset(
      jac_targets.emplace_back(SensorKey{.type = SensorKeyType::freq,
                                         .measurement_elem = measurement_elem},
                               d),
      static_cast<Size>(polyorder),
      measurement_sensor);
}
