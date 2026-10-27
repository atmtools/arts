#include <jac_polyfit.h>
#include <workspace.h>

#include "properties.h"

namespace {
//! Adds the components not already in jac_targets as targets overlapping the one that is
void add_overlapping_components(JacobianTargets&            jac_targets,
                                const std::array<AtmKey, 3> components,
                                std::string_view            field) {
  const auto orig = stdr::find_if(jac_targets.atm, [&components](const Jacobian::AtmTarget& t) {
    return stdr::any_of(components, [&t](AtmKey k) { return t.type == k; });
  });

  if (orig == jac_targets.atm.end()) arts::user_error("Cannot find original {} for component", field);

  if (Jacobian::has_transformation(*orig))
    arts::user_error("Overlapping targets cannot be transformed (e.g., to relative or logarithmic), but {} is",
                     orig->type);

  const AtmKey  original = std::get<AtmKey>(orig->type);
  const Numeric d        = orig->d;

  jac_targets.finalized = false;

  for (const AtmKey component : components) {
    if (component == original) continue;

    jac_targets.atm.push_back(Jacobian::AtmTarget{
        .type = component, .d = d, .target_pos = jac_targets.target_count(), .overlap = true, .overlap_key = original});
  }
}
}  // namespace

void jac_targetsOff(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  jac_targets.clear();
}

void jac_targetsInit(JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  jac_targetsOff(jac_targets);
}

void jac_targetsFinalize(JacobianTargets&          jac_targets,
                         const AtmField&           atm_field,
                         const SurfaceField&       surf_field,
                         const SubsurfaceField&    subsurf_field,
                         const AbsorptionBands&    abs_bands,
                         const ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  jac_targets.finalize(atm_field, surf_field, subsurf_field, abs_bands, measurement_sensor);
}

void jac_targetsAddSubsurface(JacobianTargets&                                                jac_targets,
                              const Generic<const SubsurfaceKey, const SubsurfacePropertyTag> key,
                              const Numeric&                                                  d) {
  ARTS_TIME_REPORT

  std::visit([&](const auto& selected) { jac_targets.emplace_back(SubsurfaceKeyVal{*selected}, d); }, key);
}

void jac_targetsAddSurface(JacobianTargets&                                          jac_targets,
                           const Generic<const SurfaceKey, const SurfacePropertyTag> key,
                           const Numeric&                                            d) {
  ARTS_TIME_REPORT

  std::visit([&](const auto& selected) { jac_targets.emplace_back(SurfaceKeyVal{*selected}, d); }, key);
}

void jac_targetsAddAtmosphere(JacobianTargets&                    jac_targets,
                              const Generic<const AtmKey,
                                            const QuantumLevelIdentifier,
                                            const ScatteringSpeciesProperty,
                                            const SpeciesEnum,
                                            const SpeciesIsotope> key,
                              const Numeric&                      d) {
  ARTS_TIME_REPORT

  std::visit([&](const auto& selected) { jac_targets.emplace_back(AtmKeyVal{*selected}, d); }, key);
}

void jac_targetsAddAtmosphere(JacobianTargets& jac_targets, const ScatteringSpeciesProperty& key, const Numeric& d) {
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

void jac_targetsAddMagneticField(JacobianTargets& jac_targets, const String& component, const Numeric& d) {
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

  add_overlapping_components(jac_targets, {AtmKey::mag_u, AtmKey::mag_v, AtmKey::mag_w}, "magnetic field");
}

void jac_targetsAddWindField(JacobianTargets& jac_targets, const String& component, const Numeric& d) {
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

  add_overlapping_components(jac_targets, {AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w}, "wind field");
}

void jac_targetsAddSpeciesVMR(JacobianTargets& jac_targets, const SpeciesEnum& species, const Numeric& d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{species}, d);
}

void jac_targetsAddSpeciesIsotopologueRatio(JacobianTargets&      jac_targets,
                                            const SpeciesIsotope& species,
                                            const Numeric&        d) {
  ARTS_TIME_REPORT

  jac_targets.emplace_back(AtmKeyVal{species}, d);
}

void jac_targetsAddLineParameter(JacobianTargets&          jac_targets,
                                 const QuantumIdentifier&  band,
                                 const Index&              line,
                                 const LineByLineVariable& parameter,
                                 const Numeric&            d) {
  ARTS_TIME_REPORT

  if (line < 0) arts::user_error("Line index must be non-negative, got {}", line);
  if (parameter == LineByLineVariable::unused) arts::user_error("A line parameter must be selected");

  jac_targets.emplace_back(LblLineKey{.band = band, .line = static_cast<Size>(line), .var = parameter}, d);
}

void jac_targetsAddLineShapeParameter(JacobianTargets&                 jac_targets,
                                      const QuantumIdentifier&         band,
                                      const Index&                     line,
                                      const SpeciesEnum&               species,
                                      const LineShapeModelVariable&    parameter,
                                      const LineShapeModelCoefficient& coefficient,
                                      const Numeric&                   d) {
  ARTS_TIME_REPORT

  if (line < 0) arts::user_error("Line index must be non-negative, got {}", line);
  if (parameter == LineShapeModelVariable::unused) arts::user_error("A line-shape parameter must be selected");
  if (coefficient == LineShapeModelCoefficient::unused) arts::user_error("A line-shape coefficient must be selected");

  jac_targets.emplace_back(
      LblLineKey{
          .band = band, .line = static_cast<Size>(line), .spec = species, .ls_var = parameter, .ls_coeff = coefficient},
      d);
}

void jac_targetsAddSensorFrequencyPolyOffset(JacobianTargets&          jac_targets,
                                             const ArrayOfSensorObsel& measurement_sensor,
                                             const Numeric&            d,
                                             const Index&              sensor_elem,
                                             const Index&              polyorder) {
  ARTS_TIME_REPORT

  if (polyorder < 0) arts::user_error("Polyorder must be non-negative: {}", polyorder);

  std::vector<std::pair<Index, const void*>> sensor_grid_ptrs;
  for (Size i = 0; i < measurement_sensor.size(); i++) {
    const void* ptr = reinterpret_cast<const void*>(measurement_sensor[i].f_grid_ptr().get());
    if (stdr::none_of(sensor_grid_ptrs | stdv::values, Cmp::eq(ptr))) { sensor_grid_ptrs.emplace_back(i, ptr); }
  }

  if (sensor_grid_ptrs.size() <= static_cast<Size>(sensor_elem))
    arts::user_error(
        "There are only {0} independent sensor frequency grids.  "
        "Cannot select index: {1}, please choose an index less than {0}.",
        sensor_grid_ptrs.size(),
        sensor_elem);

  const Index measurement_elem = sensor_grid_ptrs[sensor_elem].first;

  make_polyoffset(
      jac_targets.emplace_back(SensorKey{.type = SensorKeyType::freq, .measurement_elem = measurement_elem}, d),
      static_cast<Size>(polyorder),
      measurement_sensor);
}
