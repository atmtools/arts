#include <workspace.h>

void model_state_vectorInit(Vector& model_state_vector,
                            const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  model_state_vector.resize(jacobian_targets.x_size());
  model_state_vector = 0.0;
}

void atmospheric_fieldFromModelState(AtmField& atmospheric_field,
                                     const Vector& model_state_vector,
                                     const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.atm()) {
    target.update(atmospheric_field, model_state_vector);
  }
}

void model_state_vectorFromAtmosphere(Vector& model_state_vector,
                                      const AtmField& atmospheric_field,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.atm()) {
    target.update(model_state_vector, atmospheric_field);
  }
}

void surface_fieldFromModelState(SurfaceField& surface_field,
                                 const Vector& model_state_vector,
                                 const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.surf()) {
    target.update(surface_field, model_state_vector);
  }
}

void subsurface_fieldFromModelState(SubsurfaceField& subsurface_field,
                                    const Vector& model_state_vector,
                                    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.subsurf()) {
    target.update(subsurface_field, model_state_vector);
  }
}

void model_state_vectorFromSurface(Vector& model_state_vector,
                                   const SurfaceField& surface_field,
                                   const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.surf()) {
    target.update(model_state_vector, surface_field);
  }
}

void model_state_vectorFromSubsurface(Vector& model_state_vector,
                                      const SubsurfaceField& subsurface_field,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.subsurf()) {
    target.update(model_state_vector, subsurface_field);
  }
}

void absorption_bandsFromModelState(AbsorptionBands& absorption_bands,
                                    const Vector& model_state_vector,
                                    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.line()) {
    target.update(absorption_bands, model_state_vector);
  }
}

void model_state_vectorFromBands(Vector& model_state_vector,
                                 const AbsorptionBands& absorption_bands,
                                 const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.line()) {
    target.update(model_state_vector, absorption_bands);
  }
}

void measurement_sensorFromModelState(ArrayOfSensorObsel& measurement_sensor,
                                      const Vector& model_state_vector,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.sensor()) {
    target.update(measurement_sensor, model_state_vector);
  }
}

void model_state_vectorFromSensor(Vector& model_state_vector,
                                  const ArrayOfSensorObsel& measurement_sensor,
                                  const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.sensor()) {
    target.update(model_state_vector, measurement_sensor);
  }
}
