#include <workspace.h>

void model_state_vectorInit(Vector& model_state_vector,
                            const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  model_state_vector.resize(jacobian_targets.x_size());
  model_state_vector = 0.0;
}

void model_state_vectorPerturbations(Vector& model_state_vector,
                                     const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  model_state_vectorInit(model_state_vector, jacobian_targets);

  for (auto& target : jacobian_targets.atm) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }

  for (auto& target : jacobian_targets.surf) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }

  for (auto& target : jacobian_targets.subsurf) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }

  for (auto& target : jacobian_targets.line) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }

  for (auto& target : jacobian_targets.sensor) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }

  for (auto& target : jacobian_targets.error) {
    const Range r(target.x_start, target.x_size);
    model_state_vector[r] = target.d;
  }
}

////// Update the fields from the model state vector

void atm_fieldFromModelState(AtmField& atm_field,
                             const Vector& model_state_vector,
                             const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.atm) {
    target.update_model(atm_field, model_state_vector);
  }
}

void surf_fieldFromModelState(SurfaceField& surf_field,
                              const Vector& model_state_vector,
                              const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.surf) {
    target.update_model(surf_field, model_state_vector);
  }
}

void subsurf_fieldFromModelState(SubsurfaceField& subsurf_field,
                                 const Vector& model_state_vector,
                                 const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.subsurf) {
    target.update_model(subsurf_field, model_state_vector);
  }
}

void abs_bandsFromModelState(AbsorptionBands& abs_bands,
                             const Vector& model_state_vector,
                             const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.line) {
    target.update_model(abs_bands, model_state_vector);
  }
}

void measurement_sensorFromModelState(ArrayOfSensorObsel& measurement_sensor,
                                      const Vector& model_state_vector,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.sensor) {
    target.update_model(measurement_sensor, model_state_vector);
  }
}

////// Update the model state vector from the fields

void model_state_vectorFromAtmosphere(Vector& model_state_vector,
                                      const AtmField& atm_field,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.atm) {
    target.update_state(model_state_vector, atm_field);
  }
}

void model_state_vectorFromSurface(Vector& model_state_vector,
                                   const SurfaceField& surf_field,
                                   const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.surf) {
    target.update_state(model_state_vector, surf_field);
  }
}

void model_state_vectorFromSubsurface(Vector& model_state_vector,
                                      const SubsurfaceField& subsurf_field,
                                      const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.subsurf) {
    target.update_state(model_state_vector, subsurf_field);
  }
}

void model_state_vectorFromBands(Vector& model_state_vector,
                                 const AbsorptionBands& abs_bands,
                                 const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.line) {
    target.update_state(model_state_vector, abs_bands);
  }
}

void model_state_vectorFromSensor(Vector& model_state_vector,
                                  const ArrayOfSensorObsel& measurement_sensor,
                                  const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.sensor) {
    target.update_state(model_state_vector, measurement_sensor);
  }
}

///// Update the measurement Jacobian from the model state vector

void measurement_jacobianAtmosphereTransformation(
    Matrix& measurement_jacobian,
    const Vector& model_state_vector,
    const AtmField& field,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.atm) {
    target.update_jac(measurement_jacobian, model_state_vector, field);
  }
}

void measurement_jacobianSurfaceTransformation(
    Matrix& measurement_jacobian,
    const Vector& model_state_vector,
    const SurfaceField& field,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.surf) {
    target.update_jac(measurement_jacobian, model_state_vector, field);
  }
}

void measurement_jacobianSubsurfaceTransformation(
    Matrix& measurement_jacobian,
    const Vector& model_state_vector,
    const SubsurfaceField& field,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.subsurf) {
    target.update_jac(measurement_jacobian, model_state_vector, field);
  }
}

void measurement_jacobianBandTransformation(
    Matrix& measurement_jacobian,
    const Vector& model_state_vector,
    const AbsorptionBands& field,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.line) {
    target.update_jac(measurement_jacobian, model_state_vector, field);
  }
}

void measurement_jacobianSensorTransformation(
    Matrix& measurement_jacobian,
    const Vector& model_state_vector,
    const ArrayOfSensorObsel& field,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  for (auto& target : jacobian_targets.sensor) {
    target.update_jac(measurement_jacobian, model_state_vector, field);
  }
}
