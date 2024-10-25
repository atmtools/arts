#include <jacobian.h>

void model_state_vectorSize(Vector& model_state_vector,
                            const JacobianTargets& jacobian_targets) {
  model_state_vector.resize(jacobian_targets.x_size());
}

void model_state_vectorZero(Vector& model_state_vector) {
  model_state_vector = 0.0;
}

void atmospheric_fieldFromModelState(AtmField& atmospheric_field,
                                     const Vector& model_state_vector,
                                     const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.atm()) {
    target.update(atmospheric_field, model_state_vector);
  }
}

void model_state_vectorFromAtmosphere(Vector& model_state_vector,
                                      const AtmField& atmospheric_field,
                                      const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.atm()) {
    target.update(model_state_vector, atmospheric_field);
  }
}

void surface_fieldFromModelState(SurfaceField& surface_field,
                                 const Vector& model_state_vector,
                                 const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.surf()) {
    target.update(surface_field, model_state_vector);
  }
}

void model_state_vectorFromSurface(Vector& model_state_vector,
                                   const SurfaceField& surface_field,
                                   const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.surf()) {
    target.update(model_state_vector, surface_field);
  }
}

void absorption_bandsFromModelState(AbsorptionBands& absorption_bands,
                                    const Vector& model_state_vector,
                                    const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.line()) {
    target.update(absorption_bands, model_state_vector);
  }
}

void model_state_vectorFromBands(Vector& model_state_vector,
                                 const AbsorptionBands& absorption_bands,
                                 const JacobianTargets& jacobian_targets) {
  for (auto& target : jacobian_targets.line()) {
    target.update(model_state_vector, absorption_bands);
  }
}
