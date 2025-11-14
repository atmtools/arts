#include <debug.h>
#include <jac_rh.h>
#include <operators.h>
#include <workspace.h>

namespace {
void jac_targetsToggleRelativeHumidityAtmTargetImpl(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const AtmKeyVal& key,
    const NumericUnaryOperator& psat,
    const Index& nonnegative) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.atm) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<rhinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_rhfit(t, f, psat, nonnegative != 0);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}
}  // namespace

// Atm

void jac_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const AtmKey& key,
    const Index& nonnegative) {
  jac_targetsToggleRelativeHumidityAtmTargetImpl(
      jac_targets, f, key, psat, nonnegative);
}

void jac_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const SpeciesEnum& key,
    const Index& nonnegative) {
  jac_targetsToggleRelativeHumidityAtmTargetImpl(
      jac_targets, f, key, psat, nonnegative);
}

void jac_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const SpeciesIsotope& key,
    const Index& nonnegative) {
  jac_targetsToggleRelativeHumidityAtmTargetImpl(
      jac_targets, f, key, psat, nonnegative);
}

void jac_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const QuantumLevelIdentifier& key,
    const Index& nonnegative) {
  jac_targetsToggleRelativeHumidityAtmTargetImpl(
      jac_targets, f, key, psat, nonnegative);
}

void jac_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jac_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const ScatteringSpeciesProperty& key,
    const Index& nonnegative) {
  jac_targetsToggleRelativeHumidityAtmTargetImpl(
      jac_targets, f, key, psat, nonnegative);
}
