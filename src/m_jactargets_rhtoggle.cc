#include <jac_rh.h>
#include <workspace.h>

#include "debug.h"
#include "operators.h"

namespace {
void jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const AtmKeyVal& key,
    const NumericUnaryOperator& psat,
    const Index& nonnegative) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.atm()) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<pairinv<AtmField>>() != nullptr) {
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

void jacobian_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const AtmKey& key,
    const Index& nonnegative) {
  jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
      jacobian_targets, f, key, psat, nonnegative);
}

void jacobian_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const SpeciesEnum& key,
    const Index& nonnegative) {
  jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
      jacobian_targets, f, key, psat, nonnegative);
}

void jacobian_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const SpeciesIsotope& key,
    const Index& nonnegative) {
  jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
      jacobian_targets, f, key, psat, nonnegative);
}

void jacobian_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const QuantumLevelIdentifier& key,
    const Index& nonnegative) {
  jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
      jacobian_targets, f, key, psat, nonnegative);
}

void jacobian_targetsToggleRelativeHumidityAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const NumericUnaryOperator& psat,
    const ScatteringSpeciesProperty& key,
    const Index& nonnegative) {
  jacobian_targetsToggleRelativeHumidityAtmTargetImpl(
      jacobian_targets, f, key, psat, nonnegative);
}
