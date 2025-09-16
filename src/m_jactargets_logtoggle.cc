#include <debug.h>
#include <jac_log.h>
#include <workspace.h>

namespace {
void jacobian_targetsToggleLogarithmicAtmTargetImpl(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.atm) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<loginv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleLogarithmicSurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.surf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<loginv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleLogarithmicSubsurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.subsurf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<loginv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}
}  // namespace

// Atm

void jacobian_targetsToggleLogarithmicAtmTarget(
    JacobianTargets& jacobian_targets, const AtmField& f, const AtmKey& key) {
  jacobian_targetsToggleLogarithmicAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogarithmicAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const SpeciesEnum& key) {
  jacobian_targetsToggleLogarithmicAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogarithmicAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const SpeciesIsotope& key) {
  jacobian_targetsToggleLogarithmicAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogarithmicAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const QuantumLevelIdentifier& key) {
  jacobian_targetsToggleLogarithmicAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogarithmicAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const ScatteringSpeciesProperty& key) {
  jacobian_targetsToggleLogarithmicAtmTargetImpl(jacobian_targets, f, key);
}

// Surface

void jacobian_targetsToggleLogarithmicSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKey& key) {
  jacobian_targetsToggleLogarithmicSurfaceTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogarithmicSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfacePropertyTag& key) {
  jacobian_targetsToggleLogarithmicSurfaceTargetImpl(jacobian_targets, f, key);
}

// Subsurface

void jacobian_targetsToggleLogarithmicSubsurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKey& key) {
  jacobian_targetsToggleLogarithmicSubsurfaceTargetImpl(
      jacobian_targets, f, key);
}
