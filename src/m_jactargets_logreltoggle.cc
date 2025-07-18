#include <workspace.h>

#include "debug.h"
#include "jac_logrel.h"

namespace {
void jacobian_targetsToggleLogRelAtmTargetImpl(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.atm()) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logrelfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleLogRelSurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.surf()) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logrelfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleLogRelSubsurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.subsurf()) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_logrelfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}
}  // namespace

// Atm

void jacobian_targetsToggleLogRelAtmTarget(JacobianTargets& jacobian_targets,
                                           const AtmField& f,
                                           const AtmKey& key) {
  jacobian_targetsToggleLogRelAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogRelAtmTarget(JacobianTargets& jacobian_targets,
                                           const AtmField& f,
                                           const SpeciesEnum& key) {
  jacobian_targetsToggleLogRelAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogRelAtmTarget(JacobianTargets& jacobian_targets,
                                           const AtmField& f,
                                           const SpeciesIsotope& key) {
  jacobian_targetsToggleLogRelAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogRelAtmTarget(JacobianTargets& jacobian_targets,
                                           const AtmField& f,
                                           const QuantumLevelIdentifier& key) {
  jacobian_targetsToggleLogRelAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogRelAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const ScatteringSpeciesProperty& key) {
  jacobian_targetsToggleLogRelAtmTargetImpl(jacobian_targets, f, key);
}

// Surface

void jacobian_targetsToggleLogRelSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKey& key) {
  jacobian_targetsToggleLogRelSurfaceTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleLogRelSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfacePropertyTag& key) {
  jacobian_targetsToggleLogRelSurfaceTargetImpl(jacobian_targets, f, key);
}

// Subsurface

void jacobian_targetsToggleLogRelSubsurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKey& key) {
  jacobian_targetsToggleLogRelSubsurfaceTargetImpl(jacobian_targets, f, key);
}
