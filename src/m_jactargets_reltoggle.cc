#include <debug.h>
#include <jac_rel.h>
#include <workspace.h>

namespace {
void jacobian_targetsToggleRelativeAtmTargetImpl(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.atm) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<relinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_relfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleRelativeSurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.surf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<relinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_relfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}

void jacobian_targetsToggleRelativeSubsurfaceTargetImpl(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jacobian_targets.subsurf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<relinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        make_relfit(t, f);
      }
      return;
    }
  }
  ARTS_USER_ERROR("Could not find target {}", key)
}
}  // namespace

// Atm

void jacobian_targetsToggleRelativeAtmTarget(JacobianTargets& jacobian_targets,
                                             const AtmField& f,
                                             const AtmKey& key) {
  jacobian_targetsToggleRelativeAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleRelativeAtmTarget(JacobianTargets& jacobian_targets,
                                             const AtmField& f,
                                             const SpeciesEnum& key) {
  jacobian_targetsToggleRelativeAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleRelativeAtmTarget(JacobianTargets& jacobian_targets,
                                             const AtmField& f,
                                             const SpeciesIsotope& key) {
  jacobian_targetsToggleRelativeAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleRelativeAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const QuantumLevelIdentifier& key) {
  jacobian_targetsToggleRelativeAtmTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleRelativeAtmTarget(
    JacobianTargets& jacobian_targets,
    const AtmField& f,
    const ScatteringSpeciesProperty& key) {
  jacobian_targetsToggleRelativeAtmTargetImpl(jacobian_targets, f, key);
}

// Surface

void jacobian_targetsToggleRelativeSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfaceKey& key) {
  jacobian_targetsToggleRelativeSurfaceTargetImpl(jacobian_targets, f, key);
}

void jacobian_targetsToggleRelativeSurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SurfaceField& f,
    const SurfacePropertyTag& key) {
  jacobian_targetsToggleRelativeSurfaceTargetImpl(jacobian_targets, f, key);
}

// Subsurface

void jacobian_targetsToggleRelativeSubsurfaceTarget(
    JacobianTargets& jacobian_targets,
    const SubsurfaceField& f,
    const SubsurfaceKey& key) {
  jacobian_targetsToggleRelativeSubsurfaceTargetImpl(jacobian_targets, f, key);
}
