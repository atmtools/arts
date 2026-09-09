#include <debug.h>
#include <jac_rel.h>
#include <workspace.h>

namespace {
void jac_targetsToggleRelativeAtmTargetImpl(JacobianTargets& jac_targets, const AtmField& f, const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.atm) {
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

void jac_targetsToggleRelativeSurfaceTargetImpl(JacobianTargets&     jac_targets,
                                                const SurfaceField&  f,
                                                const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.surf) {
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

void jac_targetsToggleRelativeSubsurfaceTargetImpl(JacobianTargets&        jac_targets,
                                                   const SubsurfaceField&  f,
                                                   const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.subsurf) {
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

void jac_targetsToggleRelativeAtmTarget(JacobianTargets&                               jac_targets,
                                        const AtmField&                                f,
                                        const Generic<const AtmKey,
                                                      const SpeciesEnum,
                                                      const SpeciesIsotope,
                                                      const QuantumLevelIdentifier,
                                                      const ScatteringSpeciesProperty> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleRelativeAtmTargetImpl(jac_targets, f, *selected); }, key);
}

// Surface

void jac_targetsToggleRelativeSurfaceTarget(JacobianTargets&                                          jac_targets,
                                            const SurfaceField&                                       f,
                                            const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleRelativeSurfaceTargetImpl(jac_targets, f, *selected); }, key);
}

// Subsurface

void jac_targetsToggleRelativeSubsurfaceTarget(JacobianTargets&       jac_targets,
                                               const SubsurfaceField& f,
                                               const SubsurfaceKey&   key) {
  jac_targetsToggleRelativeSubsurfaceTargetImpl(jac_targets, f, key);
}
