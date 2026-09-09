#include <debug.h>
#include <jac_log.h>
#include <workspace.h>

namespace {
void jac_targetsToggleLogarithmicAtmTargetImpl(JacobianTargets& jac_targets, const AtmField& f, const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.atm) {
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

void jac_targetsToggleLogarithmicSurfaceTargetImpl(JacobianTargets&     jac_targets,
                                                   const SurfaceField&  f,
                                                   const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.surf) {
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

void jac_targetsToggleLogarithmicSubsurfaceTargetImpl(JacobianTargets&        jac_targets,
                                                      const SubsurfaceField&  f,
                                                      const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.subsurf) {
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

void jac_targetsToggleLogarithmicAtmTarget(JacobianTargets&                               jac_targets,
                                           const AtmField&                                f,
                                           const Generic<const AtmKey,
                                                         const SpeciesEnum,
                                                         const SpeciesIsotope,
                                                         const QuantumLevelIdentifier,
                                                         const ScatteringSpeciesProperty> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleLogarithmicAtmTargetImpl(jac_targets, f, *selected); }, key);
}

// Surface

void jac_targetsToggleLogarithmicSurfaceTarget(JacobianTargets&                                          jac_targets,
                                               const SurfaceField&                                       f,
                                               const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleLogarithmicSurfaceTargetImpl(jac_targets, f, *selected); },
             key);
}

// Subsurface

void jac_targetsToggleLogarithmicSubsurfaceTarget(JacobianTargets&       jac_targets,
                                                  const SubsurfaceField& f,
                                                  const SubsurfaceKey&   key) {
  jac_targetsToggleLogarithmicSubsurfaceTargetImpl(jac_targets, f, key);
}
