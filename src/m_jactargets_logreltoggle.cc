#include <debug.h>
#include <jac_logrel.h>
#include <workspace.h>

namespace {
void jac_targetsToggleLogRelAtmTargetImpl(JacobianTargets& jac_targets, const AtmField& f, const AtmKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.atm) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        ARTS_USER_ERROR_IF(Jacobian::overlaps(jac_targets.atm, t),
                           "Overlapping targets cannot be transformed (e.g., to relative or logarithmic), but {} is",
                           key)
        make_logrelfit(t, f);
      }
      return;
    }
  }
  arts::user_error("Could not find target {}", key);
}

void jac_targetsToggleLogRelSurfaceTargetImpl(JacobianTargets&     jac_targets,
                                              const SurfaceField&  f,
                                              const SurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.surf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        ARTS_USER_ERROR_IF(Jacobian::overlaps(jac_targets.surf, t),
                           "Overlapping targets cannot be transformed (e.g., to relative or logarithmic), but {} is",
                           key)
        make_logrelfit(t, f);
      }
      return;
    }
  }
  arts::user_error("Could not find target {}", key);
}

void jac_targetsToggleLogRelSubsurfaceTargetImpl(JacobianTargets&        jac_targets,
                                                 const SubsurfaceField&  f,
                                                 const SubsurfaceKeyVal& key) {
  ARTS_TIME_REPORT

  for (auto& t : jac_targets.subsurf) {
    if (t.type == key) {
      if (t.inverse_jacobian.target<logrelinv>() != nullptr) {
        t.inverse_jacobian = {};
        t.inverse_state    = {};
        t.transform_state  = {};
      } else {
        ARTS_USER_ERROR_IF(Jacobian::overlaps(jac_targets.subsurf, t),
                           "Overlapping targets cannot be transformed (e.g., to relative or logarithmic), but {} is",
                           key)
        make_logrelfit(t, f);
      }
      return;
    }
  }
  arts::user_error("Could not find target {}", key);
}
}  // namespace

// Atm

void jac_targetsToggleLogRelAtmTarget(JacobianTargets&                    jac_targets,
                                      const AtmField&                     f,
                                      const Generic<const AtmKey,
                                                    const QuantumLevelIdentifier,
                                                    const ScatteringSpeciesProperty,
                                                    const SpeciesEnum,
                                                    const SpeciesIsotope> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleLogRelAtmTargetImpl(jac_targets, f, *selected); }, key);
}

// Surface

void jac_targetsToggleLogRelSurfaceTarget(JacobianTargets&                                          jac_targets,
                                          const SurfaceField&                                       f,
                                          const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit([&](const auto& selected) { jac_targetsToggleLogRelSurfaceTargetImpl(jac_targets, f, *selected); }, key);
}

// Subsurface

void jac_targetsToggleLogRelSubsurfaceTarget(JacobianTargets&       jac_targets,
                                             const SubsurfaceField& f,
                                             const SubsurfaceKey&   key) {
  jac_targetsToggleLogRelSubsurfaceTargetImpl(jac_targets, f, key);
}
