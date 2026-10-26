#include <debug.h>
#include <jac_log.h>
#include <jac_logrel.h>
#include <jac_rel.h>
#include <jac_rh.h>
#include <operators.h>
#include <workspace.h>

namespace {
//! Removes the transformation with inverse Inv from the target, or replaces any transformation by calling make
template <typename Inv, typename T>
void toggle(std::vector<T>& targets, const decltype(T::type)& key, const auto& make) {
  ARTS_TIME_REPORT

  const auto t = stdr::find(targets, key, &T::type);
  ARTS_USER_ERROR_IF(t == targets.end(), "Could not find target {}", key)

  if (t->inverse_jacobian.template target<Inv>() != nullptr) {
    t->inverse_jacobian = {};
    t->inverse_state    = {};
    t->transform_state  = {};
  } else {
    ARTS_USER_ERROR_IF(Jacobian::overlaps(targets, *t),
                       "Overlapping targets cannot be transformed (e.g., to relative or logarithmic), but {} is",
                       key)
    make(*t);
  }
}
}  // namespace

// Relative

void jac_targetsToggleRelativeAtmTarget(JacobianTargets&                    jac_targets,
                                        const AtmField&                     f,
                                        const Generic<const AtmKey,
                                                      const QuantumLevelIdentifier,
                                                      const ScatteringSpeciesProperty,
                                                      const SpeciesEnum,
                                                      const SpeciesIsotope> key) {
  std::visit(
      [&](const auto& selected) { toggle<relinv>(jac_targets.atm, *selected, [&f](auto& t) { make_relfit(t, f); }); },
      key);
}

void jac_targetsToggleRelativeSurfaceTarget(JacobianTargets&                                          jac_targets,
                                            const SurfaceField&                                       f,
                                            const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit(
      [&](const auto& selected) { toggle<relinv>(jac_targets.surf, *selected, [&f](auto& t) { make_relfit(t, f); }); },
      key);
}

void jac_targetsToggleRelativeSubsurfaceTarget(JacobianTargets&       jac_targets,
                                               const SubsurfaceField& f,
                                               const SubsurfaceKey&   key) {
  toggle<relinv>(jac_targets.subsurf, key, [&f](auto& t) { make_relfit(t, f); });
}

// Logarithmic

void jac_targetsToggleLogarithmicAtmTarget(JacobianTargets&                    jac_targets,
                                           const AtmField&                     f,
                                           const Generic<const AtmKey,
                                                         const QuantumLevelIdentifier,
                                                         const ScatteringSpeciesProperty,
                                                         const SpeciesEnum,
                                                         const SpeciesIsotope> key) {
  std::visit(
      [&](const auto& selected) { toggle<loginv>(jac_targets.atm, *selected, [&f](auto& t) { make_logfit(t, f); }); },
      key);
}

void jac_targetsToggleLogarithmicSurfaceTarget(JacobianTargets&                                          jac_targets,
                                               const SurfaceField&                                       f,
                                               const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit(
      [&](const auto& selected) { toggle<loginv>(jac_targets.surf, *selected, [&f](auto& t) { make_logfit(t, f); }); },
      key);
}

void jac_targetsToggleLogarithmicSubsurfaceTarget(JacobianTargets&       jac_targets,
                                                  const SubsurfaceField& f,
                                                  const SubsurfaceKey&   key) {
  toggle<loginv>(jac_targets.subsurf, key, [&f](auto& t) { make_logfit(t, f); });
}

// Logarithmic relative

void jac_targetsToggleLogRelAtmTarget(JacobianTargets&                    jac_targets,
                                      const AtmField&                     f,
                                      const Generic<const AtmKey,
                                                    const QuantumLevelIdentifier,
                                                    const ScatteringSpeciesProperty,
                                                    const SpeciesEnum,
                                                    const SpeciesIsotope> key) {
  std::visit(
      [&](const auto& selected) {
        toggle<logrelinv>(jac_targets.atm, *selected, [&f](auto& t) { make_logrelfit(t, f); });
      },
      key);
}

void jac_targetsToggleLogRelSurfaceTarget(JacobianTargets&                                          jac_targets,
                                          const SurfaceField&                                       f,
                                          const Generic<const SurfaceKey, const SurfacePropertyTag> key) {
  std::visit(
      [&](const auto& selected) {
        toggle<logrelinv>(jac_targets.surf, *selected, [&f](auto& t) { make_logrelfit(t, f); });
      },
      key);
}

void jac_targetsToggleLogRelSubsurfaceTarget(JacobianTargets&       jac_targets,
                                             const SubsurfaceField& f,
                                             const SubsurfaceKey&   key) {
  toggle<logrelinv>(jac_targets.subsurf, key, [&f](auto& t) { make_logrelfit(t, f); });
}

// Relative humidity

void jac_targetsToggleRelativeHumidityAtmTarget(JacobianTargets&                    jac_targets,
                                                const AtmField&                     f,
                                                const NumericUnaryOperator&         psat,
                                                const Generic<const AtmKey,
                                                              const QuantumLevelIdentifier,
                                                              const ScatteringSpeciesProperty,
                                                              const SpeciesEnum,
                                                              const SpeciesIsotope> key,
                                                const Index&                        nonnegative) {
  std::visit(
      [&](const auto& selected) {
        toggle<rhinv>(jac_targets.atm, *selected, [&](auto& t) { make_rhfit(t, f, psat, nonnegative != 0); });
      },
      key);
}
