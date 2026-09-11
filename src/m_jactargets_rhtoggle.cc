#include <debug.h>
#include <jac_rh.h>
#include <operators.h>
#include <workspace.h>

namespace {
void jac_targetsToggleRelativeHumidityAtmTargetImpl(JacobianTargets&            jac_targets,
                                                    const AtmField&             f,
                                                    const AtmKeyVal&            key,
                                                    const NumericUnaryOperator& psat,
                                                    const Index&                nonnegative) {
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
        jac_targetsToggleRelativeHumidityAtmTargetImpl(jac_targets, f, *selected, psat, nonnegative);
      },
      key);
}
