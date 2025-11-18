#include "fwd_predef.h"

#include <atm.h>
#include <jacobian.h>
#include <rtepack.h>

namespace fwd::predef {
full::full()                           = default;
full::full(const full&)                = default;
full::full(full&&) noexcept            = default;
full& full::operator=(const full&)     = default;
full& full::operator=(full&&) noexcept = default;

void full::adapt() try {
  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  if (not data) {
    return;
  }

  if (data->empty()) {
    return;
  }
}
ARTS_METHOD_ERROR_CATCH

full::full(std::shared_ptr<AtmPoint> atm_,
           std::shared_ptr<PredefinedModelData> data_)
    : atm(std::move(atm_)), data(std::move(data_)) {
  adapt();
}

Complex full::operator()(const Numeric frequency) const {
  if (not data) {
    return {};
  }

  PropmatVector propmat_clearsky(1);
  PropmatMatrix dpropmat_clearsky_dx;
  JacobianTargets jac_targets;
  Vector f_grid{frequency};

  for (auto& [tag, mod] : *data) {
    Absorption::PredefinedModel::compute(propmat_clearsky,
                                         dpropmat_clearsky_dx,
                                         tag,
                                         f_grid,
                                         *atm,
                                         jac_targets,
                                         mod);
  }

  return propmat_clearsky[0].A();
}

void full::set_model(std::shared_ptr<PredefinedModelData> data_) {
  data = std::move(data_);
  adapt();
}

void full::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}
}  // namespace fwd::predef
