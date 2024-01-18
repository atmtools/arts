#include "fwd_predef.h"

#include <algorithm>

#include "atm.h"
#include "jacobian.h"
#include "rtepack.h"

namespace fwd::predef {
full::full(const AtmPoint& atm_point,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const std::shared_ptr<PredefinedModelData>& data)
    : P(atm_point.pressure),
      T(atm_point.temperature),
      vmrs(Absorption::PredefinedModel::VMRS(atm_point)),
      predefined_model_data(data) {
  for (auto& specs : allspecs) {
    for (auto& spec : specs) {
      if (spec.type == Species::TagType::Predefined) {
        tags.push_back(spec.Isotopologue());
      }
    }
  }
}

Complex full::at(Numeric f) const {
  PropmatVector propmat_clearsky(1);
  PropmatMatrix dpropmat_clearsky_dx;
  JacobianTargets jacobian_targets;
  Vector f_grid{f};

  for (auto& tag : tags) {
    Absorption::PredefinedModel::compute(propmat_clearsky,
                                         dpropmat_clearsky_dx,
                                         tag,
                                         f_grid,
                                         P,
                                         T,
                                         vmrs,
                                         jacobian_targets,
                                         *predefined_model_data);
  }

  return propmat_clearsky[0].A();
}

void full::at(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), abs.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector full::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}
}  // namespace fwd::predef
