#include "fwd_predef.h"

#include <algorithm>

#include "rtepack.h"

namespace fwd::predef {
full::full(Numeric p,
           Numeric t,
           const Vector& allvmrs,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const std::shared_ptr<PredefinedModelData>& data)
    : P(p),
      T(t),
      vmrs(Absorption::PredefinedModel::VMRS(allspecs, allvmrs)),
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
  ArrayOfRetrievalQuantity jacobian_quantities;
  Vector f_grid{f};

  for (auto& tag : tags) {
    Absorption::PredefinedModel::compute(propmat_clearsky,
                                         dpropmat_clearsky_dx,
                                         tag,
                                         f_grid,
                                         P,
                                         T,
                                         vmrs,
                                         jacobian_quantities,
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
