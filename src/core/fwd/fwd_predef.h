#pragma once

#include <predefined_absorption_models.h>
#include <species_tags.h>

#include <memory>

namespace fwd::predef {
class full {
  Absorption::PredefinedModel::VMRS vmrs;
  std::shared_ptr<AtmPoint> atm;
  std::shared_ptr<PredefinedModelData> data;

  void adapt();

public:
  full() = default;

  full(std::shared_ptr<AtmPoint> atm,
      std::shared_ptr<PredefinedModelData> data);

  [[nodiscard]] Complex operator()(Numeric f) const;
  void operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector operator()(const Vector& fs) const;

  void set_model(std::shared_ptr<PredefinedModelData> data);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};
}  // namespace fwd::predef
