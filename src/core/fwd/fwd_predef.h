#pragma once

#include <predefined_absorption_models.h>
#include <species_tags.h>

#include <memory>

namespace fwd::predef {
struct full {
  std::vector<Species::IsotopeRecord> tags;
  Numeric P;
  Numeric T;
  Absorption::PredefinedModel::VMRS vmrs;
  std::shared_ptr<PredefinedModelData> predefined_model_data;

  full() = default;

  full(const AtmPoint& atm_point,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::shared_ptr<PredefinedModelData>& data);

  [[nodiscard]] Complex operator()(Numeric f) const;
  void operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector operator()(const Vector& fs) const;
};
}  // namespace fwd::predef
