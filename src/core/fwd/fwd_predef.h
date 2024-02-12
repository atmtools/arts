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
  full(const full&) = default;
  full(full&&) = default;
  full& operator=(const full&) = default;
  full& operator=(full&&) = default;

  full(std::shared_ptr<AtmPoint> atm,
       std::shared_ptr<PredefinedModelData> data);

  [[nodiscard]] Complex operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<PredefinedModelData> data);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};
}  // namespace fwd::predef
