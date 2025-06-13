#pragma once

#include <predefined_absorption_models.h>
#include <species_tags.h>

#include <memory>

namespace fwd::predef {
class full {
  std::shared_ptr<AtmPoint> atm;
  std::shared_ptr<PredefinedModelData> data;

  void adapt();

 public:
  full();
  full(const full&);
  full(full&&) noexcept;
  full& operator=(const full&);
  full& operator=(full&&) noexcept;

  full(std::shared_ptr<AtmPoint> atm,
       std::shared_ptr<PredefinedModelData> data);

  [[nodiscard]] Complex operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<PredefinedModelData> data);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};
}  // namespace fwd::predef
