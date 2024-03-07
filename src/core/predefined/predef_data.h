#ifndef predefined_predef_data_h
#define predefined_predef_data_h

#include <debug.h>
#include <enums.h>
#include <isotopologues.h>
#include <matpack.h>

#include <exception>
#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "matpack_data.h"

namespace Absorption::PredefinedModel {
namespace MT_CKD400 {
struct WaterData {
  Numeric ref_press;
  Numeric ref_temp;
  Numeric ref_h2o_vmr;
  Vector self_absco_ref;
  Vector for_absco_ref;
  Vector wavenumbers;
  Vector self_texp;

  void resize(const std::vector<std::size_t>&);
  [[nodiscard]] std::vector<std::size_t> sizes() const;

  friend std::ostream& operator<<(std::ostream&, const WaterData&);
  friend std::istream& operator>>(std::istream&, WaterData&);
};
}  // namespace MT_CKD400

struct ModelName {
  static std::vector<std::size_t> sizes() { return {}; };
  static constexpr void resize(const std::vector<std::size_t>&) {}

  friend std::ostream& operator<<(std::ostream&, const ModelName&);
  friend std::istream& operator>>(std::istream&, ModelName&);
};

using ModelVariant = std::variant<ModelName, MT_CKD400::WaterData>;

template <typename T>
concept ModelVariantType = requires(T t) {
  { std::get<T>(static_cast<ModelVariant>(t)) };
};

template <typename T>
concept ModelVariantConvertible = ModelVariantType<T> or requires(T t) {
  { static_cast<ModelVariant>(t) };
};

struct Model {
  using ModelData = std::unordered_map<SpeciesIsotope, ModelVariant>;

  ModelData data{};

  friend std::ostream& operator<<(std::ostream&, const Model&);
};

std::string_view model_name(const ModelVariant&);
ModelVariant model_data(const std::string_view);
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;

#endif  // predefined_predef_data_h
