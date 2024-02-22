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

class Model {
  using ModelData = std::unordered_map<SpeciesIsotopeRecord, ModelVariant>;

  ModelData data{};

 public:
  Model() = default;
  Model(const Model& d) = default;
  Model(Model&& d) = default;
  Model& operator=(const Model& d) = default;
  Model& operator=(Model&& d) = default;

  void insert(const Model& other);

  template <ModelVariantConvertible T>
  void set(const SpeciesIsotopeRecord& tag, const T& t) {
    ARTS_USER_ERROR_IF(not is_predefined_model(tag),
                       "The tag must be of type PredefinedModel")
    data[tag] = static_cast<ModelVariant>(t);
  }

  template <ModelVariantConvertible T, Index tag>
  void set(const T& t)
    requires(tag < Species::Isotopologues.size())
  {
    set(Species::Isotopologues[tag], t);
  }

  template <ModelVariantType T>
  [[nodiscard]] const T& get(const SpeciesIsotopeRecord& tag) const try {
    ARTS_USER_ERROR_IF(not is_predefined_model(tag),
                       "The tag must be of type PredefinedModel")
    return std::get<T>(data.at(tag));
  } catch (const std::bad_variant_access&) {
    throw std::runtime_error(
        var_string("The tag ", tag, " does not have the requested type"));
  } catch (std::out_of_range&) {
    throw std::runtime_error(
        var_string("The tag ", tag, " does not exist in the model"));
  }

  template <ModelVariantType T, Index tag>
  [[nodiscard]] const T& get() const
    requires(tag < Species::Isotopologues.size())
  {
    return get<T>(Species::Isotopologues[tag]);
  }

  [[nodiscard]] const ModelVariant& at(const SpeciesIsotopeRecord& tag) const;

  [[nodiscard]] ModelVariant& at(const SpeciesIsotopeRecord& tag);

  void clear();
  [[nodiscard]] auto size() const { return data.size(); }
  [[nodiscard]] auto empty() const { return data.empty(); }
  [[nodiscard]] auto begin() const { return data.begin(); }
  [[nodiscard]] auto end() const { return data.end(); }
  friend std::ostream& operator<<(std::ostream&, const Model&);
};

std::string_view model_name(const ModelVariant&);
ModelVariant model_data(const std::string_view);
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;

#endif  // predefined_predef_data_h
