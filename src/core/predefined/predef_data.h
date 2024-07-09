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

  void resize(const std::vector<std::size_t> &);
  [[nodiscard]] std::vector<std::size_t> sizes() const;

  friend std::ostream &operator<<(std::ostream &, const WaterData &);
  friend std::istream &operator>>(std::istream &, WaterData &);
};
}  // namespace MT_CKD400

struct ModelName {
  static std::vector<std::size_t> sizes() { return {}; };
  static constexpr void resize(const std::vector<std::size_t> &) {}

  friend std::ostream &operator<<(std::ostream &, const ModelName &);
  friend std::istream &operator>>(std::istream &, ModelName &);
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

  friend std::ostream &operator<<(std::ostream &, const Model &);
};

std::string_view model_name(const ModelVariant &);
ModelVariant model_data(const std::string_view);
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;

template <>
struct std::formatter<Absorption::PredefinedModel::MT_CKD400::WaterData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(
      const Absorption::PredefinedModel::MT_CKD400::WaterData &v,
      FmtContext &ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.ref_temp,
                sep,
                v.ref_press,
                sep,
                v.ref_h2o_vmr,
                sep,
                v.self_absco_ref,
                sep,
                v.for_absco_ref,
                sep,
                v.wavenumbers,
                sep,
                v.self_texp);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<Absorption::PredefinedModel::ModelName> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Absorption::PredefinedModel::ModelName &,
                              FmtContext &ctx) const {
    return std::format_to(ctx.out(), "[]");
  }
};

template <>
struct std::formatter<PredefinedModelData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const PredefinedModelData &v,
                              FmtContext &ctx) const {
    return tags.format(ctx, v.data);
  }
};

#endif  // predefined_predef_data_h
