#pragma once

#include <ostream>
#include <unordered_map>

#include "isotopologues.h"
#include "lbl_temperature_model.h"
#include "species.h"

namespace lbl::linemixing {
struct species_data {
  temperature::data scaling{LineShapeModelType::T0, {0}};
  temperature::data beta{LineShapeModelType::T0, {0}};
  temperature::data lambda{LineShapeModelType::T0, {0}};
  temperature::data collisional_distance{LineShapeModelType::T0, {0}};

  [[nodiscard]] Numeric Q(const Rational J,
                          const Numeric T,
                          const Numeric T0,
                          const Numeric energy) const;

  [[nodiscard]] Numeric Omega(const Numeric T,
                              const Numeric T0,
                              const Numeric mass,
                              const Numeric other_mass,
                              const Numeric energy_x,
                              const Numeric energy_xm2) const;
};  // species_data

using species_data_map = std::unordered_map<SpeciesEnum, species_data>;

//! FIXME: Should behave as an unordered map of unordered maps, but isn't one because we don't understand pybind11
struct isot_map {
  std::unordered_map<SpeciesIsotope, species_data_map> data{};

  species_data_map& operator[](const SpeciesIsotope& key) { return data[key]; }
  [[nodiscard]] auto find(const SpeciesIsotope& key) const {
    return data.find(key);
  }
  [[nodiscard]] auto begin() { return data.begin(); }
  [[nodiscard]] auto end() { return data.end(); }
  [[nodiscard]] auto begin() const { return data.begin(); }
  [[nodiscard]] auto end() const { return data.end(); }
  [[nodiscard]] auto cbegin() const { return data.cbegin(); }
  [[nodiscard]] auto cend() const { return data.cend(); }
  void clear() { data.clear(); }
  void reserve(const size_t n) { data.reserve(n); }
  [[nodiscard]] std::size_t size() const { return data.size(); }
  [[nodiscard]] bool empty() const { return data.empty(); }

  friend std::ostream& operator<<(std::ostream&, const isot_map&);
};  // isot_map
}  // namespace lbl::linemixing

using LinemixingEcsData = lbl::linemixing::isot_map;

template <>
struct std::formatter<lbl::linemixing::species_data> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::linemixing::species_data& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.scaling,
                sep,
                v.beta,
                sep,
                v.lambda,
                sep,
                v.collisional_distance);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<LinemixingEcsData> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const LinemixingEcsData& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.data);
  }
};
