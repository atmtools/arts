#pragma once

#include <xml.h>

#include "functional_atm.h"
#include "functional_atm_field.h"
#include "functional_gravity.h"

using NumericTernary = std::function<Numeric(Numeric, Numeric, Numeric)>;

namespace ternary {
//! Method that returns the x-value for retrieval from a NumericTernary
ConstVectorView get_xc(const NumericTernary& f);

//! Method that returns the x-value for retrieval from a NumericTernary
VectorView get_xm(NumericTernary& f);

//! Method that returns the scaling-value for retrieval from a NumericTernary
std::vector<std::pair<Index, Numeric>> get_w(const NumericTernary& f,
                                             Numeric alt,
                                             Numeric lat,
                                             Numeric lon);

//! Method checks that the NumericTernary is a field function (i.e., requires multiple derivatives)
bool is_field_function(const NumericTernary& f);
}  // namespace ternary

template <>
struct xml_io_stream_name<NumericTernary> {
  constexpr static std::string_view name = "NumericTernary"sv;
};

template <>
struct xml_io_stream_functional<NumericTernary> {
  using func_t    = Numeric (*)(Numeric, Numeric, Numeric);
  using structs_t = std::variant<Atm::IGRF13,
                                 Atm::SchmidthLegendre,
                                 Atm::HydrostaticPressure,
                                 Atm::MagnitudeField,
                                 Atm::External,
                                 EllipsoidGravity>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <>
struct std::formatter<NumericTernary> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const NumericTernary&, FmtContext& ctx) const {
    return tags.format(ctx, "<numeric-ternary>");
  }
};
