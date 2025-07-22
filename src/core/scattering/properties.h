#ifndef ARTS_CORE_SCATTERING_PROPERTIES_H_
#define ARTS_CORE_SCATTERING_PROPERTIES_H_

#include <array.h>
#include <enums.h>
#include <xml.h>

#include <algorithm>
#include <string_view>

#include "mystring.h"

/*** ScatteringSpeciesProperty
 *
 * Used to uniquely identify an atmospheric field that holds properties
 * of a given scattering species.
 *
 * */
struct ScatteringSpeciesProperty {
  std::string species_name;
  ParticulateProperty pproperty;

  constexpr auto operator<=>(const ScatteringSpeciesProperty& other) const =
      default;

  friend std::ostream& operator<<(std::ostream& os,
                                  const ScatteringSpeciesProperty& ssp) {
  return os << ssp.species_name << "_" << ssp.pproperty;
}
};

namespace std {
template <>
struct hash<ScatteringSpeciesProperty> {
  std::size_t operator()(const ScatteringSpeciesProperty& ssp) const {
    return std::hash<String>{}(ssp.species_name +
                               std::string(toString(ssp.pproperty)));
  }
};
}  // namespace std

template <>
struct std::formatter<ScatteringSpeciesProperty> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ScatteringSpeciesProperty& v,
                              FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(
        ctx.out(), "{}{}_{}{}", quote, v.species_name, v.pproperty, quote);
  }
};

template <>
struct xml_io_stream<ScatteringSpeciesProperty> {
  static constexpr std::string_view type_name = "ScatteringSpeciesProperty"sv;

  static void write(std::ostream& os,
                    const ScatteringSpeciesProperty& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   ScatteringSpeciesProperty& x,
                   bifstream* pbifs = nullptr);
};

#endif
