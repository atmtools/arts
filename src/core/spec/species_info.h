#pragma once

#include <enums.h>
#include <format_tags.h>
#include <xml.h>

#include <boost/container_hash/hash.hpp>
#include <string>

struct SpeciesIsotopologueInfo {
  SpeciesEnum  species;
  String       code;
  double       mass;
  double       default_ratio;
  std::int64_t degeneracy;

  constexpr auto operator<=>(const SpeciesIsotopologueInfo& other) const {
    if (auto cmp = species <=> other.species; cmp != 0) return cmp;
    return code <=> other.code;
  }

  constexpr bool operator==(const SpeciesIsotopologueInfo& other) const {
    return species == other.species && code == other.code;
  }

  [[nodiscard]] String name() const;
};

namespace std {
//! Allow SpeciesTag to be used in hashes
template <> struct hash<SpeciesIsotopologueInfo> {
  static std::size_t operator()(const SpeciesIsotopologueInfo& g) {
    std::size_t seed = 0;

    boost::hash_combine(seed, std::hash<SpeciesEnum>{}(g.species));
    boost::hash_combine(seed, std::hash<String>{}(g.code));

    return seed;
  }
};
}  // namespace std

template <> struct xml_io_stream_name<SpeciesIsotopologueInfo> {
  static constexpr std::string_view name = "SpeciesIsotopologueInfo"sv;
};

template <> struct xml_io_stream<SpeciesIsotopologueInfo> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<SpeciesIsotopologueInfo>;

  static void write(std::ostream&                  os,
                    const SpeciesIsotopologueInfo& v,
                    bofstream*            = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, SpeciesIsotopologueInfo& v, bifstream* = nullptr);
};

template <> struct std::formatter<SpeciesIsotopologueInfo> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SpeciesIsotopologueInfo& v, FmtContext& ctx) const {
    return tags.format(ctx, v.species, "-"sv, v.code);
  }
};
