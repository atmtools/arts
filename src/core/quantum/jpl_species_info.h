#pragma once

#include <enums.h>
#include <isotopologues.h>

struct JplSpeciesInfo {
  Index id;
  SpeciesIsotope spec;
  bool has_qn{};

  Numeric QT0;
  Numeric T0;
};

template <>
struct std::formatter<JplSpeciesInfo> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const JplSpeciesInfo& v, FmtContext& ctx) const {
    return tags.format(ctx,
                       v.spec,
                       " "sv,
                       v.id,
                       " "sv,
                       v.has_qn,
                       " "sv,
                       v.QT0,
                       " "sv,
                       v.T0);
  }
};

template <>
struct xml_io_stream<JplSpeciesInfo> {
  static constexpr std::string_view type_name = "JplSpeciesInfo"sv;

  static void write(std::ostream& os,
                    const JplSpeciesInfo& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   JplSpeciesInfo& x,
                   bifstream* pbifs = nullptr);
};
