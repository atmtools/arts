#pragma once

#include <isotopologues.h>

struct HitranSpeciesInfo {
  SpeciesIsotope spec;
  Index          hitind;
  char           hitchar;
  Numeric        ratio;
};

template <> struct std::formatter<HitranSpeciesInfo> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const HitranSpeciesInfo& v, FmtContext& ctx) const {
    return tags.format(ctx, v.spec.FullName(), " "sv, v.hitind, " "sv, std::string_view(&v.hitchar, 1), " "sv, v.ratio);
  }
};

template <> struct xml_io_stream<HitranSpeciesInfo> {
  static constexpr std::string_view type_name = "HitranSpeciesInfo"sv;

  static void write(std::ostream&            os,
                    const HitranSpeciesInfo& x,
                    bofstream*               pbofs = nullptr,
                    std::string_view         name  = ""sv);

  static void read(std::istream& is, HitranSpeciesInfo& x, bifstream* pbifs = nullptr);
};
