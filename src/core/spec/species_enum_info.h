#pragma once

#include <format_tags.h>
#include <xml.h>

#include <string>

struct SpeciesEnumInfo {
  Size   enum_value;
  String shortname;
  String longname;

  constexpr auto operator<=>(const SpeciesEnumInfo& other) const { return enum_value <=> other.enum_value; }

  constexpr bool operator==(const SpeciesEnumInfo& other) const { return enum_value == other.enum_value; }
};

template <> struct std::hash<SpeciesEnumInfo> {
  static std::size_t operator()(const SpeciesEnumInfo& g) { return g.enum_value; }
};

template <> struct xml_io_stream_name<SpeciesEnumInfo> {
  static constexpr std::string_view name = "SpeciesEnumInfo"sv;
};

template <> struct xml_io_stream<SpeciesEnumInfo> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<SpeciesEnumInfo>;

  static void write(std::ostream& os, const SpeciesEnumInfo& v, bofstream* = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, SpeciesEnumInfo& v, bifstream* = nullptr);
};

template <> struct std::formatter<SpeciesEnumInfo> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SpeciesEnumInfo& v, FmtContext& ctx) const {
    return tags.format(ctx, v.shortname, tags.sep(), v.longname);
  }
};
