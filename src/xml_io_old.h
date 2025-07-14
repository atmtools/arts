#pragma once

#include <lbl_data.h>
#include <xml.h>

struct ArtscatMeta {
  lbl::line data;
  QuantumIdentifier quantumidentity{};
  bool bad{true};
};

using ArrayOfArtscatMeta = std::vector<ArtscatMeta>;

template <>
struct xml_io_stream<ArrayOfArtscatMeta> {
  static constexpr std::string_view type_name = "ARTSCAT"sv;
  static void write(std::ostream&,
                    const ArrayOfArtscatMeta&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);
  static void read(std::istream& is_xml,
                   ArrayOfArtscatMeta& meta,
                   bifstream* = nullptr);
};