#pragma once

#include <xml.h>

#include "atm_field.h"

template <>
struct xml_io_stream<AtmPoint> {
  static constexpr std::string_view type_name = "AtmPoint"sv;

  static void write(std::ostream& os,
                    const AtmPoint& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, AtmPoint& x, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<AtmField> {
  static constexpr std::string_view type_name = "AtmField"sv;

  static void write(std::ostream& os,
                    const AtmField& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, AtmField& x, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<AtmData> {
  static constexpr std::string_view type_name = "AtmData"sv;

  static void write(std::ostream& os,
                    const AtmData& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, AtmData& x, bifstream* pbifs = nullptr);
};
