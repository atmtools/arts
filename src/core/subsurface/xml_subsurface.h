
#pragma once

#include "subsurface_field.h"

template <>
struct xml_io_stream<SubsurfaceField> {
  static constexpr std::string_view type_name = "SubsurfaceField"sv;

  static void write(std::ostream& os,
                    const SubsurfaceField& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SubsurfaceField& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<SubsurfaceData> {
  static constexpr std::string_view type_name = "SubsurfaceData"sv;

  static void write(std::ostream& os,
                    const SubsurfaceData& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SubsurfaceData& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<SubsurfacePoint> {
  static constexpr std::string_view type_name = "SubsurfacePoint"sv;

  static void write(std::ostream& os,
                    const SubsurfacePoint& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SubsurfacePoint& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<SubsurfacePropertyTag> {
  static constexpr std::string_view type_name = "SubsurfacePropertyTag"sv;

  static void write(std::ostream& os,
                    const SubsurfacePropertyTag& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SubsurfacePropertyTag& x,
                   bifstream* pbifs = nullptr);
};
