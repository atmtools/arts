
#pragma once

#include "subsurf_field.h"
#include "subsurf_spectral_field.h"

template <> struct xml_io_stream<SubsurfaceField> {
  static constexpr std::string_view type_name = "SubsurfaceField"sv;

  static void write(std::ostream&          os,
                    const SubsurfaceField& x,
                    bofstream*             pbofs = nullptr,
                    std::string_view       name  = ""sv);

  static void read(std::istream& is, SubsurfaceField& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfaceData> {
  static constexpr std::string_view type_name = "SubsurfaceData"sv;

  static void write(std::ostream&         os,
                    const SubsurfaceData& x,
                    bofstream*            pbofs = nullptr,
                    std::string_view      name  = ""sv);

  static void read(std::istream& is, SubsurfaceData& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfacePoint> {
  static constexpr std::string_view type_name = "SubsurfacePoint"sv;

  static void write(std::ostream&          os,
                    const SubsurfacePoint& x,
                    bofstream*             pbofs = nullptr,
                    std::string_view       name  = ""sv);

  static void read(std::istream& is, SubsurfacePoint& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfacePropertyTag> {
  static constexpr std::string_view type_name = "SubsurfacePropertyTag"sv;

  static void write(std::ostream&                os,
                    const SubsurfacePropertyTag& x,
                    bofstream*                   pbofs = nullptr,
                    std::string_view             name  = ""sv);

  static void read(std::istream& is, SubsurfacePropertyTag& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfaceSpectralField> {
  static constexpr std::string_view type_name = "SubsurfaceSpectralField"sv;

  static void write(std::ostream&                  os,
                    const SubsurfaceSpectralField& x,
                    bofstream*                     pbofs = nullptr,
                    std::string_view               name  = ""sv);

  static void read(std::istream& is, SubsurfaceSpectralField& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfaceSpectralData> {
  static constexpr std::string_view type_name = "SubsurfaceSpectralData"sv;

  static void write(std::ostream&                 os,
                    const SubsurfaceSpectralData& x,
                    bofstream*                    pbofs = nullptr,
                    std::string_view              name  = ""sv);

  static void read(std::istream& is, SubsurfaceSpectralData& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfaceSpectralPropertyTag> {
  static constexpr std::string_view type_name = "SubsurfaceSpectralPropertyTag"sv;

  static void write(std::ostream&                        os,
                    const SubsurfaceSpectralPropertyTag& x,
                    bofstream*                           pbofs = nullptr,
                    std::string_view                     name  = ""sv);

  static void read(std::istream& is, SubsurfaceSpectralPropertyTag& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SubsurfaceSpectralPoint> {
  static constexpr std::string_view type_name = "SubsurfaceSpectralPoint"sv;

  static void write(std::ostream&                  os,
                    const SubsurfaceSpectralPoint& x,
                    bofstream*                     pbofs = nullptr,
                    std::string_view               name  = ""sv);

  static void read(std::istream& is, SubsurfaceSpectralPoint& x, bifstream* pbifs = nullptr);
};
