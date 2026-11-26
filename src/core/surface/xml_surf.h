#pragma once

#include <xml.h>

#include "surf_field.h"
#include "surf_spectral_field.h"

template <> struct xml_io_stream<SurfaceField> {
  static constexpr std::string_view type_name = "SurfaceField"sv;

  static void write(std::ostream& os, const SurfaceField& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, SurfaceField& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfaceData> {
  static constexpr std::string_view type_name = "SurfaceData"sv;

  static void write(std::ostream& os, const SurfaceData& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, SurfaceData& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfacePoint> {
  static constexpr std::string_view type_name = "SurfacePoint"sv;

  static void write(std::ostream& os, const SurfacePoint& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, SurfacePoint& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfacePropertyTag> {
  static constexpr std::string_view type_name = "SurfacePropertyTag"sv;

  static void write(std::ostream&             os,
                    const SurfacePropertyTag& x,
                    bofstream*                pbofs = nullptr,
                    std::string_view          name  = ""sv);

  static void read(std::istream& is, SurfacePropertyTag& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfaceSpectralField> {
  static constexpr std::string_view type_name = "SurfaceSpectralField"sv;

  static void write(std::ostream&               os,
                    const SurfaceSpectralField& x,
                    bofstream*                  pbofs = nullptr,
                    std::string_view            name  = ""sv);

  static void read(std::istream& is, SurfaceSpectralField& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfaceSpectralData> {
  static constexpr std::string_view type_name = "SurfaceSpectralData"sv;

  static void write(std::ostream&              os,
                    const SurfaceSpectralData& x,
                    bofstream*                 pbofs = nullptr,
                    std::string_view           name  = ""sv);

  static void read(std::istream& is, SurfaceSpectralData& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream<SurfaceSpectralPropertyTag> {
  static constexpr std::string_view type_name = "SurfaceSpectralPropertyTag"sv;

  static void write(std::ostream&                     os,
                    const SurfaceSpectralPropertyTag& x,
                    bofstream*                        pbofs = nullptr,
                    std::string_view                  name  = ""sv);

  static void read(std::istream& is, SurfaceSpectralPropertyTag& x, bifstream* pbifs = nullptr);
};
