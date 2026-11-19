#pragma once

#include <enumsInterpolationExtrapolation.h>
#include <format_tags.h>
#include <matpack.h>
#include <operators.h>
#include <xml_io_stream.h>

namespace surface_scattering {
using FunctionalData = NumericTernaryOperator;
using FieldData      = std::variant<SpectralGeodeticField3, Numeric, FunctionalData>;

//! Hold all surface data
struct Data {
  FieldData                  data{FunctionalData{}};
  InterpolationExtrapolation frq_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation frq_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_low{InterpolationExtrapolation::None};

  // Standard
  Data();
  Data(const Data &);
  Data(Data &&) noexcept;
  Data &operator=(const Data &);
  Data &operator=(Data &&) noexcept;

  void adjust_interpolation_extrapolation();

  // Allow copy implicitly from all types
  explicit Data(Numeric);
  Data &operator=(Numeric);

  explicit Data(SpectralGeodeticField3);
  Data &operator=(SpectralGeodeticField3);

  explicit Data(FunctionalData);
  Data &operator=(FunctionalData);

  [[nodiscard]] String data_type() const;

  template <typename T> [[nodiscard]] T get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    if (out == nullptr) arts::user_error("Does not contain correct type");
    return *out;
  }

  template <typename T> [[nodiscard]] T get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    if (out == nullptr) arts::user_error("Does not contain correct type");
    return *out;
  }

  [[nodiscard]] ConstVectorView flat_view() const;

  [[nodiscard]] VectorView flat_view();

  //! Flat weights for the positions on the surface
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> flat_weights(const Numeric &frq,
                                                                    const Numeric &lat,
                                                                    const Numeric &lon) const;

  [[nodiscard]] Numeric at(const Numeric frq, const Numeric lat, const Numeric lon) const;

  //! The data at all frequencies of a position on the surface
  [[nodiscard]] Vector at(const Vector &frq, const Numeric lat, const Numeric lon) const;

  [[nodiscard]] bool ok() const;
};
}  //namespace surface_scattering

using SurfaceSpectralData = surface_scattering::Data;

template <> struct std::formatter<SurfaceSpectralData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SurfaceSpectralData &v, FmtContext &ctx) const {
    if (tags.short_str) return tags.format(ctx, v.data);

    const std::string_view sep = tags.sep();
    if (tags.help) tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx, v.data);

    if (tags.help) {
      tags.format(ctx, sep);
      tags.add_if_bracket(ctx, "["sv);
      tags.format(ctx, v.frq_upp, sep, v.frq_low, sep, v.lat_upp, sep, v.lat_low, sep, v.lon_upp, sep, v.lon_low);
      tags.add_if_bracket(ctx, "]]"sv);
    }

    return ctx.out();
  }
};

template <> struct xml_io_stream<SurfaceSpectralData> {
  static constexpr std::string_view type_name = "SurfaceSpectralData"sv;

  static void write(std::ostream              &os,
                    const SurfaceSpectralData &x,
                    bofstream                 *pbofs = nullptr,
                    std::string_view           name  = ""sv);

  static void read(std::istream &is, SurfaceSpectralData &x, bifstream *pbifs = nullptr);
};