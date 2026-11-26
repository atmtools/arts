#pragma once

#include <format_tags.h>
#include <xml.h>

#include <variant>

#include "general_tro_spectral.h"
#include "scattering_medium.h"

namespace scattering {
/** The callbacks of a CallbackScatteringSpecies, one type per representation of the bulk properties
 *
 * Each takes the points as a MediumPoint and a SpectralMediumPoint and deals
 * with the media it supports.
 * A new representation is a new alternative here, with its method.
 */
using CallbackScatterer = std::variant<ScatteringGeneralSpectralTRO>;

/** A scattering species given by a user's callback
 *
 * A method passes the points of the medium it is called in to the callback
 * held, whose own method checks the shape of what the user's function gives.
 * A call for a representation the callback does not give throws.
 */
struct CallbackScatteringSpecies {
  CallbackScatterer callback{};

  CallbackScatteringSpecies() = default;

  //! Any callback is a species
  template <typename T> requires(std::constructible_from<CallbackScatterer, T> and
                                 not std::same_as<std::remove_cvref_t<T>, CallbackScatteringSpecies>)
  CallbackScatteringSpecies(T&& x) : callback(std::forward<T>(x)) {}

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(const Point&         point,
                                                         const SpectralPoint& spectral_point,
                                                         const Vector&        f_grid,
                                                         Index                degree) const {
    return std::visit(
        [&]<typename T>(const T& c) -> ScatteringTroSpectralVector {
          if constexpr (requires { c.tro_spectral(MediumPoint{point}, SpectralMediumPoint{spectral_point}, f_grid, degree); }) {
            return c.tro_spectral(MediumPoint{point}, SpectralMediumPoint{spectral_point}, f_grid, degree);
          } else {
            not_given<T>("its TRO Legendre series");
          }
        },
        callback);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(
      const Point&         point,
      const SpectralPoint& spectral_point,
      const Vector&        f_grid,
      const Vector&        za_inc_grid,
      const Vector&        za_scat_grid,
      Index                max_mode) const {
    return std::visit(
        [&]<typename T>(const T& c) -> BulkScatteringProperties<Format::ARO, Representation::Fourier> {
          if constexpr (requires {
                          c.aro_fourier(MediumPoint{point},
                                        SpectralMediumPoint{spectral_point},
                                        f_grid,
                                        za_inc_grid,
                                        za_scat_grid,
                                        max_mode);
                        }) {
            return c.aro_fourier(
                MediumPoint{point}, SpectralMediumPoint{spectral_point}, f_grid, za_inc_grid, za_scat_grid, max_mode);
          } else {
            not_given<T>("its azimuthal Fourier modes");
          }
        },
        callback);
  }

 private:
  template <typename T> [[noreturn]] static void not_given(std::string_view what) {
    arts::user_error("The {} callback does not give {}", xml_io_stream<T>::type_name, what);
  }
};
}  // namespace scattering

using CallbackScatteringSpecies = scattering::CallbackScatteringSpecies;

template <> struct std::formatter<CallbackScatteringSpecies> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const CallbackScatteringSpecies& x, FmtContext& ctx) const {
    return tags.format(ctx, x.callback);
  }
};

template <> struct xml_io_stream<CallbackScatteringSpecies> {
  static constexpr std::string_view type_name = "CallbackScatteringSpecies";

  static void write(std::ostream&, const CallbackScatteringSpecies&, bofstream* = nullptr, std::string_view = ""sv);

  static void read(std::istream&, CallbackScatteringSpecies&, bifstream* = nullptr);
};
