#pragma once

#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <format>
#include <iosfwd>
#include <string_view>

#include "vdisort.h"

/** The workspace types of VDISORT, the polarized counterparts of those of DISORT (disort.h).
 *
 * They mirror DisortSettings and DisortRadiance member for member, and the
 * fluxes are the SpectralFluxField of disort.h.  Where DISORT holds a scalar,
 * VDISORT holds a Stokes vector, and the phase function's Legendre moments
 * become those of the six elements of the scattering-plane matrix of randomly
 * oriented particles.
 */
using VDisortBDRF         = vdisort::BDRF;
using MatrixOfVDisortBDRF = matpack::data_t<VDisortBDRF, 2>;

/** The radiance field of VDISORT on the grids of DisortRadiance.
 *
 * data is [freq_grid, layer (alt_grid[1:]), azi_grid, zen_grid] Stokes
 * vectors.  zen_grid and azi_grid are those of the propagation direction, as
 * the line of sight of a PropagationPathPoint.
 */
struct VDisortRadiance {
  AscendingGrid  freq_grid;
  DescendingGrid alt_grid;
  AziGrid        azi_grid;
  ZenGrid        zen_grid;
  StokvecTensor4 data;

  /** Replace all grids and resize the radiance tensor to their Cartesian product. */
  void resize(AscendingGrid freq_grid, DescendingGrid alt_grid, AziGrid azi_grid, ZenGrid zen_grid);

  /** Reorder the zenith dimension from the solver's stream order, as DisortRadiance::sort. */
  void sort(const Vector& solver_mu);
};

/** The settings of a VDISORT calculation, member for member those of DisortSettings.
 *
 * The phase matrix is given by the Legendre moments of the six elements
 * [F11, F12, F22, F33, F34, F44] of the scattering-plane matrix of randomly
 * oriented particles: element e is sum_l (2 l + 1) legendre_coefficients[f,
 * layer, l, e] P_l(cos Theta), normalised so that the F11 moment of degree 0
 * is 1 (DISORT's moments are those of F11).  VDISORT's Fourier coefficients
 * on its streams are made from them for each frequency in set().
 */
struct VDisortSettings {
  Index quadrature_dimension{0};
  Index legendre_polynomial_dimension{0};
  Index fourier_mode_dimension{0};

  // Grids
  AscendingGrid  freq_grid{};
  DescendingGrid alt_grid{};  // levels not layers

  // freq_grid.size()
  Vector solar_azimuth_angle{};

  // freq_grid.size()
  Vector solar_zenith_angle{};

  // freq_grid.size(), the Stokes irradiance of the beam on the horizontal plane
  StokvecVector solar_source{};

  // freq_grid.size() x nbrdf
  MatrixOfVDisortBDRF bidirectional_reflectance_distribution_functions{};

  // freq_grid.size() x [alt_grid.size() - 1]
  Matrix optical_thicknesses{};

  // freq_grid.size() x [alt_grid.size() - 1]
  Matrix single_scattering_albedo{};

  // freq_grid.size() x [alt_grid.size() - 1] x nsrc
  StokvecTensor3 source_polynomial{};

  // freq_grid.size() x [alt_grid.size() - 1] x legendre_polynomial_dimension_full x 6
  Tensor4 legendre_coefficients{};

  // freq_grid.size() x 2 (cosine, sine) x fourier_mode_dimension x quadrature_dimension / 2
  StokvecTensor4 upward_boundary_condition{};

  // freq_grid.size() x 2 (cosine, sine) x fourier_mode_dimension x quadrature_dimension / 2
  StokvecTensor4 downward_boundary_condition{};

  /** Construct empty VDISORT settings. */
  VDisortSettings() = default;

  /** Replace grids and dimensions, then resize every frequency- and layer-dependent input. */
  void resize(Index          quadrature_dimension,
              Index          legendre_polynomial_dimension,
              Index          fourier_mode_dimension,
              AscendingGrid  freq_grid,
              DescendingGrid alt_grid_);

  /** Return the number of configured frequencies. */
  [[nodiscard]] Index frequency_count() const { return freq_grid.size(); }
  /** Return the number of atmospheric layers. */
  [[nodiscard]] Index layer_count() const { return alt_grid.size() - 1; }

  /** Allocate a VDISORT solver with the dimensions represented by these settings. */
  [[nodiscard]] vdisort::main_data init() const;
  /** Populate and solve an existing VDISORT solver for frequency index iv. */
  vdisort::main_data& set(vdisort::main_data&, Index iv) const;
  /** Return whether every member agrees with the grids and dimensions. */
  [[nodiscard]] bool ok() const;
  /** Validate settings dimensions and grids. */
  void check() const;
};

/** VDISORT's azimuth [rad, in [0, 2 pi)] of a propagation direction of ARTS azimuth aa [deg].  ARTS counts azimuth
 * clockwise seen from above and VDISORT counterclockwise (vdisort_arts.h), so phi = -aa.  A PropagationPathPoint's
 * line of sight is the propagation direction of its radiation, so the radiance field's azimuth grid, like its zenith
 * grid, is of propagation directions.  The settings' solar azimuth is that of the line of sight toward the sun, so the
 * beam propagates at the solar azimuth + 180 deg. */
[[nodiscard]] Numeric vdisort_azimuth(Numeric aa);

/** The zenith grid [deg] of the streams and their quadrature weights, as DISORT's gridded_weights() */
[[nodiscard]] ZenGriddedField1 vdisort_quadrature_weights(const vdisort::main_data& dis);

template <> struct std::formatter<VDisortBDRF> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const VDisortBDRF&, FmtContext& ctx) const {
    return tags.format(ctx, "VDisortBDRF"sv);
  }
};

template <> struct xml_io_stream<VDisortBDRF> {
  static constexpr std::string_view type_name = "VDisortBDRF"sv;

  static void write(std::ostream& os, const VDisortBDRF& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, VDisortBDRF& x, bifstream* pbifs = nullptr);
};

template <> struct xml_io_stream_name<VDisortRadiance> {
  constexpr static std::string_view name = "VDisortRadiance"sv;
};

template <> struct xml_io_stream_aggregate<VDisortRadiance> {
  constexpr static bool value = true;
};

template <> struct std::formatter<VDisortRadiance> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext> FmtContext::iterator format(const VDisortRadiance& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(ctx,
                         "freq_grid.shape: "sv,
                         v.freq_grid.shape(),
                         "\n"sv,
                         "alt_grid.shape:  "sv,
                         v.alt_grid.shape(),
                         "\n"sv,
                         "azi_grid.shape:  "sv,
                         v.azi_grid.shape(),
                         "\n"sv,
                         "zen_grid.shape:  "sv,
                         v.zen_grid.shape(),
                         "\n"sv,
                         "data.shape:      "sv,
                         v.data.shape());
    }

    const auto sep = tags.sep();
    return tags.format(ctx, v.freq_grid, sep, v.alt_grid, sep, v.azi_grid, sep, v.zen_grid, sep, v.data);
  }
};

template <> struct std::formatter<VDisortSettings> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const VDisortSettings& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(ctx,
                         "quadrature_dimension:          "sv,
                         v.quadrature_dimension,
                         "\nlegendre_polynomial_dimension: "sv,
                         v.legendre_polynomial_dimension,
                         "\nfourier_mode_dimension:        "sv,
                         v.fourier_mode_dimension,
                         "\nfreq_grid.size():              "sv,
                         v.freq_grid.size(),
                         "\nalt_grid.size():               "sv,
                         v.alt_grid.size(),
                         "\n\nsolar_source.shape():                                     "sv,
                         v.solar_source.shape(),
                         " - should be freq_grid.size().\nsolar_zenith_angle.shape():                               "sv,
                         v.solar_zenith_angle.shape(),
                         " - should be freq_grid.size().\nsolar_azimuth_angle.shape():                              "sv,
                         v.solar_azimuth_angle.shape(),
                         " - should be freq_grid.size().\nbidirectional_reflectance_distribution_functions.shape(): "sv,
                         v.bidirectional_reflectance_distribution_functions.shape(),
                         " - should be freq_grid.size() x nbrdf.\noptical_thicknesses.shape():                              "sv,
                         v.optical_thicknesses.shape(),
                         " - should be freq_grid.size() x (alt_grid.size() - 1).\nsingle_scattering_albedo.shape():                         "sv,
                         v.single_scattering_albedo.shape(),
                         " - should be freq_grid.size() x (alt_grid.size() - 1).\nsource_polynomial.shape():                                "sv,
                         v.source_polynomial.shape(),
                         " - should be freq_grid.size() x (alt_grid.size() - 1) x nsrc.\nlegendre_coefficients.shape():                            "sv,
                         v.legendre_coefficients.shape(),
                         " - should be freq_grid.size() x (alt_grid.size() - 1) x (>= legendre_polynomial_dimension) x 6.\nupward_boundary_condition.shape():                        "sv,
                         v.upward_boundary_condition.shape(),
                         " - should be freq_grid.size() x 2 x fourier_mode_dimension x (quadrature_dimension / 2).\ndownward_boundary_condition.shape():                      "sv,
                         v.downward_boundary_condition.shape(),
                         " - should be freq_grid.size() x 2 x fourier_mode_dimension x (quadrature_dimension / 2)."sv);
    }

    const auto sep = tags.sep();
    return tags.format(ctx,
                       v.quadrature_dimension,
                       sep,
                       v.legendre_polynomial_dimension,
                       sep,
                       v.fourier_mode_dimension,
                       sep,
                       v.freq_grid,
                       sep,
                       v.alt_grid,
                       sep,
                       v.solar_azimuth_angle,
                       sep,
                       v.solar_zenith_angle,
                       sep,
                       v.solar_source,
                       sep,
                       v.bidirectional_reflectance_distribution_functions,
                       sep,
                       v.optical_thicknesses,
                       sep,
                       v.single_scattering_albedo,
                       sep,
                       v.source_polynomial,
                       sep,
                       v.legendre_coefficients,
                       sep,
                       v.upward_boundary_condition,
                       sep,
                       v.downward_boundary_condition);
  }
};

template <> struct xml_io_stream<VDisortSettings> {
  static constexpr std::string_view type_name = "VDisortSettings"sv;

  static void write(std::ostream&          os,
                    const VDisortSettings& x,
                    bofstream*             pbofs = nullptr,
                    std::string_view       name  = ""sv);

  static void read(std::istream& is, VDisortSettings& x, bifstream* pbifs = nullptr);
};
