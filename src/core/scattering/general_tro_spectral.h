#pragma once

#include <atm.h>
#include <operators.h>
#include <rtepack.h>

#include "bulk_scattering_properties.h"

/** The bulk Legendre series of totally randomly oriented scatterers, as solvers take it
 *
 * phase_matrix [f, l] holds, for every frequency, the coefficient a_l on the
 * orthonormal Y_l0 of each element of the scattering-plane phase matrix (see
 * tro_legendre.h): a compact_planar_muelmat per degree l.  extinction_matrix
 * [f] and absorption_vector [f] hold K11 as A and a1 as I.
 */
struct ScatteringTroSpectralVector {
  using general_t = scattering::BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Spectral>;
  using gridded_t = scattering::BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>;

  std::optional<CompactPlanarMuelmatMatrix> phase_matrix;
  CompactPlanarPropmatVector                extinction_matrix;
  CompactPlanarStokvecVector                absorption_vector;

  ScatteringTroSpectralVector& operator+=(const ScatteringTroSpectralVector& other);

  [[nodiscard]]
  static scattering::PhaseMatrixData<scattering::Format::TRO, scattering::Representation::Spectral> to_general(
      const CompactPlanarMuelmatMatrix&, const std::shared_ptr<Vector>& f);

  [[nodiscard]]
  static scattering::ExtinctionMatrixData<scattering::Format::TRO, scattering::Representation::Spectral> to_general(
      const CompactPlanarPropmatVector&, const std::shared_ptr<Vector>& f);

  [[nodiscard]]
  static scattering::AbsorptionVectorData<scattering::Format::TRO, scattering::Representation::Spectral> to_general(
      const CompactPlanarStokvecVector&, const std::shared_ptr<Vector>&);

  [[nodiscard]]
  general_t to_general(const std::shared_ptr<Vector>& f = nullptr) const;

  [[nodiscard]] gridded_t to_lab_frame(std::shared_ptr<const Vector>                      za_inc_grid,
                                       std::shared_ptr<const Vector>                      delta_aa_grid,
                                       std::shared_ptr<const scattering::ZenithAngleGrid> za_scat_grid_new) const;
};

using ScatteringGeneralSpectralTROFunc =
    CustomOperator<ScatteringTroSpectralVector, const AtmPoint&, const Vector&, Index>;

/** A scattering species given by a user's function of its TRO Legendre series
 *
 * f(atm_point, f_grid, degree) must give the bulk Legendre series to degree
 * at f_grid (see ScatteringTroSpectralVector), with the extinction and
 * absorption at f_grid; get_bulk_scattering_properties_tro_spectral throws
 * if it gives anything else.
 *
 * degree is where the series is truncated when the caller does not choose
 * it, as for the azimuthal Fourier modes; it is unset while negative.
 */
struct ScatteringGeneralSpectralTRO {
  ScatteringGeneralSpectralTROFunc f{};
  Index                            degree{-1};

  [[nodiscard]]
  ScatteringTroSpectralVector get_bulk_scattering_properties_tro_spectral(const AtmPoint& atm_point,
                                                                          const Vector&   f_grid,
                                                                          Index           degree) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame phase matrix of the series to degree
   *
   * See PhaseMatrixData<TRO, Spectral>::to_lab_frame_fourier_modes; the
   * series is evaluated exactly at the given zenith angles [deg].
   */
  [[nodiscard]]
  scattering::BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
  get_bulk_scattering_properties_aro_fourier(const AtmPoint& atm_point,
                                             const Vector&   f_grid,
                                             const Vector&   za_inc_grid,
                                             const Vector&   za_scat_grid,
                                             Index           max_mode) const;
};

template <> struct std::formatter<ScatteringTroSpectralVector> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const ScatteringTroSpectralVector& v, FmtContext& ctx) const {
    const auto sep = tags.sep();

    if (v.phase_matrix.has_value()) { tags.format(ctx, v.phase_matrix.value(), sep); }

    return tags.format(ctx, v.extinction_matrix, sep, v.absorption_vector);
  }
};

template <> struct std::formatter<ScatteringGeneralSpectralTRO> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const ScatteringGeneralSpectralTRO&, FmtContext& ctx) const {
    return tags.format(ctx, "GeneralSpectralTRO"sv);
  }
};

template <> struct xml_io_stream<ScatteringGeneralSpectralTRO> {
  static constexpr std::string_view type_name = "ScatteringGeneralSpectralTRO";

  static void write(std::ostream&, const ScatteringGeneralSpectralTRO&, bofstream* = nullptr, std::string_view = ""sv);

  static void read(std::istream&, ScatteringGeneralSpectralTRO&, bifstream* = nullptr);
};
