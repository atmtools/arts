#pragma once

#include <atm.h>
#include <operators.h>
#include <rtepack.h>

#include "bulk_scattering_properties.h"

struct ScatteringTroSpectralVector {
  using general_t = scattering::BulkScatteringProperties<
      scattering::Format::TRO,
      scattering::Representation::Spectral>;
  using gridded_t =
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Gridded>;

  std::optional<SpecmatMatrix> phase_matrix;
  PropmatVector extinction_matrix;
  StokvecVector absorption_vector;

  ScatteringTroSpectralVector& operator+=(
      const ScatteringTroSpectralVector& other);

  [[nodiscard]]
  static scattering::PhaseMatrixData<Numeric,
                                     scattering::Format::TRO,
                                     scattering::Representation::Spectral>
  to_general(const SpecmatMatrix&,
             const std::shared_ptr<Vector>& f);

  [[nodiscard]]
  static scattering::ExtinctionMatrixData<Numeric,
                                          scattering::Format::TRO,
                                          scattering::Representation::Spectral>
  to_general(const PropmatVector&, const std::shared_ptr<Vector>& f);

  [[nodiscard]]
  static scattering::AbsorptionVectorData<Numeric,
                                          scattering::Format::TRO,
                                          scattering::Representation::Spectral>
  to_general(const StokvecVector&, const std::shared_ptr<Vector>&);

  [[nodiscard]]
  general_t to_general(const std::shared_ptr<Vector>& f = nullptr) const;

  [[nodiscard]] gridded_t to_lab_frame(
      std::shared_ptr<const Vector> za_inc_grid,
      std::shared_ptr<const Vector> delta_aa_grid,
      std::shared_ptr<const scattering::ZenithAngleGrid> za_scat_grid_new)
      const;

  template <scattering::Format format>
  [[nodiscard]]
  scattering::BulkScatteringProperties<format,
                                       scattering::Representation::Spectral>
  to_spectral() const {
    return to_general().to_spectral();
  }

  template <scattering::Format format>
  [[nodiscard]]
  scattering::BulkScatteringProperties<format,
                                       scattering::Representation::Spectral>
  to_spectral(Index degree, Index order) const {
    return to_general().to_spectral(degree, order);
  };
};

using ScatteringGeneralSpectralTROFunc =
    CustomOperator<ScatteringTroSpectralVector,
                   const AtmPoint&,
                   const Vector&,
                   Index>;

struct ScatteringGeneralSpectralTRO {
  ScatteringGeneralSpectralTROFunc f{};

  [[nodiscard]]
  ScatteringTroSpectralVector get_bulk_scattering_properties_tro_spectral(
      const AtmPoint& atm_point, const Vector& f_grid, Index degree) const;
};

template <>
struct std::formatter<ScatteringTroSpectralVector> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ScatteringTroSpectralVector& v,
                              FmtContext& ctx) const {
    const auto sep = tags.sep();

    if (v.phase_matrix.has_value()) {
      tags.format(ctx, v.phase_matrix.value(), sep);
    }

    return tags.format(ctx, v.extinction_matrix, sep, v.absorption_vector);
  }
};

template <>
struct std::formatter<ScatteringGeneralSpectralTRO> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ScatteringGeneralSpectralTRO&,
                              FmtContext& ctx) const {
    return tags.format(ctx, "GeneralSpectralTRO"sv);
  }
};
