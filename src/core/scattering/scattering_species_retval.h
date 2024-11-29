#pragma once

#include <rtepack.h>

struct ScatteringTroSpectralVector {
  std::optional<ComplexMuelmatMatrix> phase_matrix;
  PropmatVector extinction_matrix;
  StokvecVector absorption_vector;

  ScatteringTroSpectralVector& operator+=(
      const ScatteringTroSpectralVector& other);
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