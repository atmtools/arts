#pragma once

#include <enumsSpectralRadianceUnitType.h>
#include <path_point.h>
#include <rtepack.h>

struct SpectralRadianceTransformOperator {
  using Op = CustomOperator<void,
                            StokvecVector&,
                            StokvecMatrix&,
                            const AscendingGrid&,
                            const PropagationPathPoint&>;

  Op f{};
  std::string type{"CustomOperator"};

  SpectralRadianceTransformOperator(SpectralRadianceUnitType);
  SpectralRadianceTransformOperator(const std::string_view&);

  SpectralRadianceTransformOperator();
  SpectralRadianceTransformOperator(const SpectralRadianceTransformOperator& x);
  SpectralRadianceTransformOperator(
      SpectralRadianceTransformOperator&& x) noexcept;
  SpectralRadianceTransformOperator& operator=(
      const SpectralRadianceTransformOperator& x);
  SpectralRadianceTransformOperator& operator=(
      SpectralRadianceTransformOperator&& x) noexcept;

  SpectralRadianceTransformOperator(
      SpectralRadianceTransformOperator::Op&& x) noexcept;
  SpectralRadianceTransformOperator(
      const SpectralRadianceTransformOperator::Op& x);

  void operator()(StokvecVector& spectral_radiance,
                  StokvecMatrix& spectral_radiance_jacobian,
                  const AscendingGrid& frequency_grid,
                  const PropagationPathPoint& ray_path_point) const;
};

SpectralRadianceTransformOperator::Op from(SpectralRadianceUnitType x);

template <>
struct std::formatter<SpectralRadianceTransformOperator> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SpectralRadianceTransformOperator& op,
                              FmtContext& ctx) const {
    return std::format_to(
        ctx.out(), "<SpectralRadianceTransformOperator::{0}>"sv, op.type);
  }
};

template <>
struct xml_io_stream<SpectralRadianceTransformOperator> {
  static constexpr std::string_view type_name =
      "SpectralRadianceTransformOperator"sv;

  static void write(std::ostream& os,
                    const SpectralRadianceTransformOperator& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SpectralRadianceTransformOperator& x,
                   bifstream* pbifs = nullptr);
};
