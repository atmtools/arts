#pragma once

#include <isotopologues.h>
#include <species.h>
#include <xml.h>

#include <unordered_map>

#include "lbl_temperature_model.h"

namespace lbl::linemixing {
struct species_data {
  temperature::data scaling{LineShapeModelType::T0, {0}};
  temperature::data beta{LineShapeModelType::T0, {0}};
  temperature::data lambda{LineShapeModelType::T0, {0}};
  temperature::data collisional_distance{LineShapeModelType::T0, {0}};

  [[nodiscard]] Numeric Q(const Rational J,
                          const Numeric T,
                          const Numeric T0,
                          const Numeric energy) const;

  [[nodiscard]] Numeric Omega(const Numeric T,
                              const Numeric T0,
                              const Numeric mass,
                              const Numeric other_mass,
                              const Numeric energy_x,
                              const Numeric energy_xm2) const;
};  // species_data
}  // namespace lbl::linemixing

using LinemixingSingleEcsData = lbl::linemixing::species_data;
using LinemixingSpeciesEcsData =
    std::unordered_map<SpeciesEnum, LinemixingSingleEcsData>;
using LinemixingEcsData =
    std::unordered_map<SpeciesIsotope, LinemixingSpeciesEcsData>;

template <>
struct std::formatter<LinemixingSingleEcsData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::linemixing::species_data &v,
                              FmtContext &ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.scaling,
                sep,
                v.beta,
                sep,
                v.lambda,
                sep,
                v.collisional_distance);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct xml_io_stream<LinemixingSingleEcsData> {
  static constexpr std::string_view type_name = "LinemixingSingleEcsData"sv;

  static void write(std::ostream &os,
                    const LinemixingSingleEcsData &x,
                    bofstream *pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream &is,
                   LinemixingSingleEcsData &x,
                   bifstream *pbifs = nullptr);
};
