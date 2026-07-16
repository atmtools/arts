#pragma once

#include <jpl_species.h>
#include <matpack.h>
#include <quantum.h>

#include "lbl_data.h"

namespace lbl {
struct jpl_record {
  static constexpr Numeric T0 = 300.0;

  Jpl::LineDataMod jpl_id;  // ID of the species

  // JPL format in order [F13.4,2F8.4,I2,F10.4,I3,I7,I4,12I2]
  Numeric f0;     // Central frequency
  Numeric df;     // Error central frequency
  Numeric s;      // Line intensity
  Index   dr;     // Degree of freedom
  Numeric E;      // Energy of lower state
  Index   g_upp;  // Upper state degeneracy
  Index   qnfmt;  // Quantum number format

  // Create a line with fixed pressure broadening and no local quantum state.
  [[nodiscard]] line from() const;
};
using jpl_data = std::vector<jpl_record>;

jpl_data read_jpl_lines(std::istream& file);
jpl_data read_jpl_lines(std::istream&& file);
}  // namespace lbl

template <>
struct std::formatter<lbl::jpl_record> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::jpl_record& v, FmtContext& ctx) const {
    const auto sep = tags.sep();
    tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx, v.jpl_id, sep, v.f0, sep, v.df, sep, v.s, sep, v.dr, sep, v.E, sep, v.g_upp, sep, v.qnfmt);
    tags.add_if_bracket(ctx, "]"sv);
    return ctx.out();
  }
};
