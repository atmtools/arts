#pragma once

#include <enumsHitranLineStrengthOption.h>
#include <enumsHitranZeemanOption.h>
#include <matpack.h>
#include <quantum.h>

#include "lbl_data.h"

namespace lbl {
struct hitran_record {
  QuantumIdentifier qid;
  Numeric           f0;
  Numeric           S;
  Numeric           A;
  Numeric           gamma_air;
  Numeric           gamma_self;
  Numeric           E;
  Numeric           n;
  Numeric           delta;
  Numeric           g_upp;
  Numeric           g_low;

  /** Convert to an absorption line
   *
   * @param[in] ls Whether the line strength or the Einstein coefficient is used
   * @param[in] local The local quantum numbers of the line
   * @param[in] zeeman How the Zeeman g-factors are set
   */
  [[nodiscard]] line from(HitranLineStrengthOption ls, QuantumState&& local, HitranZeemanOption zeeman) const;
};
using hitran_data = std::vector<hitran_record>;

hitran_data read_hitran_par(std::istream&                            file,
                            const std::vector<HitranFileFormatType>& format_order,
                            const Vector2&                           frequency_range);
hitran_data read_hitran_par(std::istream&&                           file,
                            const std::vector<HitranFileFormatType>& format_order,
                            const Vector2&                           frequency_range);
}  // namespace lbl

template <> struct std::formatter<lbl::hitran_record> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const lbl::hitran_record& v, FmtContext& ctx) const {
    const auto sep = tags.sep();
    tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx,
                v.qid,
                sep,
                v.f0,
                sep,
                v.S,
                sep,
                v.A,
                sep,
                v.gamma_air,
                sep,
                v.gamma_self,
                sep,
                v.E,
                sep,
                v.n,
                sep,
                v.delta,
                sep,
                v.g_upp,
                sep,
                v.g_low);
    tags.add_if_bracket(ctx, "]"sv);
    return ctx.out();
  }
};
