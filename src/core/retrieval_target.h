#pragma once

#include <unordered_map>

#include "covariance_matrix.h"
#include <jacobian.h>

struct PairOfBlockMatrix {
  BlockMatrix first, second;
};

using JacobianTargetsDiagonalCovarianceMatrixMap =
    std::unordered_map<JacobianTargetType, PairOfBlockMatrix>;

template <>
struct std::formatter<PairOfBlockMatrix> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const PairOfBlockMatrix& x,
                              FmtContext& ctx) const {
    return tags.format(ctx, x.first, tags.sep(), x.second);
  }
};

template <>
struct std::formatter<JacobianTargetsDiagonalCovarianceMatrixMap> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const JacobianTargetsDiagonalCovarianceMatrixMap&,
                              FmtContext& ctx) const {
    return tags.format(ctx, "NotImplementedYet"sv);
  }
};
