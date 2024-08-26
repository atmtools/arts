#pragma once

#include <unordered_map>

#include "atm.h"
#include "covariance_matrix.h"
#include "enums.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "surf.h"

struct JacobianTargetsDiagonalCovarianceMatrixMap {
  std::unordered_map<JacobianTargetType, std::pair<BlockMatrix, BlockMatrix>>
      map;

  [[nodiscard]] auto begin() { return map.begin(); }
  [[nodiscard]] auto end() { return map.end(); }
  [[nodiscard]] auto begin() const { return map.begin(); }
  [[nodiscard]] auto end() const { return map.end(); }
  [[nodiscard]] auto cbegin() const { return map.cbegin(); }
  [[nodiscard]] auto cend() const { return map.cend(); }
  [[nodiscard]] auto size() const { return map.size(); }
  [[nodiscard]] auto empty() const { return map.empty(); }
  [[nodiscard]] auto clear() { return map.clear(); }
  [[nodiscard]] auto find(const JacobianTargetType& x) const { return map.find(x); }

  void set(const JacobianTargetType& type,
           const BlockMatrix& matrix,
           const BlockMatrix& inverse = BlockMatrix{});

  void set(const AtmKeyVal& type,
           const BlockMatrix& matrix,
           const BlockMatrix& inverse = BlockMatrix{});
           
  void set(const SurfaceKeyVal& type,
           const BlockMatrix& matrix,
           const BlockMatrix& inverse = BlockMatrix{});
           
  void set(const LblLineKey& type,
           const BlockMatrix& matrix,
           const BlockMatrix& inverse = BlockMatrix{});

  friend std::ostream& operator<<(
      std::ostream& os, const JacobianTargetsDiagonalCovarianceMatrixMap&) {
    return os << "NotImplementedYet";
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
