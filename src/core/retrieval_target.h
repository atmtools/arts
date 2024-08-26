#pragma once

#include <unordered_map>

#include "covariance_matrix.h"
#include "jacobian.h"

struct JacobianTargetsDiagonalCovarianceMatrixMap {
  std::unordered_map<JacobianTargetType, std::pair<BlockMatrix, BlockMatrix>>
      map;

  void set(const JacobianTargetType& type,
           const BlockMatrix& forward,
           const BlockMatrix& inverse = BlockMatrix{});
};
