#include "retrieval_target.h"

void JacobianTargetsDiagonalCovarianceMatrixMap::set(
    const JacobianTargetType& type,
    const BlockMatrix& forward,
    const BlockMatrix& inverse) {
  map[type] = {forward, inverse};
}
