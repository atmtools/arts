#include "retrieval_target.h"

void JacobianTargetsDiagonalCovarianceMatrixMap::set(
    const JacobianTargetType& type,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  map[type] = {matrix, inverse};
}

void JacobianTargetsDiagonalCovarianceMatrixMap::set(
    const AtmKeyVal& type,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  set(JacobianTargetType{type}, matrix, inverse);
}

void JacobianTargetsDiagonalCovarianceMatrixMap::set(
    const SurfaceKeyVal& type,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  set(JacobianTargetType{type}, matrix, inverse);
}

void JacobianTargetsDiagonalCovarianceMatrixMap::set(
    const LblLineKey& type,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  set(JacobianTargetType{type}, matrix, inverse);
}
