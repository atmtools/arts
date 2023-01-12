#include "matpack/matpack_eigen2.h"
#include "matpack/matpack_math2.h"

#include <algorithm>

void mult_fast(ExhaustiveMatrixView&& A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(MatrixView&& A, const ConstMatrixView& B, const ConstMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveComplexMatrixView&& A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(ComplexMatrixView&& A, const ConstComplexMatrixView& B, const ConstComplexMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveVectorView&& A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(VectorView&& A, const ConstMatrixView& B, const ConstVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveComplexVectorView&& A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(ComplexVectorView&& A, const ConstComplexMatrixView& B, const ConstComplexVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

Vector uniform_grid(Numeric x0, Index N, Numeric dx) {
  Vector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto out=x; x+=dx; return out;});
  return out;
}

ComplexVector uniform_grid(Complex x0, Index N, Complex dx) {
  ComplexVector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto out=x; x+=dx; return out;});
  return out;
}