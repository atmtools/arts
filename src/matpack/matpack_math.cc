#include "matpack_eigen.h"
#include "matpack_math.h"

#include <algorithm>

void mult_fast(ExhaustiveMatrixView A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(MatrixView A, const ConstMatrixView& B, const ConstMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveComplexMatrixView A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(ComplexMatrixView A, const ConstComplexMatrixView& B, const ConstComplexMatrixView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveVectorView A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(VectorView A, const ConstMatrixView& B, const ConstVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_fast(ExhaustiveComplexVectorView A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

void mult_slow(ComplexVectorView A, const ConstComplexMatrixView& B, const ConstComplexVectorView& C) {
  using namespace matpack::eigen;
  as_eigen(A).noalias() = B * C;
}

Vector uniform_grid(Numeric x0, Index N, Numeric dx) {
  Vector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto xd=x; x+=dx; return xd;});
  return out;
}

ComplexVector uniform_grid(Complex x0, Index N, Complex dx) {
  ComplexVector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto xd=x; x+=dx; return xd;});
  return out;
}

Vector diagonal(const ConstMatrixView& A) {
  using namespace matpack::eigen;
  return as_eigen(A).diagonal();
}

void cross3(VectorView c, const ConstVectorView& a, const ConstVectorView& b) {
  ARTS_ASSERT(a.nelem() == 3, a.nelem(), " vs 3");
  ARTS_ASSERT(b.nelem() == 3, b.nelem(), " vs 3");
  ARTS_ASSERT(c.nelem() == 3, c.nelem(), " vs 3");

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void cross3(ComplexVectorView c, const ConstComplexVectorView& a, const ConstComplexVectorView& b) {
  ARTS_ASSERT(a.nelem() == 3, a.nelem(), " vs 3");
  ARTS_ASSERT(b.nelem() == 3, b.nelem(), " vs 3");
  ARTS_ASSERT(c.nelem() == 3, c.nelem(), " vs 3");

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
