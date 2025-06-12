#include "matpack_mdspan_helpers_eigen.h"

namespace matpack::eigen {
#define EIGEN_STRIDED_MAT(U, T)                                            \
  Eigen::Map<                                                              \
      Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,   \
      Eigen::Unaligned,                                                    \
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>                       \
  mat(T x) {                                                               \
    return Eigen::Map<                                                     \
        Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, \
        Eigen::Unaligned,                                                  \
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(                    \
        const_cast<U *>(x.data_handle()),                                  \
        x.extent(0),                                                       \
        x.extent(1),                                                       \
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{x.stride(0),         \
                                                      x.stride(1)});       \
  }
EIGEN_STRIDED_MAT(Numeric, StridedMatrixView &)
EIGEN_STRIDED_MAT(Numeric, StridedConstMatrixView &)
EIGEN_STRIDED_MAT(Numeric, const StridedMatrixView &)
EIGEN_STRIDED_MAT(Numeric, const StridedConstMatrixView &)
EIGEN_STRIDED_MAT(Complex, StridedComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, StridedConstComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, const StridedComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, const StridedConstComplexMatrixView &)
#undef EIGEN_STRIDED_MAT

#define EIGEN_MAT(U, T)                                                     \
  Eigen::Map<                                                               \
      Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>    \
  mat(T x) {                                                                \
    return Eigen::Map<                                                      \
        Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>( \
        const_cast<U *>(x.data_handle()), x.extent(0), x.extent(1));        \
  }
EIGEN_MAT(Numeric, Matrix &)
EIGEN_MAT(Numeric, MatrixView &)
EIGEN_MAT(Numeric, ConstMatrixView &)
EIGEN_MAT(Numeric, const Matrix &)
EIGEN_MAT(Numeric, const MatrixView &)
EIGEN_MAT(Numeric, const ConstMatrixView &)
EIGEN_MAT(Complex, ComplexMatrix &)
EIGEN_MAT(Complex, ComplexMatrixView &)
EIGEN_MAT(Complex, ConstComplexMatrixView &)
EIGEN_MAT(Complex, const ComplexMatrix &)
EIGEN_MAT(Complex, const ComplexMatrixView &)
EIGEN_MAT(Complex, const ConstComplexMatrixView &)
#undef EIGEN_MAT

#define EIGEN_STRIDED_VEC(U, T)                                         \
  Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>,                       \
             Eigen::Unaligned,                                          \
             Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>             \
  col_vec(T x) {                                                        \
    return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>,              \
                      Eigen::Unaligned,                                 \
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(   \
        const_cast<U *>(x.data_handle()),                               \
        1,                                                              \
        x.extent(0),                                                    \
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)}); \
  }                                                                     \
  Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>,                       \
             Eigen::Unaligned,                                          \
             Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>             \
  row_vec(T x) {                                                        \
    return Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>,              \
                      Eigen::Unaligned,                                 \
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(   \
        const_cast<U *>(x.data_handle()),                               \
        x.extent(0),                                                    \
        1,                                                              \
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)}); \
  }
EIGEN_STRIDED_VEC(Numeric, StridedVectorView &)
EIGEN_STRIDED_VEC(Numeric, StridedConstVectorView &)
EIGEN_STRIDED_VEC(Numeric, const StridedVectorView &)
EIGEN_STRIDED_VEC(Numeric, const StridedConstVectorView &)
EIGEN_STRIDED_VEC(Complex, StridedComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, StridedConstComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, const StridedComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, const StridedConstComplexVectorView &)
#undef EIGEN_STRIDED_VEC

#define EIGEN_VEC(U, T)                                          \
  Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>> col_vec(T x) { \
    return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>>(      \
        const_cast<U *>(x.data_handle()), 1, x.extent(0));       \
  }                                                              \
  Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>> row_vec(T x) { \
    return Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>>(      \
        const_cast<U *>(x.data_handle()), x.extent(0), 1);       \
  }
EIGEN_VEC(Numeric, Vector &)
EIGEN_VEC(Numeric, VectorView &)
EIGEN_VEC(Numeric, ConstVectorView &)
EIGEN_VEC(Numeric, const Vector &)
EIGEN_VEC(Numeric, const VectorView &)
EIGEN_VEC(Numeric, const ConstVectorView &)
EIGEN_VEC(Complex, ComplexVector &)
EIGEN_VEC(Complex, ComplexVectorView &)
EIGEN_VEC(Complex, ConstComplexVectorView &)
EIGEN_VEC(Complex, const ComplexVector &)
EIGEN_VEC(Complex, const ComplexVectorView &)
EIGEN_VEC(Complex, const ConstComplexVectorView &)
#undef EIGEN_VEC
}  // namespace matpack::eigen
