// --------------//
//  Arts Vector  //
// ------------- //

// Included by arts_wrapper.h: definitions must be safe in multiple translation units.
#include <matpack.h>

inline auto ArtsVector::rows() const -> Index { return this->size(); }

inline auto ArtsVector::operator()(Index i) const -> Numeric {
  return this->elem_at(i);
}

inline auto ArtsVector::operator()(Index i) -> Numeric & { return this->elem_at(i); }

inline auto ArtsVector::data_pointer() -> Numeric * { return this->data_handle(); }

inline auto ArtsVector::data_pointer() const -> const Numeric * {
  return this->data_handle();
}

inline auto ArtsVector::accumulate(const ArtsVector &w) -> void {
  this->operator+=(w);
}

inline auto ArtsVector::subtract(const ArtsVector &w) -> void { this->operator-=(w); }

inline auto ArtsVector::scale(Numeric c) -> void { this->operator*=(c); }

inline auto ArtsVector::norm() const -> Numeric {
  return std::sqrt(dot(*this, *this));
}

inline Numeric dot(const ArtsVector &v, const ArtsVector &w) {
  Numeric x{};
  for (Size i = 0; i < v.size(); i++) {
    x += v[i] * w[i];
  }
  return x;
}

//-----------------//
//   Arts Matrix   //
//-----------------//

inline ArtsMatrix::ArtsMatrix(const Matrix &A) : Matrix(A) {
  // Nothing to do here.
}

template <typename ArtsType>
inline ArtsMatrix::ArtsMatrix(const ArtsMatrixReference<ArtsType> &A)
    : Matrix(static_cast<const ArtsType &>(A)) {
  // Nothing to do here.
}

inline auto ArtsMatrix::rows() const -> Index { return this->nrows(); }

inline auto ArtsMatrix::cols() const -> Index { return this->ncols(); }

inline auto ArtsMatrix::operator()(Index i, Index j) -> RealType & {
  return this->operator[](i, j);
}

inline auto ArtsMatrix::operator()(Index i, Index j) const -> RealType {
  return this->operator[](i, j);
}

inline auto ArtsMatrix::data_pointer() -> Numeric * { return this->data_handle(); }

inline void ArtsMatrix::accumulate(const MatrixType &B) { this->operator+=(B); }

inline void ArtsMatrix::accumulate(const ArtsCovarianceMatrixWrapper &B) {
  if (B.is_inverse()) {
    ::add_inv(*this, B);
  } else {
    MatrixView{ *this} += ConstMatrixView{B};
  }
}

inline auto ArtsMatrix::multiply(const ArtsCovarianceMatrixWrapper &B) -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(this->nrows(), B.cols());
  if (B.is_inverse()) {
    ::mult_inv(C, *this, B);
  } else {
    ::mult(C, *this, B.get_covmat());
  }
  return C;
}

inline void ArtsMatrix::subtract(const ArtsMatrix &B) { this->operator-=(B); }

inline auto ArtsMatrix::multiply(const ArtsMatrix &B) const -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(this->nrows(), B.ncols());
  mult(MatrixView{C}, MatrixView{*this}, MatrixView{B});
  return C;
}

inline auto ArtsMatrix::multiply(const ArtsVector &v) const -> ArtsVector {
  ArtsVector w;
  w.resize(this->nrows());
  mult(VectorView{w}, MatrixView{*this}, VectorView{v});
  return w;
}

inline auto ArtsMatrix::transpose_multiply(const ArtsMatrix &B) const -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(this->ncols(), B.ncols());
  ::mult(MatrixView{C}, matpack::transpose(MatrixView{*this}), MatrixView{B});
  return C;
}

inline auto ArtsMatrix::transpose_multiply(const ArtsVector &v) const -> ArtsVector {
  ArtsVector w;
  w.resize(this->ncols());
  ::mult(VectorView{w}, matpack::transpose(MatrixView{*this}), VectorView{v});
  return w;
}

inline auto ArtsMatrix::transpose_multiply_block(const ArtsVector &v,
                                          unsigned int start,
                                          unsigned int extent) const
    -> ArtsVector {
  ArtsVector w;
  w.resize(this->ncols());
  ConstVectorView v_view = v[Range(start, extent)];
  ::mult(VectorView{w}, matpack::transpose(ConstMatrixView{*this}), v_view);
  return w;
}

inline auto ArtsMatrix::solve(const VectorType &v) const -> ArtsVector {
  VectorType w;
  w.resize(this->nrows());
  ::solve(w, *this, v);
  return w;
}

inline auto ArtsMatrix::invert() const -> ArtsMatrix {
  ArtsMatrix B;
  B.resize(this->nrows(), this->ncols());
  ::inv(B, *this);
  return B;
}

inline void ArtsMatrix::scale(Numeric c) { this->operator*=(c); }

inline auto ArtsMatrix::transpose() const -> ArtsMatrix {
  ArtsMatrix B;
  B.Matrix::operator=(matpack::transpose(ConstMatrixView{*this}));
  return B;
}

inline auto ArtsMatrix::transpose_view() const & -> StridedConstMatrixView {
  return matpack::transpose(ConstMatrixView{*this});
}

//---------------------------//
//   Arts Matrix Reference   //
//---------------------------//

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::rows() const -> Index {
  return A.get().nrows();
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::cols() const -> Index {
  return A.get().ncols();
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::operator()(unsigned int i,
                                               unsigned int j) const
    -> RealType {
  return A.ro(i, j);
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::multiply(const ArtsMatrix &B) const
    -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(A.get().nrows(), B.ncols());
  mult(StridedMatrixView{C}, StridedMatrixView{A}, StridedConstMatrixView{B});
  return C;
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::multiply(const ArtsVector &v) const
    -> ArtsVector {
  ArtsVector w;
  w.resize(A.get().nrows());
  mult(StridedVectorView{w}, StridedMatrixView{A}, StridedConstVectorView{v});
  return w;
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::transpose_multiply(
    const ArtsVector &v) const -> ArtsVector {
  ArtsVector w;
  w.resize(A.get().ncols());
  ::mult(VectorView{w}, matpack::transpose(A.get()), ConstVectorView{v});
  return w;
}

template <>
inline auto ArtsMatrixReference<const Sparse>::transpose_multiply(
    const ArtsVector &v) const -> ArtsVector {
  ArtsVector w;
  w.resize(A.get().ncols());
  ::transpose_mult(w, A.get(), v);
  return w;
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::transpose_multiply(
    const ArtsMatrix &B) const -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(A.get().ncols(), B.ncols());
  ::mult(MatrixView{C}, matpack::transpose(A.get()), ConstMatrixView{B});
  return C;
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::multiply_add(
    const ArtsMatrix &B, const ArtsCovarianceMatrixWrapper &C) const -> ArtsMatrix {
  ArtsMatrix result;
  result.resize(A.get().nrows(), B.ncols());
  static_cast<Matrix&>(result) = 0;
  if (C.is_inverse()) {
    ::add_inv(result, C);
  } else {
    StridedMatrixView{result} += C.get_covmat();
  }
  ::mult(StridedMatrixView{result}, StridedConstMatrixView{A.get()},
         StridedConstMatrixView{B}, 1.0, 1.0);
  return result;
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::transpose() const -> ArtsMatrix {
  return ArtsMatrix{Matrix{matpack::transpose(A.get())}};
}

template <typename ArtsType>
inline auto ArtsMatrixReference<ArtsType>::transpose_view() const -> StridedConstMatrixView
  requires requires(const ArtsType& value) { StridedConstMatrixView{value}; } {
  return matpack::transpose(StridedConstMatrixView{A.get()});
}

//---------------------------//
//   Arts Covariance Matrix  //
//---------------------------//

inline auto ArtsCovarianceMatrixWrapper::rows() const -> Index {
  return covmat_.nrows();
}

inline auto ArtsCovarianceMatrixWrapper::cols() const -> Index {
  return covmat_.ncols();
}

inline auto ArtsCovarianceMatrixWrapper::multiply(const ArtsVector &v) const
    -> ArtsVector {
  ArtsVector w;
  w.resize(covmat_.nrows());
  if (is_inverse_) {
    ::mult_inv(w.view_as(w.size(), 1), covmat_, v.view_as(v.size(), 1));
  } else {
    ::mult(w, covmat_, v);
  }
  return w;
}

inline auto ArtsCovarianceMatrixWrapper::multiply(const ArtsMatrix &B) const
    -> ArtsMatrix {
  return multiply(StridedConstMatrixView{B});
}

inline auto ArtsCovarianceMatrixWrapper::multiply(StridedConstMatrixView B) const
    -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(covmat_.nrows(), B.ncols());
  if (is_inverse_) {
    ::mult_inv(C, covmat_, B);
  } else {
    ::mult(C, covmat_, B);
  }
  return C;
}

inline auto ArtsCovarianceMatrixWrapper::transpose_multiply(const ArtsVector &v) const
    -> ArtsVector {
  ArtsVector w;
  w.resize(covmat_.ncols());
  if (is_inverse_) {
    ::mult_inv(w.view_as(w.size(), 1), covmat_, v.view_as(v.size(), 1));
  } else {
    ::mult(w, covmat_, v);
  }
  return w;
}

inline auto ArtsCovarianceMatrixWrapper::transpose_multiply(const ArtsMatrix &B) const
    -> ArtsMatrix {
  ArtsMatrix C;
  C.resize(covmat_.ncols(), B.ncols());
  if (is_inverse_) {
    ::mult_inv(C, covmat_, B);
  } else {
    ::mult(C, covmat_, B);
  }
  return C;
}
