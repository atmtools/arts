// --------------//
//  Arts Vector  //
// ------------- //

#include "matpack_data.h"

auto ArtsVector::rows() const
    -> Index
{
    return this->nelem();
}

auto ArtsVector::operator()(Index i) const
    -> Numeric
{
    return this->elem_at(i);
}

auto ArtsVector::operator()(Index i)
    -> Numeric &
{
    return this->elem_at(i);
}

auto ArtsVector::data_pointer()
    -> Numeric *
{
    return this->data_handle();
}

auto ArtsVector::data_pointer() const
    -> const Numeric *
{
    return this->data_handle();
}

auto ArtsVector::accumulate(const ArtsVector& w)
    -> void
{
    this->operator+=(w);
}

auto ArtsVector::subtract(const ArtsVector& w)
    -> void
{
    this->operator-=(w);
}

auto ArtsVector::scale(Numeric c)
    -> void
{
    this->operator*=(c);
}

auto ArtsVector::norm() const
    -> Numeric
{
    return sqrt(ExhaustiveConstVectorView{*this} * ExhaustiveConstVectorView{*this});
}

Numeric dot(const ArtsVector& v, const ArtsVector& w)
{
    return ExhaustiveConstVectorView{v} * ExhaustiveConstVectorView{w};
}

//-----------------//
//   Arts Matrix   //
//-----------------//

ArtsMatrix::ArtsMatrix(const Matrix &A)
    : Matrix(A)
{
    // Nothing to do here.
}

template<typename ArtsType>
ArtsMatrix::ArtsMatrix(
    const ArtsMatrixReference<ArtsType> & A)
    : Matrix(static_cast<const ArtsType &>(A))
{
    // Nothing to do here.
}

auto ArtsMatrix::rows() const
    -> Index
{
    return this->nrows();
}

auto ArtsMatrix::cols() const
    -> Index
{
    return this->ncols();
}


auto ArtsMatrix::operator()(Index i, Index j)
    -> RealType &
{
    return this->elem_at(i, j);
}

auto ArtsMatrix::operator()(Index i, Index j) const
    -> RealType
{
    return this->elem_at(i, j);
}

auto ArtsMatrix::data_pointer()
    -> Numeric *
{
    return this->data_handle();
}

void ArtsMatrix::accumulate(const MatrixType & B)
{
    this->operator+=(B);
}

void ArtsMatrix::accumulate(const ArtsCovarianceMatrixWrapper & B)
{
    if (B.is_inverse()) {
        ::add_inv(*this, B);
    } else {
        ExhaustiveMatrixView{*this} += ConstMatrixView{B};
    }
}

auto ArtsMatrix::multiply(const ArtsCovarianceMatrixWrapper & B)
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(this->nrows(), B.cols());
    if (B.is_inverse()) {
        ::mult_inv(C, *this, B);
    } else {
        ::mult(C, *this, B.get_covmat());
    }
    return C;
}

void ArtsMatrix::subtract(const ArtsMatrix& B)
{
    this->operator-=(B);
}

auto ArtsMatrix::multiply(const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(this->nrows(), B.ncols());
    ::mult(ExhaustiveMatrixView{C},ExhaustiveMatrixView{*this}, ExhaustiveMatrixView{B});
    return C;
}

auto ArtsMatrix::multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->nrows());
    ::mult(ExhaustiveVectorView{w}, ExhaustiveMatrixView{*this}, ExhaustiveVectorView{v});
    return w;
}

auto ArtsMatrix::transpose_multiply(const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(this->ncols(), B.ncols());
    ::mult(ExhaustiveMatrixView{C}, ::transpose(ExhaustiveMatrixView{*this}), ExhaustiveMatrixView{B});
    return C;
}

auto ArtsMatrix::transpose_multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->ncols());
    ::mult(ExhaustiveVectorView{w}, ::transpose(ExhaustiveMatrixView{*this}), ExhaustiveVectorView{v});
    return w;
}

auto ArtsMatrix::transpose_multiply_block(const ArtsVector &v,
                                          unsigned int start,
                                          unsigned int extent) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->ncols());
    ConstVectorView v_view = v[Range(start, extent)];
    ::mult(ExhaustiveVectorView{w}, ::transpose(ExhaustiveConstMatrixView{*this}), v_view);
    return w;
}

auto ArtsMatrix::solve(const VectorType& v) const
    -> ArtsVector
{
    VectorType w; w.resize(this->nrows());
    ::solve(w, *this, v);
    return w;
}

auto ArtsMatrix::invert() const
    -> ArtsMatrix
{
    ArtsMatrix B; B.resize(this->nrows(), this->ncols());
    ::inv(B, *this);
    return B;
}

void ArtsMatrix::scale(Numeric c)
{
    this->operator*=(c);
}

auto ArtsMatrix::transpose() const
    -> ArtsMatrix
{
    ArtsMatrix B;
    B.Matrix::operator=(::transpose(ExhaustiveConstMatrixView{*this}));
    return B;
}

//---------------------------//
//   Arts Matrix Reference   //
//---------------------------//

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::rows() const
    -> Index
{
    return A.get().nrows();
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::cols() const
    -> Index
{
    return A.get().ncols();
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::operator()(
    unsigned int i,
    unsigned int j) const
    -> RealType
{
    return A.ro(i, j);
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::multiply(
    const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(A.get().nrows(), B.ncols());
    ::mult(ExhaustiveMatrixView{C}, MatrixView{A}, ExhaustiveConstMatrixView{B});
    return C;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.get().nrows());
    ::mult(ExhaustiveVectorView{w}, MatrixView{A}, ExhaustiveConstVectorView{v});
    return w;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::transpose_multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.get().ncols());
    ::mult(ExhaustiveVectorView{w}, ::transpose(A.get()), ExhaustiveConstVectorView{v});
    return w;
}

template<>
auto ArtsMatrixReference<const Sparse>::transpose_multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.get().ncols());
    ::transpose_mult(w, A.get(), v);
    return w;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::transpose_multiply(
    const ArtsMatrix & B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(A.get().ncols(), B.ncols());
    ::mult(ExhaustiveMatrixView{C}, ::transpose(A.get()), ExhaustiveConstMatrixView{B});
    return C;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::transpose() const
    -> ConstMatrixView
{
    return ::transpose(A.get());
}

//---------------------------//
//   Arts Covariance Matrix  //
//---------------------------//

auto ArtsCovarianceMatrixWrapper::rows() const
    -> Index
{
    return covmat_.nrows();
}

auto ArtsCovarianceMatrixWrapper::cols() const
    -> Index
{
    return covmat_.ncols();
}

auto ArtsCovarianceMatrixWrapper::multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(covmat_.nrows());
    if (is_inverse_) {
        ::mult_inv(MatrixView{w}, covmat_, ConstMatrixView{v});
    } else {
        ::mult(w, covmat_, v);
    }
    return w;
}

auto ArtsCovarianceMatrixWrapper::multiply(
    const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(covmat_.nrows(), B.ncols());
    if (is_inverse_) {
        ::mult_inv(C, covmat_, B);
    } else {
        ::mult(C, covmat_, B);
    }
    return C;
}

auto ArtsCovarianceMatrixWrapper::transpose_multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(covmat_.ncols());
    if (is_inverse_) {
        ::mult_inv(MatrixView{w}, covmat_, ConstMatrixView{v});
    } else {
        ::mult(w, covmat_, v);
    }
    return w;
}

auto ArtsCovarianceMatrixWrapper::transpose_multiply(
    const ArtsMatrix & B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(covmat_.ncols(), B.ncols());
    if (is_inverse_) {
        ::mult_inv(C, covmat_, B);
    } else {
        ::mult(C, covmat_, B);
    }
    return C;
}
