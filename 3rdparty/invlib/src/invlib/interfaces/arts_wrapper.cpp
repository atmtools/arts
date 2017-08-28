// --------------//
//  Arts Vector  //
// ------------- //

ArtsVector::ArtsVector(ArtsVector &&v)
{
    this->mrange = v.mrange;
    this->mdata  = v.mdata;
    v.mdata      = nullptr;
}

auto ArtsVector::operator=(ArtsVector &&v)
    -> ArtsVector &
{
    delete[] this->mdata;
    this->mrange  = v.mrange;
    this->mdata   = v.mdata;
    v.mdata       = nullptr;
    return *this;
}

auto ArtsVector::rows() const
    -> Index
{
    return this->nelem();
}

auto ArtsVector::operator()(Index i) const
    -> Numeric
{
    return this->get(i);
}

auto ArtsVector::operator()(Index i)
    -> Numeric &
{
    return this->get(i);
}

auto ArtsVector::data_pointer()
    -> Numeric *
{
    return this->mdata;
}

auto ArtsVector::data_pointer() const
    -> const Numeric *
{
    return this->mdata;
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
    return sqrt(operator*(*this, *this));
}

Numeric dot(const ArtsVector& v, const ArtsVector& w)
{
    return v * w;
}

//-----------------//
//   Arts Matrix   //
//-----------------//

ArtsMatrix::ArtsMatrix(const Matrix &A)
    : Matrix(A)
{
    // Nothing to do here.
}

ArtsMatrix::ArtsMatrix(ArtsMatrix &&A)
{
    this->mrr = A.mrr;
    this->mcr = A.mcr;
    this->mdata  = A.mdata;
    A.mdata = nullptr;
}

auto ArtsMatrix::operator=(ArtsMatrix &&A)
    -> ArtsMatrix &
{
    delete[] this->mdata;
    this->mcr  = A.mcr;
    this->mrr  = A.mrr;
    this->mdata   = A.mdata;
    A.mdata = nullptr;
    return *this;
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
    return this->get(i,j);
}

auto ArtsMatrix::operator()(Index i, Index j) const
    -> RealType
{
    return this->get(i,j);
}

auto ArtsMatrix::data_pointer()
    -> Numeric *
{
    return this->mdata;
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
        ::operator+=(*this, B);
    }
}

void ArtsMatrix::subtract(const ArtsMatrix& B)
{
    this->operator-=(B);
}

auto ArtsMatrix::multiply(const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(this->nrows(), B.ncols());
    ::mult(C, *this, B);
    return C;
}

auto ArtsMatrix::multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->nrows());
    ::mult(w, *this, v);
    return w;
}

auto ArtsMatrix::transpose_multiply(const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(this->ncols(), B.ncols());
    ::mult(C, ::transpose(*this), B);
    return C;
}

auto ArtsMatrix::transpose_multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->ncols());
    ::mult(w, ::transpose(*this), v);
    return w;
}

auto ArtsMatrix::transpose_multiply_block(const ArtsVector &v,
                                          unsigned int start,
                                          unsigned int extent) const
    -> ArtsVector
{
    ArtsVector w; w.resize(this->ncols());
    ConstVectorView v_view = v[Range(start, extent)];
    ::mult(w, ::transpose(*this), v_view);
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
    B.Matrix::operator=(::transpose(*this));
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
    ::mult(C, A, B);
    return C;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.get().nrows());
    ::mult(w, A, v);
    return w;
}

template<typename ArtsType>
auto ArtsMatrixReference<ArtsType>::transpose_multiply(
    const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.get().ncols());
    ::mult(w, ::transpose(A.get()), v);
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
    ::mult(C, ::transpose(A.get()), B);
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
        ::mult_inv(w, covmat_, v);
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
        ::mult_inv(w, covmat_, v);
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
