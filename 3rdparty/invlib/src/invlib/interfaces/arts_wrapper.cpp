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

void ArtsMatrix::accumulate(const ArtsMatrix& B)
{
    this->operator+=(B);
}

void ArtsMatrix::accumulate(const ArtsSparse& B)
{
    ArtsMatrix C = B.operator ArtsMatrix();
    this->operator+=(C);
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

//-----------------//
//   Arts Sparse   //
//-----------------//

auto ArtsSparse::rows() const
    -> Index
{
    return A.nrows();
}

auto ArtsSparse::cols() const
    -> Index
{
    return A.ncols();
}

auto ArtsSparse::multiply(const ArtsMatrix &B) const
    -> ArtsMatrix
{
    ArtsMatrix C; C.resize(A.nrows(), B.ncols());
    ::mult(C, A, B);
    return C;
}

auto ArtsSparse::multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.nrows());
    ::mult(w, A, v);
    return w;
}

auto ArtsSparse::transpose_multiply(const ArtsVector &v) const
    -> ArtsVector
{
    ArtsVector w; w.resize(A.ncols());
    ::transpose_mult(w, A, v);
    return w;
}

ArtsSparse::operator ArtsMatrix() const
{
    ArtsMatrix B;
    B.Matrix::operator=(A.operator Matrix());
    return B;
}
