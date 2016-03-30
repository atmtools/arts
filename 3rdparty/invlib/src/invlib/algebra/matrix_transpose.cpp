template<typename T1>
MatrixTranspose<T1>::MatrixTranspose(T1 A_)
    : A(A_)
{
    // Nothing to do here.
}

template<typename T1>
auto MatrixTranspose<T1>::multiply(const VectorType &v) const
    -> VectorType
{
    MatrixType C = A;
    return C.transpose_multiply(v);
}

template<typename T1>
auto MatrixTranspose<T1>::multiply(const MatrixType &B) const
    -> MatrixType
{
    MatrixType C = A;
    return C.transpose_multiply(B);
}

template<typename T1>
auto MatrixTranspose<T1>::invert() const
    -> MatrixType
{
    MatrixType C = ((MatrixType) A).transpose();
    return C.invert();
}

template<typename T1>
auto MatrixTranspose<T1>::solve(const VectorType &v) const
    -> VectorType
{
    MatrixType C = ((MatrixType) A).transpose();
    return C.solve(v);
}

template<typename T1>
    template<typename T2>
auto MatrixTranspose<T1>::operator+(T2 &&B) const
    -> Sum<T2>
{
    return Sum<T2>(*this, B);
}

template<typename T1>
    template<typename T2>
auto MatrixTranspose<T1>::operator*(T2 &&B) const
    -> Product<T2>
{
    return Product<T2>(*this, B);
}

template<typename T1>
MatrixTranspose<T1>::operator ResultType() const
{
    ResultType B = A.transpose();
    return B;
}
