template<typename T1>
MatrixTranspose<T1>::MatrixTranspose(T1 A_)
    : A(A_)
{
    // Nothing to do here.
}

template<typename T1>
template<typename T2>
auto MatrixTranspose<T1>::multiply(const T2 &v) const
    -> typename T2::ResultType
{
    return remove_reference_wrapper(A).transpose_multiply(v);
}

template<typename T1>
auto MatrixTranspose<T1>::multiply(const VectorType &v) const
    -> VectorType
{
    return remove_reference_wrapper(A).transpose_multiply(v);
}

template<typename T1>
auto MatrixTranspose<T1>::multiply(const MatrixType &B) const
    -> MatrixType
{
    return remove_reference_wrapper(A).transpose_multiply(B);
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
auto MatrixTranspose<T1>::diagonal() const
    -> VectorType
{
    return A.diagonal();
}

template<typename T1>
auto MatrixTranspose<T1>::row(size_t i) const
    -> VectorType
{
    return remove_reference_wrapper(A).col(i);
}

template<typename T1>
auto MatrixTranspose<T1>::col(size_t i) const
    -> VectorType
{
    return remove_reference_wrapper(A).row(i);
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
    ResultType B = remove_reference_wrapper(A).transpose();
    return B;
}
