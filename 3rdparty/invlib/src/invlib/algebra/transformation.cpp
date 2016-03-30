template<typename T1, typename Transform>
Transformation<T1, Transform>::Transformation(const T1 A_, Transform t_)
        : A(A_), t(std::forward<Transform>(t_))
{
    // Nothing to do here.
}

template<typename T1, typename Transform>
Transformation<T1, Transform>::operator MatrixType() const
{
    MatrixType B = A;
    t.apply_matrix(B);
    return B;
}

template<typename T1, typename Transform>
Transformation<T1, Transform>::operator VectorType() const
{
    VectorType v = A;
    t.apply_vector(v);
    return v;
}

template<typename T1, typename Transform>
auto Transformation<T1, Transform>::multiply(const VectorType & v) const
    -> VectorType
{
    MatrixType B = *this;
    return B.multiply(v);
}

template<typename T1, typename Transform>
auto Transformation<T1, Transform>::multiply(const MatrixType & B) const
    -> MatrixType
{
    MatrixType C = *this;
    return C.multiply(B);
}

template<typename T1, typename Transform>
auto Transformation<T1, Transform>::invert() const
    -> MatrixType
{
    MatrixType B = *this;
    return B.invert();
}

template<typename T1, typename Transform>
auto Transformation<T1, Transform>::solve(const VectorType & v) const
    -> VectorType
{
    MatrixType B = *this;
    return B.solve(v);
}

template<typename T1, typename Transform>
    template<typename T2>
auto Transformation<T1, Transform>::operator+(T2 &&B) const
    ->Sum<T2>
{
    return Sum<T2>(*this, B);
}

template<typename T1, typename Transform>
    template<typename T2>
auto Transformation<T1, Transform>::operator-(T2 &&B) const
    ->Difference<T2>
{
    return Difference<T2>(*this, B);
}

template<typename T1, typename Transform>
    template<typename T2>
auto Transformation<T1, Transform>::operator*(T2 &&B) const
    ->Product<T2>
{
    return Product<T2>(*this, B);
}

template <typename T1>
constexpr auto Identity::apply(T1 &&t)
    -> decltype(std::forward<T1>(t))
{
    return std::forward<T1>(t);
}

template <typename MatrixType>
    template <typename T1>
void NormalizeDiagonal<MatrixType>::apply_matrix(T1 &B) const
{
    unsigned int m, n;
    m = (unsigned int) B.rows();
    n = (unsigned int) B.cols();

    for (unsigned int i = 0; i < m; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            B(i,j) *= 1.0 / sqrt(A(i,i) * A(j,j));
        }
    }
}

template <typename MatrixType>
    template <typename T1>
void NormalizeDiagonal<MatrixType>::apply_vector(T1 &v) const
{
    unsigned int m;
    m = (unsigned int) v.rows();

    for (unsigned int i = 0; i < m; i++)
    {
        v(i) *= 1.0 / sqrt(A(i,i));
    }
}

template <typename MatrixType>
    template <typename T1>
auto NormalizeDiagonal<MatrixType>::apply(T1&& B)
    -> Transform<T1>
{
    return Transform<T1>(std::forward<T1>(B), *this);
}
