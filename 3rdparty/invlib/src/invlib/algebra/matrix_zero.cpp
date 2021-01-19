template <typename T1>
inline auto MatrixZero::multiply(T1 &&A) const
    -> T1 &&
{
    A.scale(0.0);
    return std::forward<T1>(A);
}

inline auto MatrixZero::invert() const
    -> MatrixZero
{
    return MatrixZero();
}

template <typename T1>
inline T1 && MatrixZero::solve(T1 && v) const
{
    v.scale(0.0);
    return std::forward<T1>(v);
}

template<typename T1>
inline auto MatrixZero::operator+(T1 &&B) const
    -> T1 &&
{
    return std::forward<T1>(B);
}

template <typename T1>
inline auto MatrixZero::operator*(const T1 &A) const
    -> Product<T1>
{
    return Product<T1>(*this, A);
}

inline MatrixZero inv(const MatrixZero &)
{
    return MatrixZero{};
}

