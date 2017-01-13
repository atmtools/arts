#ifndef PRECISION_MATRIX_H
#define PRECISION_MATRIX_H

// -------------------- //
// Forward Declarations //
// -------------------- //

namespace invlib
{

template
<
typename Matrix
>
class PrecisionMatrix;

template
<
typename T1,
typename T2
>
class MatrixProduct;

template
<
typename T1,
typename T2
>
class MatrixSum;

template
<
typename Matrix
>
const Matrix & inv(PrecisionMatrix<Matrix> &A);

template
<
typename Matrix
>
const Matrix & inv(const PrecisionMatrix<Matrix> &A);

// --------------------- //
// Class PrecisionMatrix //
// --------------------- //

/**
 * \brief Wrapper to represent precision matrices.
 *
 * The computation of MAP estimators involves covariance matrices mainly
 * through their inverse. In some cases it may thus be beneficial to provide
 * the inverses directly. Using the PrecisionMatrix type to instantiate the
 * MAP class will interpret the matrices given to the MAP constructor as
 * precision matrices.
 *
 * Technically this is achieved by overloading the inv() function as well
 * as overriding all relevant member functions used for the matrix algebra.
 * The class provides the overloaded inv() function to return a reference to
 * the actual precision matrix provided. Operations involving the wrapper
 * directly, however, require inverting the matrix (or at least solving the
 * the correspongind linear system).
 *
 * Also provides an overloaded transp() function, that does nothing, since the
 * a precision matrix is by definition symmetric.
 *
 * \todo Maybe remove in favor of passing inv(S) as parameter to MAP directly.
 */
template
<
typename Matrix
>
class PrecisionMatrix
{

public:

    /*! The basic scalar type. */
    using RealType   = typename decay<Matrix>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<Matrix>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<Matrix>::MatrixType;
    /*! The type of the result of the expression */
    using ResultType = typename decay<Matrix>::ResultType;

    PrecisionMatrix(const Matrix& A_)
        : A(A_) {}

    // ------------------- //
    // Algebraic Operators //
    // ------------------- //

    template <typename T>
    using Product = MatrixProduct<PrecisionMatrix, T>;

    
    template <typename T>
    Product<T> operator*(T &&B) const
    {
        return Product<T>(*this, B);
    }

    template <typename T>
    using Sum = MatrixSum<PrecisionMatrix, T>;

    template <typename T>
    Product<T> operator+(T &&B) const
    {
        return Sum<T>(*this, B);
    }

    ResultType multiply(const ResultType& B) const
    {
        MatrixType tmp = A.invert();
        return tmp * B;
    }

    VectorType multiply(const VectorType& v) const
    {
        return A.solve(v);
    }

    operator ResultType() const
    {
        ResultType tmp = A.invert();
        return tmp;
    }

    friend const Matrix & inv<Matrix>(PrecisionMatrix &A);
    friend const Matrix & inv<Matrix>(const PrecisionMatrix &A);

private:

    const Matrix &A;

};

/** \brief Inversion of precision matrices.
 *
 * Simply returns the reference to the precision matrix that is contained in
 * the PrecisionMatrix wrapper.
 *
 * \tparam Matrix The type of the underlying precision matrix.
 * \param A The precision matrix to invert.
 *
 */
template
<
typename Matrix
>
const Matrix & inv(PrecisionMatrix<Matrix> &A)
{
    return A.A;
}

template
<
typename Matrix
>
const Matrix & inv(const PrecisionMatrix<Matrix> &A)
{
    return A.A;
}

}      // namespace invlib

#endif // PRECISION_MATRIX_H
