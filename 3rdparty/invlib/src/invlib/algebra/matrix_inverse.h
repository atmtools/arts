/** \file algebra/matrix_inverse.h
 *
 * \brief Proxy class for computing matrix inverses.
 *
 * Contains the class invlib::MatrixInverse which is a proxy class
 * for computing the inverse of a matrix or solving the corresponding
 * system of equations.
 *
 * Also contains the generic inv() function which creates a MatrixInverse
 * object from a given algebraic expression.
 *
 */
#ifndef ALGEBRA_MATRIX_INVERSE
#define ALGEBRA_MATRIX_INVERSE

#include "invlib/traits.h"

namespace invlib
{

// -------------------- //
// Forward Declarations //
// -------------------- //

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

// --------------------- //
//  MatrixInverse Class  //
// --------------------- //

/** \brief Proxy class for computing matrix inverses.
 *
 * The MatrixDifference class template provides a template proxy class
 * for computing the inverses of matrices or solving the corresponding
 * linear system of equations.
 *
 * The key purpose of the class is to delay the computation of the inverse
 * until it is either multiplied by a vector or used in another algebraic
 * operation (including multiplication by another matrix). If the MatrixInverse
 * object is multiplied by a vector only the corresponding linear system
 * has to be solved. All other operations require the inversion of the matrix.
 *
 * A MatrixInverse object can be used in another operation or converted to
 * an object of type MatrixType. To this end, the T1 type is required to provide
 * the member functions
 *
 * - VectorType solve(const typename VectorType&)
 * - MatrixType invert()
 *
 * for solving the associated linear system of equations and inverting the
 * corresponding matrix, respectively.
 *
 * \tparam T1 The type of the algebraic expression to invert.
 * \tparam MatrixType The underlying matrix type.
 *
 */
template
<
typename T1
>
class MatrixInverse
{
public:

    /*! The basic scalar type. */
    using RealType   = typename decay<T1>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<T1>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<T1>::MatrixType;
    /*! The type of the result of the expression */
    using ResultType = typename decay<T1>::ResultType;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Create an arithmetic expression representing the inverse of another
     * arithmetic expression.
     *
     * \param A The arithmetic expression to invert.
     *
     */
    MatrixInverse(const T1 &A);

    MatrixInverse(const MatrixInverse &) = default;
    MatrixInverse(MatrixInverse &&)      = default;
    MatrixInverse & operator=(const MatrixInverse &) = default;
    MatrixInverse & operator=(MatrixInverse &&)      = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    /* Compute product of this matrix inverse and given vector.
     *
     * Does not evaluate this matrix inverse but instead solves the
     * corresponding linear system which is less costly.
     *
     * \param v The vector to multiply the inverse by.
     * \return The vector \f$w = A^{-1} v\f$.
     */
    VectorType multiply(const VectorType &v) const;

    /* Compute product of this matrix inverse and a given matrix.
     *
     * Evaluates matrix inverse and multiplies result by the given matrix.
     *
     * \param B The matrix to multiply this inverse by.
     * \return The matrix \f$C = A^{-1} B\f$.
     *
     * \todo Could be optimized by combining inversion and product.
     */
    template <typename T2, typename U = disable_if<has_solve<MatrixType, T2>>>
    typename T2::ResultType multiply(const T2 &B) const;

    /*! Evaluates the algebraic expression and multiplies it
     * from the right by the given vector.
     */
    VectorType solve(const VectorType &v) const;

    /*! The inverse of the inverse of an expression is just the result of the
     * expression itself.
     */
    MatrixType invert() const;

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    /*!
     * Proxy type template for a product of this inverse expression and
     * another algebraic expression of given type.
     */
    template<typename T2>
    using Product = MatrixProduct<MatrixInverse, T2>;

    /*! Create algebraic expression for the product of this inverse and another
     * given algebraic expression.
     *
     * \tparam The type of the algebraic expression to multiply this matrix
     * by.
     * \param B The algebraic expression to multiply this matrix inverse by.
     * \return An algebraic expression representing the product of this matrix
     * inverse and the algebraic expression B.
     */
    template <typename T2>
    Product<T2> operator*(T2 &&B) const;

    /*!
     * Proxy type template for the sum of this inverse expression and
     * another algebraic expression of given type.
     */
    template<typename T2>
    using Sum = MatrixSum<MatrixInverse, T2>;

    /*! Create algebraic expression for the sum of this inverse and another
     * given algebraic expression.
     *
     * \tparam The type of the algebraic expression add to this matrix inverse.
     *
     * \param B The algebraic expression to add to this matrix inverse.
     * \return An algebraic expression representing the sum of this matrix
     * inverse and the algebraic expression B.
     */
    template <typename T2>
    Sum<T2> operator+(T2 &&B) const;

    /*!
     * Proxy type template for the difference of this inverse expression and
     * another algebraic expression of given type.
     */
    template<typename T2>
    using Difference = MatrixDifference<MatrixInverse, T2>;

    /*! Create algebraic expression for the difference of this inverse and another
     * given algebraic expression.
     *
     * \tparam The type of the algebraic expression subtact from this matrix inverse.
     *
     * \param B The algebraic expression to subtract from this matrix inverse.
     * \return An algebraic expression representing the difference of this matrix
     * inverse and the algebraic expression B.
     */
    template <typename T2>
    Difference<T2> operator-(T2 &&B) const;

    // -----------------//
    //     Evaluation   //
    // ---------------- //

    /*! Evaluate matrix inverse.
     *
     * Evaluates the matrix inverse by first converting the expression
     * to MatrixType and then calling the member function invert() of
     * the fundamental matrix type.
     */
    operator ResultType() const;

private:

    T1 A;

};

template<typename Matrix>
class PrecisionMatrix;
/** \brief Inverse of an algebraic expression.
 *
 * Creates a proxy object of type MatrixInverse<T2> will evaluate to
 * either inversion of the corresponding matrix or solution of the
 * corresponding linear system.
 *
 * If the resulting MatrixInverse object is multiplied from the right
 * with a vector the corresponding linear system will only be
 * solved. All other operations will result in the inversion of the
 * corresponding system.
 *
 * \tparam T2 The type of the algebraic expression.
 *
 * \param A The algebraic expression to be inverted.
 *
 * \return The MatrixInverse object representing the inverted algebraic
 * expression.
 *
 */
template <typename T1, typename T2 = enable_if<is_invertible<T1>>>
MatrixInverse<T1> inv(T1 &&A);

#include "matrix_inverse.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_INVERSE

