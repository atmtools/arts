/** \file algebra/matrix_transpose.h
 *
 * \brief Proxy class for computing the transpose of a matrix.
 *
 * Contains the invlib::MatrixTranspose class template with provides a proxy
 * class for transposing an algebraic expression.
 *
 * Also provides the function template transp(), to compute the transpose of
 * an algebraic expression.
 *
 */

#ifndef ALGEBRA_MATRIX_TRANSPOSE_H
#define ALGEBRA_MATRIX_TRANSPOSE_H

namespace invlib
{

/** \brief Proxy class for computing the transpose of a matrix.
 *
 * The MatrixTranspose class template provides a template proxy class
 * for computing the transpose of an algebraic expression.
 *
 * The class assumes that the algebraic expression provides the following
 * member functions:
 *
 * - MatrixType transpose()
 * - MatrixType transpose_multiply(const MatrixType&)
 * - MatrixType transpose_multiply(const VectorType&)
 * - MatrixType transpose_add(const VectorType&)
 *
 * \tparam T1 The type of the algebraic expression to invert.
 * \tparam MatrixType The underlying matrix type.
 *
 */
template
<
typename T1
>
class MatrixTranspose
{

public:

    /*! The basic scalar type. */
    using RealType   = typename decay<T1>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<T1>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<T1>::MatrixType;
    /*! The basic matrix type. */
    using ResultType = MatrixType;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    MatrixTranspose(T1 A);

    MatrixTranspose(const MatrixTranspose &) = default;
    MatrixTranspose(MatrixTranspose &&)      = default;

    MatrixTranspose & operator=(const MatrixTranspose &) = default;
    MatrixTranspose & operator=(MatrixTranspose &&) = default;

    size_t cols() const {return remove_reference_wrapper(A).rows();}
    size_t rows() const {return remove_reference_wrapper(A).cols();}

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    template <typename T2>
    auto multiply(const T2 &t) const -> typename T2::ResultType;

    /*! Multiply this transposed expression by a vector.
     *
     * Evaluates the algebraic expression and uses the
     * tranpose_multiply(const VectorType &) member function of the
     * MatrixType class to multiply this expression by the supplied argument v.
     *
     * \param v The vector to multiply this transposed expression by.
     * \return The result \f$w = A^T2 v\f$
     * \todo Make optional.
     */
    VectorType multiply(const VectorType &v) const;

    /*! Multiply this transposed expression by matrix.
     *
     * Evaluates the algebraic expression and uses the
     * tranpose_multiply(const MatrixType &) member function of the
     * MatrixType class to multiply this expression by the supplied argumnt B.
     *
     * \param B The matrix to multiply this transposed expression by.
     * \return The result \f$C = A^T2 B\f$
     * \todo Make optional.
     */
    MatrixType multiply(const MatrixType &B) const;

    /*! Invert matrix corresponding to this transposed expression.
     *
     * Evaluate the argument to a matrix and uses its invert(const Vector&) member
     * function to invert it.
     *
     * \return The inverse \f$(A^T2)^{-1}\f$ of the matrix \f$A\f$ corresponding to
     * this algebraic expression.
     */
    MatrixType invert() const;

    /*! Solve linear system corresponding to this transposed expression.
     *
     * Evaluate the argument to a matrix and uses its solve(const Vector&) member
     * function to solve the corresponding linear system.
     *
     * \param v The right hand side vector of the linear system.
     * \return The solution \f$x\f$ of the linear system \f$x = A^T2 v\f$
     */
    VectorType solve(const VectorType& v) const;

    /*! Diagonal of this transposed expression.
     *
     * This is an identity operation, so the call is simply forwarded to the underlying
     * expression.
     *
     * \return A VectorType object representing the diagonal of this arithmetic expression.
     */
    VectorType diagonal() const;

    /*! Row i of this transposed expression.
     *
     * The ith row of the transposed expression is the ith column of the non-transposed
     * expression, thus the call is simply fowarded to the column() member function of
     * the underlying expression.
     *
     * \return A VectorType object representing the ith row of this transposed arithmetic
     *  expression.
     */
    VectorType row(size_t i) const;

    /*! Column i of this transposed expression.
     *
     * The ith column of the transposed expression is the ith column of the non-transposed
     * expression, thus the call is simply fowarded to the row() member function of
     * the underlying expression.
     *
     * \return A VectorType object representing the ith column of this transposed
     * arithmetic expression.
     */
    VectorType col(size_t i) const;

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    template <typename T2>
    using Sum = MatrixSum<MatrixTranspose, T2>;

    template<typename T2>
    Sum<T2> operator+(T2 &&B) const;

    template <typename T2>
    using Product = MatrixProduct<MatrixTranspose, T2>;

    template<typename T2>
    Product<T2> operator*(T2 &&A) const;

    operator ResultType() const;

private:

    CopyWrapper<T1> A;

};

/** \brief Transpose of an algebraic expression
 *
 * Returns a MatrixTranspose object instantiated with the type of the given
 * algebraic expression
 *
 * \tparam T2 The type of the algebraic expression to invert.
 * \return The MatrixTranspose proxy object representing the transpose
 * of the given algebraic expression.
 *
 */
template <typename T> void foo();
template <typename T1>
MatrixTranspose<T1> transp(T1 &&A)
{
    return MatrixTranspose<T1>(std::forward<T1>(A));
}

#include "matrix_transpose.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_TRANSPOSE_H
