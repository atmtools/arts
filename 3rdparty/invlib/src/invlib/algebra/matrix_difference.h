/** \file algebra/matrix_difference.h
 *
 * \brief Proxy class for computing differences of matrices.
 *
 * Contains the class invlib::MatrixDifference which is a proxy class
 * for computing the difference of marices and vectors.
 *
 */

#ifndef ALGEBRA_MATRIX_DIFFERENCE_H
#define ALGEBRA_MATRIX_DIFFERENCE_H

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

/** \brief Proxy class for computing differences of matrices.
 *
 * The MatrixDifference class template provides a template proxy class
 * for computing the difference of two algebraic expressions.
 *
 * For given type @T1, @T2, for the difference to be evaluatable to
 * type MatrixType, @T1 must be convertible to the @MatrixType type and provide
 * a subtract(const T2&) member function.
 *
 * \tparam T1 The type of the left operand of the difference.
 * \tparam T2 The type of the right operand of the difference.
 * \tparam MatrixType The underlying matrix type.
 *
 */
template
<
typename T1,
typename T2
>
class MatrixDifference
{

public:

    /*! The basic scalar type. */
    using RealType   = typename decay<T2>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<T2>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<T2>::MatrixType;
    /*! The type of the result of the expression */
    using ResultType = typename decay<T2>::ResultType;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    /*! Create arithmetic expression for the difference of the two
     * operands \$fA\f$ and \$fB\$f.
     *
     * \param A The left hand operand of the difference.
     * \param B The right hand operand of the difference.
     * \return Arithmetic expression representing the difference \f$A - B\f$.
     */
    MatrixDifference(T1 A, T2 B);

    /*! Default copy constructor.
     *
     * Proxy objects are lightweight and can be efficiently copied. However,
     * note that there is no caching in place to reuse the results from identical
     * arithmetic operations.
     *
     * \todo Enable caching.
     *
     */
    MatrixDifference(const MatrixDifference &) = default;

    /*! Default move operator. */
    MatrixDifference(MatrixDifference &&) = default;

    /*! Default assignment operator. */
    MatrixDifference & operator=(const MatrixDifference &) = default;

    /*! Default move assignment operator. */
    MatrixDifference & operator=(MatrixDifference &&) = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    /*! Compute product of this difference expression and a vector.
     *
     * Avoids explicit evaluation of the difference by passing the vector down
     * the expression tree and returning the difference of the matrix-vector products
     * of both operands \f$A\f$ and \f$B\f$ and the given vector \f$v\f$. This
     * is done by calling the multiply(const VectorType&) member function of the
     * summands.
     *
     * \param v The vector v to multiply the sum by.
     * \return The matrix-vector product \f$(A - B)v\f$
     */
    VectorType multiply(const VectorType &v) const;

    /*! Compute product of this difference expression and a matrix.
     *
     * Evaluates the difference \f$D = A - B\$f represented by this
     * object and multiplies the result \f$D\f$ by the given matrix
     * \f$C\f$.
     *
     * \param C The matrix \f$C\f$ to multiply the difference by.
     * \return The matrix-matrix product \f$(A - B)C\f$
     */
    MatrixType multiply(const MatrixType &C) const;

    /*! Solve this linear system.
     *
     * Solve the linear system \f$Cx = v\f$ for the linear system
     * corresponding to this matrix difference \f$C = A - B\f$. The solution
     * requires the evaluation of the matrix difference and thus the
     * allocation of an additional matrix. The solution of the linear
     * system is delegated to the fundamental matrix type by calling
     * the solve(const VectorType&) member function of the result of
     * this sum.
     *
     * \param The right hand side vector \f$v\f$ \return The solution
     * vector of the linear system \f$x\f$
     */
    VectorType solve(const VectorType &v) const;

    /*! Invert this matrix difference.
     *
     * Invert the linear system represented by this difference algebraic
     * expression \f$C = A - B\f$. The inversion requires the evaluation of the
     * matrix sum and thus the allocation of an additional matrix. The
     * solution of the linear system is delegated to the fundamental
     * matrix type by calling the invert() member function of the
     * result of this sum.
     *
     * \return The inverse of the matrix represented by this algebraic
     * expression.
     */
    MatrixType invert() const;

    /*! Diagonal of this matrix difference.
     *
     * Compute the diagonal of this matrix difference as the difference of the
     * diagonals of the two operands of the difference. Delegates the computation
     * of the operand diagonals to the diagonal() member functions of the operands.
     *
     * \result A VectorType object representing the diagonal of this matrix difference.
     */
    VectorType diagonal() const;

    /*! Extract row.
     *
     * Returns the ith row of this matrix difference as an object of VectorType. The row
     * is computed as the difference of the corresponding rows of the operands of the
     * difference. Calls the row member functions of the operands.
     *
     * \return A VectorType object representing the ith row of this matrix expression.
     */
    VectorType row(size_t i) const;

    /*! Extract columns.
     *
     * Returns the ith columns of this matrix difference as an object of VectorType. The columns
     * is computed as the difference of the corresponding columns of the operands of the
     * difference. Calls the columns member functions of the operands.
     *
     * \return A VectorType object representing the ith columns of this matrix expression.
     */
    VectorType col(size_t i) const;

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    /*!
     * Proxy type template for a product of this difference expression
     * and another algebraic expression of given type.
     */
    template <typename T3>
    using Product = MatrixProduct<MatrixDifference, T3>;

    /*! Create product arithmetic expression.
     *
     * \tparam T3 The type of the object to multiply this difference with.
     * \param C   The object to multiply this difference with.
     * \return An algebraic expression object representing the product of
     * this difference and the provided argument.
     *
     * \todo Add static asserts.
     */
    template<typename T3>
    Product<T3> operator*(T3 &&C) const;

    /*!
     * Proxy type template for a sum of this difference expression
     * and another algebraic expression of given type.
     */
    template <typename T3>
    using Sum = MatrixSum<MatrixDifference, T3>;

    /*! Create sum arithmetic expression.
     *
     * \tparam T3 The type of the object to add to this difference with.
     * \param C   The object to add to this difference with.
     * \return An algebraic expression object representing the sum of
     * this difference and the provided argument.
     *
     * \todo Add static asserts.
     */
    template <typename T3>
    Sum<T3> operator+(T3 &&C) const;

    /*!
     * Proxy type template for a nested difference of this difference and
     * another algebraic expression of given type.
     */
    template <typename T3>
    using Difference = MatrixDifference<MatrixDifference, T3>;

    /*! Create difference arithmetic expression.
     *
     * \tparam T3 The type of the object to subtract from this difference.
     * \param C   The object to subtract from this difference.
     * \return An algebraic expression object representing the difference of
     * this difference and the provided argument.
     *
     * \todo Add static asserts.
     */
    template <typename T3>
    Difference<T3> operator-(T3 &&C) const;

    operator ResultType() const;

private:

    T1 A;
    T2 B;

};

#include "matrix_difference.cpp"

}

#endif // ALGEBRA_MATRIX_DIFFERENCE_H
