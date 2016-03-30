/** \file algebra/matrix_sum.h
 *
 * \brief Proxy class for computing the sum of two algebraic expressions.
 *
 */

#ifndef ALGEBRA_SUM_H
#define ALGEBRA_SUM_H

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
class MatrixDifference;

/** \brief Proxy class for computing the sum of two algebraic expressions.
 *
 * The MatrixSum class template provides a proxy class template for computing the
 * sum of two algebraic expressions. The class expects the left hand side
 * operand to be convertible to either a vector or a matrix, which then must
 * provide a member function accumulate, which can be called with the right
 * hand operand as only argument.
 *
 * \tparam T1 The type of the left hand side operand
 * \tparam T2 the type of the right hand side operand
 * \tparam MatrixType The underlying matrix type used.
 *
 */
template
<
typename T1,
typename T2
>
class MatrixSum
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

    /*!
     * Create a proxy object for the two given algebraic expressions.
     */
    MatrixSum(T1 Op1, T2 Op2);

    /*! Default copy constructor.
     *
     * Proxy objects are lightweight and can be efficiently copied. However,
     * note that there is no caching in place to reuse the results from identical
     * arithmetic operations.
     *
     * \todo Enable caching.
     *
     */
    MatrixSum(const MatrixSum &) = default;

    /*! Default move constructor. */
    MatrixSum(MatrixSum &&) = default;

    /*! Default assignment operator. */
    MatrixSum & operator=(const MatrixSum &) = default;
    MatrixSum & operator=(MatrixSum &&)      = delete;

    /*! Default destructor. */
    ~MatrixSum() = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    /*! Multiply sum expression by a vector.
     *
     * Avoids explicit evaluation of the sum by passing the vector down
     * the expression tree and returning the sum of the matrix-vector products
     * of both summands \f$A\f$ and \f$B\f$ and the given vector \f$v\f$. This
     * is done by calling the multiply(const VectorType&) member function of the
     * summands.
     *
     * \param v The vector v to multiply the sum by.
     * \return The matrix-vector product \f$(A + B)v\f$
     */
    VectorType multiply(const VectorType &v) const;

    /*! Multiply sum expression by a matrix.
     *
     * Evaluates the sum \f$D = A + B\$f represented by this object and multiplies
     * the result \f$D\f$ by the given matrix \f$C\f$.
     *
     * \param C The matrix \f$C\f$ to multiply the sum by.
     * \return The matrix-matrix product \f$(A + B)C\f$
     */
    MatrixType multiply(const MatrixType &C) const;

    /*! Solve this linear system.
     *
     * Solve the linear system \f$Cx = v\f$ for the linear system
     * corresponding to this matrix sum \f$C = A + B\f$. The solution
     * requires the evaluation of the matrix sum and thus the
     * allocation of an additional matrix. The solution of the linear
     * system is delegated to the fundamental matrix type by calling
     * the solve(const VectorType&) member function of the result of
     * this sum.
     *
     * \param The right hand side vector \f$v\f$
     * \return The solution vector of the linear system \f$x\f$
     */
    VectorType solve(const VectorType &v) const;

    /*! Invert this matrix sum.
     *
     * Inver the sum \f$C = A + B\f$ represented by this algebraic
     * expression.  The inversion requires the evaluation of the
     * matrix sum and thus the allocation of an additional matrix. The
     * solution of the linear system is delegated to the fundamental
     * matrix type by calling the invert() member function of the
     * result of this sum.
     *
     * \return The inverse of the matrix represented by this algebraic
     * expression.
     */
    MatrixType invert() const;

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    /*!
     * Proxy type template for a product of this sum expression and another algebraic
     * expression of given type.
     */
    template <typename T3>
    using Product = MatrixProduct<MatrixSum, T3>;

    /*! Create product arithmetic expression.
     *
     * \tparam T3 The type of the object to multiply this sum with.
     * \param C   The object to multiply this sum with.
     * \return An algebraic expression object representing the product of
     * this sum and the provided argument.
     *
     * \todo Add static asserts.
     */
    template<typename T3>
    Product<T3> operator*(T3 &&C) const;

    /*!
     * Proxy type template for a nested sum of this sum and another algebraic
     * expression of given type.
     */
    template <typename T3>
    using Sum = MatrixSum<MatrixSum, T3>;

    /*! Create sum arithmetic expression.
     *
     * \tparam T3 The type of the object to add to this sum.
     * \param C   The object to add to this sum with.
     * \return An algebraic expression object representing the sum of
     * this sum and the provided argument.
     *
     * \todo Add static asserts.
     */
    template <typename T3>
    Sum<T3> operator+(T3 &&C) const;

    /*!
     * Proxy type template for the difference sum of this sum and another algebraic
     * expression of given type.
     */
    template <typename T3>
    using Difference = MatrixDifference<MatrixSum, T3>;

    /*! Create difference arithmetic expression.
     *
     * \tparam T3 The type of the object to subtract from this sum with.
     * \param   C The object to subtract from this sum.
     * \return An algebraic expression object representing the difference of
     * this sum and the provided argument.
     *
     * \todo Add static asserts.
     */
    template <typename T3>
    Difference<T3> operator-(T3 &&C) const;

    // -----------------//
    //     Evaluation   //
    // ---------------- //

    /*! Evaluate sum.
     *
     * Evaluate the sum \f$C = A + B\f$ represented by this algebraic
     * expression. Evaluates the left hand operand \f$A\f$ to the
     * ResultType and the uses the accumulate member function to add
     * the right hand operand \f$B\f$ to the result.
     *
     * \return The evaluated sum represented by this algebraic expression.
     */
    operator ResultType() const;

private:

    T1 A;
    T2 B;

};

#include "matrix_sum.cpp"

}

#endif // ALGEBRA_SUM_H
