/** \file algebra/transformation.h
 *
 * \brief Contains the Transformation class template representing generic
 * coordinate transformation.
 *
 * Also provides two concrete transformation classes, one being the identity
 * transformation and the other the scaling by the reciprocal of the square
 * root of the diagonal elements of a matrix.
 *
 */

#ifndef ALGEBRA_TRANSFORMATION
#define ALGEBRA_TRANSFORMATION

#include <cmath>
#include <utility>
#include <iostream>

namespace invlib
{

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

/**
 * \brief Proxy class for transfomations.
 *
 * This class is used to defer applications of transformations to matrix
 * algebraic expression. Holds a reference to the transformation object and
 * applies either the void apply_vector(Vector&) or the void apply_matrix(Matrix&)
 * member methods, depending if the algebraic expression is evaluated to a
 * vector or to a matrix. The vector type Vector must be provided by the matrix
 * type via .
 *
 * \tparam T1 The type of the algebraic expression.
 * \tparam Transform The tansformation type implementing the transformation.
 */
template
<
typename T1,
typename Transform
>
class Transformation
{

public:

    /*! The basic scalar type. */
    using RealType   = typename decay<T1>::RealType;
    /*! The basic vector type  */
    using VectorType = typename decay<T1>::VectorType;
    /*! The basic matrix type. */
    using MatrixType = typename decay<T1>::MatrixType;
    /*! The result type of the arithmetic expression to transform. */
    using ResultType = typename decay<T1>::ResultType;

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    Transformation(T1 A_, Transform t_);

    Transformation(const Transformation&) = default;
    Transformation(Transformation&&)      = default;

    Transformation & operator=(const Transformation &) = delete;
    Transformation & operator=(Transformation &&)      = delete;

    VectorType multiply(const VectorType &v) const;
    MatrixType multiply(const MatrixType &v) const;
    MatrixType invert() const;
    VectorType solve(const VectorType & v) const;

    // -------------------- //
    // Arithmetic Operators //
    // -------------------- //

    template <typename T2>
    using Sum = MatrixSum<Transformation, T2>;

    template<typename T2>
    Sum<T2> operator+(T2 &&B) const;

    template <typename T2>
    using Difference = MatrixDifference<Transformation, T2>;

    template <typename T2>
    Difference<T2> operator-(T2 &&C) const;

    template <typename T2>
    using Product = MatrixProduct<Transformation, T2>;

    template<typename T2>
    Product<T2> operator*(T2 &&B) const;

    // ------------------- //
    //     Evaluation      //
    // ------------------- //

    operator MatrixType() const;
    operator VectorType() const;

private:

    T1 A;
    Transform t;
};

/**
 * \brief The identity transformation.
 */
class Identity
{
public:

    /**
    * \brief Apply identity.
    */
    template <typename T1>
    constexpr auto apply(T1 &&t)
        -> decltype(std::forward<T1>(t));

};

/**
 * \brief Transform to normalize diagonal of given matrix.
 *
 * When applied to a vector, this transformation scales each component
 * by the reciprocal square root of the absolute value of the corresponding
 * diagonal element  of the matrix @A_. When applied to a matrix, each row and
 * each column are scaled by the reciprocal of the square root of the diagonal
 * elements. When applied to the matrix A itself, the resulting matrix will
 * have only +/- 1.0 on the diagonal.
 *
 * \tparam MatrixType The type of the matrix @A_.
 */
template
<
typename MatrixType
>
class NormalizeDiagonal
{

public:

    NormalizeDiagonal(const MatrixType &A_)
        : A(A_) {}

    template <typename T1>
    void apply_matrix(T1 &) const;

    template <typename T1>
    void apply_vector(T1 &) const;

    template <typename T1>
    using Transform = Transformation<T1, NormalizeDiagonal&>;

    template <typename T1>
    Transform<T1> apply(T1&& A);

private:

    const MatrixType& A;

};

#include "transformation.cpp"

}      // namespace invlib

#endif // ALGEBRA_TRANSFORMATION
