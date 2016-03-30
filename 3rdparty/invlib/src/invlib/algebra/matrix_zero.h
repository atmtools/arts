#ifndef ALGEBRA_MATRIX_ZERO_H
#define ALGEBRA_MATRIX_ZERO_H

#include <utility>
#include "matrix_inverse.h"

/** \file algebra/matrix_zero.h
 *
 * Contains MatrixZero class representing a generic zero matrix.
 * Also provides and overload of the inv() function for zero matrices.
 *
 */
namespace invlib
{

 /*! Represents the zero matrix, i.e. the identity matrix of matrix addition.
 *
 * \tparam Matrix The underlying matrix type.
 *
 */
class MatrixZero
{

public:

    // ------------------------------- //
    //  Constructors and Destructors   //
    // ------------------------------- //

    MatrixZero() = default;

    MatrixZero(const MatrixZero&) = default;
    MatrixZero(MatrixZero&&)      = default;

    MatrixZero& operator= (const MatrixZero&) = default;
    MatrixZero& operator= (MatrixZero&&)      = default;

    // --------------------- //
    //   Nested Evaluation   //
    // --------------------- //

    template <typename T1>
    T1 && multiply(T1 &&A) const;

    MatrixZero invert() const;

    template <typename T1>
        T1 && solve(T1 && v) const;

    // --------------------- //
    // Arithmetic Operators  //
    // --------------------- //

    template <typename T1>
        using Sum = T1;

    template<typename T1>
    T1 && operator+(T1 &&B) const;

    template <typename T1>
    using Product = MatrixProduct<MatrixZero, T1>;

    template <typename T1>
    Product<T1> operator*(const T1 &A) const;

};

 /*! Inverse of inverse matrix, here defined as zero matrix, so that
  * terms involving the inverse of a zero matrix are simply ignored.
  *
  * \return A zero matrix object.
  */
MatrixZero inv(const MatrixZero &A);

#include "matrix_zero.cpp"

}      // namespace invlib

#endif // ALGEBRA_MATRIX_ZERO_H
