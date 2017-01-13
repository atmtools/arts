/** \file algebra/preconditioners.h
 *
 * \brief Preconditioners for the preconditioned conjugate gradient method.
 *
 */

#ifndef ALGEBRA_PRECONDITIONERS_H
#define ALGEBRA_PRECONDITIONERS_H

#include <utility>

#include "invlib/algebra/matrix.h"
#include "invlib/algebra/vector.h"

namespace invlib
{

// ------------------------- //
//  Jacobian Preconditioner  //
// ------------------------- //
/**
 * \brief Jacobian Preconditioner for the OEM.
 *
 * Jacobian preconditioner for the linear system
 * \f[
 *   \left(\mathbf{K}^T\mathbf{S}_\epsilon\mathbf{K} + \mathbf{S}_a\right) \mathbf{x}
 *    = \mathbf{b}
 * \f]
 * This class acts as a functor applying the action of the inverse of the preconditioner
 * \f$\mathbf{M} = \mathbf{C}^T\mathbf{C}\f$ with
 * \f[
 *  \mathbf{C} = diag
 *  \left(\mathbf{K}^T\mathbf{S}_\epsilon\mathbf{K} + \mathbf{S}_a\right)^{-1/2},
 * \f]
 * , i.e. scales the vector by the element-wise inverse of the diagonal of the matrix
 * defining the linear system (c.f. Chapter 5 in \cite nocedal).
 */
template
<
    typename VectorType
>
class JacobianPreconditioner
{
public:

    /*! Construct a Jacobian preconditioner from a given algebraic expression.
     *
     * Computes the diagonal of the expression. All underlying types must thus
     * provide the necessary diagonal(), respectively row(...) and col(...)
     * functions.
     */
    template <typename T1>
    JacobianPreconditioner(const T1 &);

    /** Apply the action of the inverse preconditioner to the given vector. */
    VectorType operator()(const VectorType &) const;

private:

    VectorType diag;

};

#include "preconditioners.cpp"

}       // namespace invlib

# endif // ALGEBRA_MATRIX_H
