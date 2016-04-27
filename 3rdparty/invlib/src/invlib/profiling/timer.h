/** \file profiling/timer.h
 *
 * \brief Generic class to time arithmetic operations.
 *
 */

#ifndef PROFILING_TIMER_H
#define PROFILING_TIMER_H

#include <vector>
#include <fstream>
#include <chrono>

using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include "invlib/traits.h"

namespace invlib
{

duration<double> multiply_mv_time;
duration<double> multiply_mm_time;
duration<double> solve_time;
duration<double> invert_time;

/*! Generic timer class that wraps around a generic vector or matrix type and
 *  records the time spent in the function for Matrix-Vector multiplication,
 *  matrix-matrix multiplication, solution of linear systems and inversion of
 *  matrices in the global variablese multiply_mv_time, multiply_mm_time,
 *  solve_time, invert_time. Note that thos variables don't distinguis between
 *  underlying types. Thus if different types provide for example matrix-vector
 *  multiplication and are timed, the total time spent in those operations
 *  will be accumulated into the same global variable.
 *
 * \tparam Base The type of the type to be timed.
 */
template
<
typename Base
>
class Timer : public Base
{

public:

    using BaseType   = Base;
    /*! The basic scalar type. */
    using RealType   = typename Base::RealType;
    /*! The basic vector type  */
    using VectorType = Timer<typename Base::VectorType>;
    /*! The basic matrix type. */
    using MatrixType = Timer<typename Base::MatrixType>;
    /*!
     * Result type of an algebraic expression with Matrix as right hand
     * operator
     */
    using ResultType = Timer<typename Base::ResultType>;

    // -------------------------- //
    //  Construction & Assignment //
    // -------------------------- //

    Timer() = default;

    template
    <
    typename T,
    typename = is_constructible<Base, T>
    >
    Timer(T &&);

    template
    <
    typename T,
    typename = is_assignable<Base, T>
    >
    Timer & operator=(T &&);

    template
    <
    typename T,
    typename = is_copy_assignable<T>
    >
    Timer & operator=(const T &);

    // ------------------------ //
    //   Arithmetic Operations  //
    // ------------------------ //

    VectorType multiply(const VectorType &) const;
    VectorType transpose_multiply(const VectorType &) const;

    MatrixType multiply(const MatrixType &) const;
    MatrixType transpose_multiply(const MatrixType &) const;

    VectorType solve(const VectorType &) const;
    MatrixType invert() const;
};

#include "timer.cpp"

}      // namespace invlib

#endif // PROFILING_TIMER_H
