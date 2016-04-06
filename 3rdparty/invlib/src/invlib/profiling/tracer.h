/** \file profiling/tracer.h
 *
 * \brief Generic class to trace the allocation of matrix and vector objects.
 *
 */

#ifndef PROFILING_TRACER_H
#define PROFILING_TRACER_H

#include <vector>
#include <fstream>
#include "invlib/traits.h"

namespace invlib
{

// ------------------------ //
//    Forward Declarations  //
// ------------------------ //

template
<
typename RealType,
typename Solver
>
class GaussNewton;

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
class LevenbergMarquardt;

// --------------- //
//   Tracer Class  //
// --------------  //

template
<
typename Base,
const char *file_suffix
>
class Tracer : public Base
{
public:

    /*! The basic scalar type. */
    using RealType   = typename Base::RealType;
    /*! The basic vector type  */
    using VectorType = Tracer<typename Base::VectorType, file_suffix>;
    /*! The basic matrix type. */
    using MatrixType = Tracer<typename Base::MatrixType, file_suffix>;
    /*!
     * Result type of an algebraic expression with Matrix as right hand
     * operator
     */
    using ResultType = Tracer<typename Base::ResultType, file_suffix>;

    Tracer();

    Tracer(const Tracer &);
    Tracer(const Base &);

    Tracer(Tracer &&);
    Tracer(Base &&);

    Tracer & operator=(const Tracer &);
    Tracer & operator=(const Base &);

    Tracer & operator=(Tracer &&);
    Tracer & operator=(Base &&);

    ~Tracer();

    template<typename T1>
    void resize(T1 m, T1 n);

    template<typename T1>
    void resize(T1 m);

    static void start_tracing();
    static void stop_tracing(const std::string &filename);

    unsigned int cols() const;

private:

    static unsigned int object_count, size, total_size;
    static std::vector<unsigned int> object_counts, total_sizes;

    template <typename>
    unsigned int cols_base(...) const;

    template <typename T1>
    unsigned int cols_base(decltype(&T1::cols)) const;

};

#include "tracer.cpp"

}      // namespace invlib

#endif // PROFILING_TRACER_H
