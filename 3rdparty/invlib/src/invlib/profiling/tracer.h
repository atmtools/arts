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

// ----------------------- //
//  TypeName Helper Class  //
// ----------------------  //

template <typename T>
struct TypeName
{
    static constexpr auto name = "void";
};

template <typename T>
struct TypeName<Vector<T>>
{
    static constexpr auto name = "Vector";
};

template <typename T>
struct TypeName<Matrix<T>>
{
    static constexpr auto name = "Matrix";
};

// --------------- //
//   Tracer Class  //
// --------------  //

template
<
typename Base
>
class Tracer : public Base
{
public:

    Tracer() = default;

    template<typename T>
    Tracer(T &&);

    template<typename T>
    Tracer & operator=(T &&);

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
