/** \file profiling/tracer.h
 *
 * \brief Generate trace of allocation of Matrix and Vector object during
 * map computation.
 *
 */

#include "invlib/algebra.h"
#include "invlib/profiling/tracer.h"
#include "invlib/algebra/solvers.h"
#include "invlib/map.h"
#include "invlib/optimization.h"

#include "invlib/archetypes/matrix_archetype.h"
#include "invlib/archetypes/vector_archetype.h"

#include <utility>

constexpr char vector_suffix[] = "vector";
constexpr char matrix_suffix[] = "matrix";

int main()
{
    using VectorTracer = invlib::Tracer<VectorArchetype<double>, vector_suffix>;
    using MatrixTracer = invlib::Tracer<MatrixArchetype<double>, matrix_suffix>;
    using MatrixType = invlib::Matrix<MatrixTracer>;
    using VectorType = invlib::Vector<VectorTracer>;

    MatrixTracer::start_tracing();
    VectorTracer::start_tracing();

    MatrixType A, B, C, D; A.resize(100,100);

    B = A * A; B.resize(0,0);

    MatrixTracer::stop_tracing("test");
    VectorTracer::stop_tracing("test");
}
