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

int main()
{
    using MatrixTracer = invlib::Tracer<MatrixArchetype<double>>;
    using VectorTracer = invlib::Tracer<VectorArchetype<double>>;
    using MatrixType = invlib::Matrix<MatrixTracer>;
    using VectorType = invlib::Vector<VectorTracer>;

    MatrixTracer::start_tracing();
    MatrixTracer::start_tracing();
    MatrixType M; M.resize(10,10);
    MatrixTracer::stop_tracing("test");
    MatrixTracer::stop_tracing("test");
}
