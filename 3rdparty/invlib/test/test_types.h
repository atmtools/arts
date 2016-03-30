#ifndef TEST_TEST_TYPES_H
#define TEST_TEST_TYPES_H

#include <boost/mpl/list.hpp>
#include <invlib/algebra.h>
#include <invlib/archetypes/matrix_archetype.h>

constexpr double EPS = 1e-8;
constexpr unsigned int ntests = 100;

using Archetype = invlib::Matrix<MatrixArchetype<double>>;
using matrix_types = boost::mpl::list<Archetype>;

#endif // TEST_TEST_TYPES_H
