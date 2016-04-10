#ifndef TEST_TEST_TYPES_H
#define TEST_TEST_TYPES_H

#include <boost/mpl/list.hpp>
#include <invlib/algebra.h>
#include <invlib/interfaces/arts_wrapper.h>
#include <invlib/archetypes/matrix_archetype.h>

constexpr double EPS = 1e-8;
constexpr unsigned int ntests = 100;

using Archetype = invlib::Matrix<MatrixArchetype<double>>;
using Arts = invlib::Matrix<ArtsMatrix>;
using matrix_types = boost::mpl::list<Arts>;

#endif // TEST_TEST_TYPES_H
