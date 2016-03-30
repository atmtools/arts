#define BOOST_TEST_MODULE "INVLIB Unit Tests"

#include <boost/test/included/unit_test.hpp>

#include "test_types.h"

#include "algebra/identities.cpp"
#include "algebra/solvers.cpp"
#include "algebra/transformation.cpp"

#include "optimization/exact.cpp"
#include "optimization/test_functions.cpp"

#include "forward_models/linear.cpp"
#include "forward_models/sphere.cpp"
