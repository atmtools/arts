#include <iostream>

#include <scattering/maths.h>

namespace scattering {

using math::Tensor;

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

/** Generator to fill tensor with element index.
 *
 * @tparam Scalar the scalar type of the tensor.
 * @tparam rank The rank of the tensor to fill.
 */
template <typename Scalar, long int rank>
struct IndexGenerator {
    IndexGenerator(int dimension) : dimension_(dimension) {}
    Scalar operator()(const std::array<long int, rank> &coordinates) const {
        return static_cast<Scalar>(coordinates[dimension_]);
    }
    int dimension_;
};


/** Create a tensor and fill with indices.
 *
 * @tparam rank The rank of the tensor to create.
 * @param sizes An array specifying the dimensions of the tensor.
 * @param fill_along An integer specifying the dimension with whose indices
 * to fill the tensor.
 * @return The tensor of the requested rank and size with each element
 * containing the corresponding index  along dimension 'fill_along'
 *
 */
template<int rank>
Tensor<double, rank> make_tensor(
    std::array<long int, rank> sizes,
    int fill_along
    ) {
    Tensor<double, rank> tensor{sizes};
    return tensor.generate(IndexGenerator<double, rank>(fill_along));
}

/** Create a tensor of random size and fill with indices.
 *
 * @tparam rank The rank of the tensor to create.
 * @param sizes An array specifying the dimensions of the tensor.
 * @param fill_along An integer specifying the dimension with whose indices
 * to fill the tensor.
 * @return The tensor of the requested rank and size with each element
 * containing the corresponding index  along dimension 'fill_along'
 *
 */
template<int rank>
Tensor<double, rank> make_tensor(int fill_along) {
    std::array<long int, rank> dimensions{};
    for (int i = 0; i < rank; ++i) {
        dimensions[i] = 1 + rand() % 10;
    }
    return make_tensor<rank>(dimensions, fill_along);
}

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////

}

bool test_tensor_indexing(int n_tests) {
    for (int trial = 0; trial < n_tests; ++trial) {
      auto tensor = scattering::make_tensor<4>(0);

      // Test indexing along first dimension.
      for (int i = 0; i < tensor.dimensions()[0]; ++i) {
          auto part = scattering::math::tensor_index<1>(tensor, {i});
          double diff = std::abs(part(0, 0, 0) - static_cast<double>(i));
          if (!scattering::math::small(diff)) {
              return false;
          }
      }

      // Test indexing along second dimension.
      tensor = scattering::make_tensor<4>(1);
      for (int i = 0; i < tensor.dimensions()[1]; ++i) {
          auto part = scattering::math::tensor_index<2>(tensor, {0, i});
          double diff = std::abs(part(0, 0) - static_cast<double>(i));
          if (!scattering::math::small(diff)) {
              return false;
          }
      }

      // Test indexing along third dimension.
      tensor = scattering::make_tensor<4>(2);
      for (int i = 0; i < tensor.dimensions()[2]; ++i) {
          auto part = scattering::math::tensor_index<3>(tensor, {0, 0, i});
          double diff = std::abs(part(0) - static_cast<double>(i));
          if (!scattering::math::small(diff)) {
              return false;
          }
      }

    }
    return true;
}


// Test creation of matrix maps from tensors.
bool test_to_map() {
    auto tensor = scattering::make_tensor<2>(0);
    auto matrix_map = scattering::math::to_matrix_map(tensor);
    for (int i = 0; i < tensor.dimensions()[0]; ++i) {
        double diff = std::abs(matrix_map(i, 0) - static_cast<double>(i));
        if (!scattering::math::small(diff)) {
            return false;
        }
    }
    matrix_map(0, 0) = 99.0;
    if (matrix_map(0, 0) != tensor(0, 0)) {
        return false;
    }

    auto vector_map = scattering::math::to_vector_map(
        scattering::math::tensor_index<1>(tensor, {0})
        );
    if (vector_map[0] != tensor(0, 0)) {
        return false;
    }
    return true;
}

// Test extraction of submatrices from tensors.
bool test_submatrix() {
    auto tensor = scattering::make_tensor<4>(1);

    auto matrix_map = scattering::math::get_submatrix<1, 3>(tensor, {0, 0});
    for (int i = 0; i < matrix_map.rows(); ++i) {
        double diff = std::abs(matrix_map(i, 0) - static_cast<double>(i));
        if (!scattering::math::small(diff)) {
            return false;
        }
    }
    matrix_map(1, 2) = 99.0;
    if (tensor(0, 1, 0, 2) != 99.0) {
        return false;
    }
    return true;
}


// Test extraction of submatrices from tensors.
bool test_subvector() {
    auto tensor = scattering::make_tensor<4>(1);

    auto vector_map = scattering::math::get_subvector<2>(tensor, {0, 0, 0});
    for (int i = 0; i < vector_map.rows(); ++i) {
        double diff = std::abs(vector_map(i, 0) - static_cast<double>(i));
        if (!scattering::math::small(diff)) {
            return false;
        }
    }
    vector_map(2) = 99.0;
    if (tensor(0, 0, 2, 0) != 99.0) {
        return false;
    }
    return true;
}

/// Test dimension iterator.
bool test_dimension_iterator() {
    auto tensor = scattering::make_tensor<4>(1);
    auto ctr = scattering::math::DimensionCounter<4>(tensor.dimensions());

    for (;ctr; ++ctr) {
        double diff = std::abs(
            tensor(ctr.coordinates) -
            static_cast<double>(ctr.coordinates[1])
            );
        if (!scattering::math::small(diff)) {
            return false;
        }
    }
    return true;
}

int main(int /*nargs*/, char **/*argv*/) {

    auto passed = test_tensor_indexing(10);
    std::cout << "test_tensor_indexing: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }

    passed = test_to_map();
    std::cout << "test_tensor_map: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }

    passed = test_submatrix();
    std::cout << "test_submatrix: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }

    passed = test_subvector();
    std::cout << "test_subvector: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }

    passed = test_dimension_iterator();
    std::cout << "test_dimension_iterator: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }

    return 1;
}
