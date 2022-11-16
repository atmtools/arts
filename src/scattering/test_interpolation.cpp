#include <iostream>

#include <scattering/maths.h>
#include <scattering/interpolation.h>

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

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////

template <math::Index  Rank, math::Index Axis>
/// Tests regridding along a specific axis.
bool test_interpolation() {

    bool passed = true;

    math::Vector<double> grid = math::Vector<double>::LinSpaced(10, 0, 9);
    std::array<math::Vector<double>, 1> old_grid{grid};
    std::array<math::Vector<double>, 1> new_grid{grid};
    math::Vector<double> interp_coords(3);
    interp_coords[0] = -1.0;
    interp_coords[1] = 4.5;
    interp_coords[2] = 10.0;
    new_grid[0] = interp_coords;

    auto regridder = RegularRegridder<double, Axis>{
        old_grid,
        new_grid,
            true
    };


    std::array<long int, Rank> sizes{};
    sizes.fill(100);

    auto input_tensor = make_tensor<Rank>(
        sizes,
        static_cast<int>(Axis)
        );

    auto output_tensor = regridder.regrid(input_tensor);

    auto coords = std::array<math::Index, Rank>{};
    for (int i = 0; i < 3; ++i) {
        coords[Axis] = i;
        auto ref = interp_coords[i];
        auto val = output_tensor(coords);
        if (!math::small(std::abs(ref - val))) {
                passed = false;
            }
    }
    return passed;
}
}

int main(int /*nargs*/, char **/*argv*/) {
    bool passed = true;
    passed = scattering::test_interpolation<4, 0>();
    if (!passed) {
        std::cout << "Failed: test_interpolation<4, 0>" << std::endl;
    }
    passed = scattering::test_interpolation<4, 1>();
    if (!passed) {
        std::cout << "Failed: test_interpolation<4, 0>" << std::endl;
    }
    passed = scattering::test_interpolation<4, 2>();
    if (!passed) {
        std::cout << "Failed: test_interpolation<4, 0>" << std::endl;
    }
    passed = scattering::test_interpolation<4, 3>();
    if (!passed) {
        std::cout << "Failed: test_interpolation<4, 0>" << std::endl;
    }
}
