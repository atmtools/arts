#include "utils.h"

namespace scattering {

GridPos find_interp_weights(StridedConstVectorView grid,
                            Numeric x_new) {
    // Ensure the input vector is sorted
    if (grid.empty() || grid.size() == 1) {
        throw std::invalid_argument("Input vector must have at least two sorted elements.");
    }

    // Handle boundary cases
    if (x_new <= grid.front()) {
        return GridPos(0, {{1.0, 0.0}}); // Entire weight to the first point
    }
    if (x_new >= grid.back()) {
        return GridPos(grid.size() - 1, {{0.0, 1.0}}); // Entire weight to the last point
    }

    // Find the interval using std::lower_bound
    auto it = std::lower_bound(grid.begin(), grid.end(), x_new);
    size_t index = std::distance(grid.begin(), it);

    if (*it == x_new) {
        return GridPos(index, {{1.0, 0.0}});
    }

    Numeric x0 = grid[index - 1];
    Numeric x1 = grid[index];

    Numeric t = (x_new - x0) / (x1 - x0);
    return GridPos(index, {{1.0 - t, t}});

}

}
