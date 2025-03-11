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
            return GridPos(0, {{0.0, 1.0}}); // Entire weight to the first point
        }
        if (x_new >= grid.back()) {
            return GridPos(grid.size() - 2, {{1.0, 0.0}}); // Entire weight to the last point
        }

        // Find the interval using std::lower_bound
        auto it = std::lower_bound(grid.begin(), grid.end(), x_new);
        size_t index = std::distance(grid.begin(), it);

        if (*it == x_new) {
            return GridPos(index, {{0.0, 1.0}});
        }

        Numeric x0 = grid[index - 1];
        Numeric x1 = grid[index];

        Numeric t = (x_new - x0) / (x1 - x0);
        return GridPos(index - 1, {{t, 1.0 - t}});

    }

    /**! Find bin index for a given value.
     *
     * Returns the 0-based bin index of a given value.
     *
     * @param boundaries: A vector containing the n + 1 bin-boundaries of the n bins in increasing order.
     * @param value: The value to sort into the bins.
     *
     * @return  The index of the bin or:
     *     - -1 if value is smaller than any of the bin boundaries
     *     - n + 1 if value is larger than any of the bin boundaries.
     */
    Index digitize(const Vector& boundaries, Numeric value) {

        // Handle boundary cases
        if (value <= boundaries.front()) {
            return -1;
        }
        if (value >= boundaries.back()) {
            return boundaries.size();
        }

        // Find the interval using std::lower_bound
        auto it = std::lower_bound(boundaries.begin(), boundaries.end(), value);
        Index index = std::distance(boundaries.begin(), it);
        return index - 1;
    }

}
