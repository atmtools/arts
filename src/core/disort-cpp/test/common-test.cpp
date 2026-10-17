#include <common.h>

#include <stdexcept>

namespace {

void expect_invalid_omega(const Numeric omega) {
  try {
    disort_common::check_layer_input_values(AscendingGrid{1.0}, Vector{omega}, Vector{0.1}, 0.5, 0.0, true);
  } catch (const std::exception&) { return; }

  throw std::runtime_error("DISORT accepted a single-scattering albedo outside [0, 1]");
}

}  // namespace

int main() {
  disort_common::check_layer_input_values(AscendingGrid{1.0}, Vector{0.0}, Vector{0.1}, 0.5, 0.0, true);
  disort_common::check_layer_input_values(AscendingGrid{1.0}, Vector{1.0}, Vector{0.1}, 0.5, 0.0, true);
  expect_invalid_omega(-1.0e-12);
  expect_invalid_omega(1.0 + 1.0e-12);
}
