#include <antenna_pattern.h>

#include <cmath>
#include <string_view>

namespace {
void assert_stokvec(Stokvec actual,
                    Stokvec expected,
                    Numeric tol,
                    std::string_view name) {
  for (Size i = 0; i < 4; ++i) {
    ARTS_USER_ERROR_IF(std::abs(actual[i] - expected[i]) > tol,
                       "{} mismatch at {}: actual {} expected {}",
                       name,
                       i,
                       actual,
                       expected)
  }
}

void test_gaussian_initialization_uses_antenna_frame_offsets() {
  const Stokvec peak_weight{2.0, 1.0, 0.0, 0.0};
  const sensor::GaussianAntenna ant{
      ZenGrid{{0.0, 1.0}}, AziGrid{{0.0, 180.0}}, 1.0, 2.0, peak_weight};

  assert_stokvec(ant.data[0, 0], peak_weight, 1e-12, "gaussian peak weight");

  const Numeric scale    = std::exp(-0.5);
  const Stokvec expected = scale * peak_weight;
  assert_stokvec(ant.data[1, 0], expected, 1e-12, "gaussian local [1, 0]");
  assert_stokvec(ant.data[1, 1], expected, 1e-12, "gaussian local [1, 180]");

  const sensor::GaussianAntenna def{ZenGrid{{0.0}}, AziGrid{{0.0}}, 1.0, 1.0};
  assert_stokvec(def.data[0, 0],
                 {1.0, 0.0, 0.0, 0.0},
                 1e-12,
                 "default gaussian peak weight");
}
}  // namespace

int main() {
  test_gaussian_initialization_uses_antenna_frame_offsets();
  return 0;
}