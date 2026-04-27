#include <antenna_pattern.h>

#include <cmath>
#include <string_view>

namespace {
Numeric angle_error(Numeric actual, Numeric expected) {
  return std::abs(std::remainder(actual - expected, 360.0));
}

void assert_los(Vector2 actual,
                Vector2 expected,
                Numeric tol,
                std::string_view name) {
  ARTS_USER_ERROR_IF(std::abs(actual[0] - expected[0]) > tol or
                         angle_error(actual[1], expected[1]) > tol,
                     "{} mismatch: actual {} expected {}",
                     name,
                     actual,
                     expected)
}

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

sensor::AntennaPattern pattern(ZenGrid zen_grid, AziGrid azi_grid) {
  sensor::AntennaPattern out;
  out.data.grid<0>() = std::move(zen_grid);
  out.data.grid<1>() = std::move(azi_grid);
  out.data.resize(out.data.grid<0>().size(), out.data.grid<1>().size());

  for (Size izen = 0; izen < out.data.grid<0>().size(); ++izen) {
    for (Size iazi = 0; iazi < out.data.grid<1>().size(); ++iazi) {
      out.data[izen, iazi] = {Numeric(10 * izen + iazi + 1), 0.0, 0.0, 0.0};
    }
  }

  return out;
}

void test_local_wraparound_maps_across_bore() {
  auto ant = pattern(ZenGrid{{1.0}}, AziGrid{{0.0, 180.0}});

  const auto samples       = ant({45.0, 30.0});
  const auto first_weight  = ant.data[0, 0];
  const auto second_weight = ant.data[0, 1];

  ARTS_USER_ERROR_IF(
      samples.size() != 2, "Expected 2 antenna samples, got {}", samples.size())
  ARTS_USER_ERROR_IF(samples[0].first != first_weight,
                     "First Stokvec changed during LOS mapping")
  ARTS_USER_ERROR_IF(samples[1].first != second_weight,
                     "Second Stokvec changed during LOS mapping")

  assert_los(samples[0].second, {46.0, 30.0}, 1e-12, "local [1, 0]");
  assert_los(samples[1].second, {44.0, 30.0}, 1e-12, "local [1, 180]");
}

void test_bore_at_zenith_keeps_defined_local_azimuth() {
  auto ant = pattern(ZenGrid{{1.0}}, AziGrid{{0.0, 90.0, 180.0, 270.0}});

  const auto samples = ant({0.0, 15.0});

  ARTS_USER_ERROR_IF(
      samples.size() != 4, "Expected 4 antenna samples, got {}", samples.size())

  assert_los(samples[0].second, {1.0, 15.0}, 1e-12, "pole local [1, 0]");
  assert_los(samples[1].second, {1.0, 105.0}, 1e-12, "pole local [1, 90]");
  assert_los(samples[2].second, {1.0, 195.0}, 1e-12, "pole local [1, 180]");
  assert_los(samples[3].second, {1.0, 285.0}, 1e-12, "pole local [1, 270]");
}

void test_pencil_beam_defaults_and_custom_weight() {
  const sensor::PencilBeamAntenna def{};

  ARTS_USER_ERROR_IF(
      def.data.grid<0>().size() != 1 or def.data.grid<1>().size() != 1,
      "Pencil beam must be 1x1")
  ARTS_USER_ERROR_IF(
      def.data.grid<0>()[0] != 0.0 or def.data.grid<1>()[0] != 0.0,
      "Pencil beam grid must be centered at [0, 0]")
  assert_stokvec(
      def.data[0, 0], {1.0, 0.0, 0.0, 0.0}, 0.0, "default pencil beam weight");

  const auto def_samples = def({45.0, 30.0});
  ARTS_USER_ERROR_IF(def_samples.size() != 1,
                     "Default pencil beam should return one sample")
  assert_stokvec(def_samples[0].first,
                 {1.0, 0.0, 0.0, 0.0},
                 0.0,
                 "default pencil beam mapped weight");
  assert_los(
      def_samples[0].second, {45.0, 30.0}, 1e-12, "default pencil beam LOS");

  const Stokvec custom_weight{0.5, 0.25, 0.0, 0.0};
  const sensor::PencilBeamAntenna custom{custom_weight};
  assert_stokvec(
      custom.data[0, 0], custom_weight, 0.0, "custom pencil beam weight");
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
  test_local_wraparound_maps_across_bore();
  test_bore_at_zenith_keeps_defined_local_azimuth();
  test_pencil_beam_defaults_and_custom_weight();
  test_gaussian_initialization_uses_antenna_frame_offsets();
  return 0;
}