#include <antenna_pattern.h>
#include <arts_conversions.h>
#include <planet_data.h>

#include <cmath>
#include <format>
#include <stdexcept>
#include <string_view>

namespace {
void assert_stokvec(Stokvec actual,
                    Stokvec expected,
                    Numeric tol,
                    std::string_view name) {
  for (Size i = 0; i < 4; ++i) {
    if (std::abs(actual[i] - expected[i]) > tol)
      throw std::runtime_error(
          std::format("{} mismatch at {}:\nactual {} expected {}",
                      name,
                      i,
                      actual,
                      expected));
  }
}

Stokvec sum_column_weights(const sensor::Obsel& obsel, Size ifreq) try {
  Stokvec sum{0.0, 0.0, 0.0, 0.0};
  for (Size isample = 0; isample < obsel.poslos_grid().size(); ++isample) {
    sum += obsel.weight_matrix()[isample, ifreq];
  }
  return sum;
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "sum_column_weights failed for ifreq {}:\n{}", ifreq, e.what()));
}

void test_gaussian_initialization_uses_los_offset() try {
  const Stokvec peak_weight{2.0, 1.0, 0.0, 0.0};
  const sensor::GaussianAntenna ant{ZenGrid{{0.0, 1.0}}, 1.0, 2, peak_weight};

  assert_stokvec(sum(ant.data.data), 1.0, 1e-12, "gaussian peak weight");

  const Numeric scale    = std::exp(-0.5);
  const Stokvec expected = scale / (2 + 2 * scale);
  assert_stokvec(ant.data[1, 0], expected, 1e-12, "gaussian local [1, 0]");
  assert_stokvec(ant.data[1, 1], expected, 1e-12, "gaussian local [1, 1]");

  const sensor::GaussianAntenna def{ZenGrid{{0.0}}, 1.0, 1};
  assert_stokvec(def.data[0, 0],
                 {1.0, 0.0, 0.0, 0.0},
                 1e-12,
                 "default gaussian peak weight");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_gaussian_initialization_uses_los_offset failed:\n{}", e.what()));
}

Numeric gaussian_airy_expected_gain(Numeric zenith_deg,
                                    Numeric frequency,
                                    Numeric aperture_diameter) try {
  constexpr Numeric gaussian_airy_hwhm_factor =
      Constant::bessel_j_n1_k1_zero / Constant::pi;
  const Numeric wavelength = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg = Conversion::rad2deg(gaussian_airy_hwhm_factor *
                                               wavelength / aperture_diameter);
  const Numeric ratio    = zenith_deg / hwhm_deg;

  return std::exp(-std::numbers::ln2 * ratio * ratio);
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "gaussian_airy_expected_gain failed for zenith {} frequency {} aperture {}:\n{}",
      zenith_deg,
      frequency,
      aperture_diameter,
      e.what()));
}

void test_gaussian_airy_is_frequency_dependent() try {
  const Stokvec peak_weight{2.0, 0.0, 0.0, 0.0};
  const sensor::GaussianAiryAntenna ant{
      ZenGrid{{0.0, 0.2}}, 1.0, 1, peak_weight};
  const sensor::BoxChannel channel{AscendingGrid{100.0e9, 200.0e9}};

  const auto obsel = ant(channel,
                         {600e3, 10.0, 20.0},
                         {45.0, 30.0},
                         {Body::Earth::a, Body::Earth::b});

  const auto low_weight  = obsel.weight_matrix()[1, 0];
  const auto high_weight = obsel.weight_matrix()[1, 1];
  if (low_weight[0] <= high_weight[0])
    throw std::runtime_error(
        "Higher frequency Gaussian Airy sample must be narrower");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_gaussian_airy_is_frequency_dependent failed:\n{}", e.what()));
}

void test_degenerate_azimuth_ring_is_collapsed() try {
  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0, 3};
  const auto gaussian_obsel = gaussian(sensor::DiracChannel{100.0e9},
                                       {600e3, 10.0, 20.0},
                                       {45.0, 30.0},
                                       {Body::Earth::a, Body::Earth::b});

  if (gaussian_obsel.poslos_grid().size() != 4) {
    throw std::runtime_error(std::format(
        "Expected zero-zenith azimuth collapse to leave 4 LOS samples, got {}",
        gaussian_obsel.poslos_grid().size()));
  }

  assert_stokvec(
      gaussian_obsel.weight_matrix()[0, 0],
      {gaussian.data[0, 0] / (1 - 2 * gaussian.data[0, 0]), 0.0, 0.0, 0.0},
      1e-12,
      "collapsed gaussian zero-zenith ring weight");

  const sensor::GaussianAiryAntenna airy{ZenGrid{{0.0, 0.2}}, 1.0, 3};
  const auto airy_obsel        = airy(sensor::DiracChannel{100.0e9},
                                      {600e3, 10.0, 20.0},
                                      {45.0, 30.0},
                                      {Body::Earth::a, Body::Earth::b});
  const Numeric gain           = gaussian_airy_expected_gain(0.2, 100.0e9, 1.0);
  const Numeric retained_scale = 3.0 * (1.0 + gain) / (1.0 + 3.0 * gain);

  if (airy_obsel.poslos_grid().size() != 4) {
    throw std::runtime_error(std::format(
        "Expected Gaussian Airy zero-zenith azimuth collapse to leave 4 LOS samples, got {}",
        airy_obsel.poslos_grid().size()));
  }

  assert_stokvec(
      airy_obsel.weight_matrix()[0, 0],
      {retained_scale / (3.0 * (1.0 + gain)), 0.0, 0.0, 0.0},
      1e-6,
      "collapsed gaussian airy center is renormalized with retained points");
  assert_stokvec(
      airy_obsel.weight_matrix()[1, 0],
      {retained_scale * gain / (3.0 * (1.0 + gain)), 0.0, 0.0, 0.0},
      1e-6,
      "collapsed gaussian airy off-axis weight is renormalized uniformly");
  assert_stokvec(sum_column_weights(airy_obsel, 0),
                 {1.0, 0.0, 0.0, 0.0},
                 1e-6,
                 "collapsed gaussian airy column renormalizes to unity");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_degenerate_azimuth_ring_is_collapsed failed:\n{}", e.what()));
}

void test_gaussian_does_not_double_count_degenerate_azimuth() try {
  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0};
  const auto obsel = gaussian(sensor::DiracChannel{100.0e9},
                              {600e3, 10.0, 20.0},
                              {45.0, 0.0},
                              {Body::Earth::a, Body::Earth::b});

  if (obsel.poslos_grid().size() != 3) {
    throw std::runtime_error(std::format(
        "Expected degenerate azimuth collapse to leave 3 LOS samples, got {}",
        obsel.poslos_grid().size()));
  }

  const Numeric wing =
      0.5 * std::exp(-0.5) / (1 + std::exp(-0.5)) / (1 - gaussian.data[0, 0]);
  assert_stokvec(
      obsel.weight_matrix()[0, 0],
      {gaussian.data[0, 0] / (1 - gaussian.data[0, 0]), 0.0, 0.0, 0.0},
      1e-12,
      "gaussian center weight should be counted once");
  assert_stokvec(obsel.weight_matrix()[1, 0],
                 {wing, 0.0, 0.0, 0.0},
                 1e-12,
                 "gaussian wing at azimuth 0");
  assert_stokvec(obsel.weight_matrix()[2, 0],
                 {wing, 0.0, 0.0, 0.0},
                 1e-12,
                 "gaussian wing at azimuth 180");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_gaussian_does_not_double_count_degenerate_azimuth failed:\n{}",
      e.what()));
}

void test_gaussian_airy_rejects_nonpositive_frequencies() try {
  const sensor::GaussianAiryAntenna ant{ZenGrid{{0.0}}, 1.0, 1};

  bool threw = false;
  try {
    static_cast<void>(ant(sensor::DiracChannel{},
                          {600e3, 10.0, 20.0},
                          {0.0, 0.0},
                          {Body::Earth::a, Body::Earth::b}));
  } catch (const std::runtime_error&) {
    threw = true;
  }

  if (not threw)
    throw std::runtime_error(
        "Gaussian Airy antenna must reject nonpositive channel frequencies");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_gaussian_airy_rejects_nonpositive_frequencies failed:\n{}",
      e.what()));
}
}  // namespace

int main() try {
  test_gaussian_initialization_uses_los_offset();
  test_gaussian_airy_is_frequency_dependent();
  test_degenerate_azimuth_ring_is_collapsed();
  test_gaussian_does_not_double_count_degenerate_azimuth();
  test_gaussian_airy_rejects_nonpositive_frequencies();
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Test failed:\n{}", e.what());
  return EXIT_FAILURE;
}
