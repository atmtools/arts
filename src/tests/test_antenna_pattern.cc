#include <antenna_pattern.h>

#include <arts_conversions.h>

#include <cmath>
#include <stdexcept>
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

Stokvec sum_column_weights(const sensor::Obsel& obsel, Size ifreq) {
  Stokvec sum{0.0, 0.0, 0.0, 0.0};
  for (Size isample = 0; isample < obsel.poslos_grid().size(); ++isample) {
    sum += obsel.weight_matrix()[isample, ifreq];
  }
  return sum;
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

Numeric gaussian_airy_expected_gain(Numeric zenith_deg,
                                    Numeric frequency,
                                    Numeric aperture_diameter) {
  constexpr Numeric gaussian_airy_hwhm_factor =
      3.8317059702075123156 / Constant::pi;
  const Numeric wavelength = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg   = Conversion::rad2deg(
      gaussian_airy_hwhm_factor * wavelength / aperture_diameter);
  const Numeric ratio = zenith_deg / hwhm_deg;

  return std::exp(-std::log(2.0) * ratio * ratio);
}

Numeric gaussian_airy_expected_std(Numeric frequency,
                                   Numeric aperture_diameter) {
  constexpr Numeric gaussian_airy_hwhm_factor =
      3.8317059702075123156 / Constant::pi;
  const Numeric wavelength = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg   = Conversion::rad2deg(
      gaussian_airy_hwhm_factor * wavelength / aperture_diameter);

  return hwhm_deg / std::sqrt(2.0 * std::log(2.0));
}

void test_gaussian_airy_is_frequency_dependent() {
  const Stokvec peak_weight{2.0, 0.0, 0.0, 0.0};
  const sensor::GaussianAiryAntenna ant{
      ZenGrid{{0.0, 0.2}}, AziGrid{{0.0}}, 1.0, peak_weight};
  const sensor::BoxChannel channel{AscendingGrid{100.0e9, 200.0e9}};

  const auto obsel = ant(channel, {600e3, 10.0, 20.0}, {45.0, 30.0});

  const Numeric low_gain  = gaussian_airy_expected_gain(0.2, 100.0e9, 1.0);
  const Numeric high_gain = gaussian_airy_expected_gain(0.2, 200.0e9, 1.0);

    const Numeric low_norm  = 1.0 + low_gain;
    const Numeric high_norm = 1.0 + high_gain;

    assert_stokvec(obsel.weight_matrix()[0, 0],
           (0.5 / low_norm) * peak_weight,
                 1e-12,
                 "gaussian airy bore low frequency");
  assert_stokvec(obsel.weight_matrix()[0, 1],
           (0.5 / high_norm) * peak_weight,
                 1e-12,
                 "gaussian airy bore high frequency");

  assert_stokvec(obsel.weight_matrix()[1, 0],
           (0.5 * low_gain / low_norm) * peak_weight,
                 1e-12,
                 "gaussian airy off-axis low frequency");
  assert_stokvec(obsel.weight_matrix()[1, 1],
           (0.5 * high_gain / high_norm) * peak_weight,
                 1e-12,
                 "gaussian airy off-axis high frequency");

    assert_stokvec(
      sum_column_weights(obsel, 0), 0.5 * peak_weight, 1e-12, "gaussian airy normalized low frequency");
    assert_stokvec(
      sum_column_weights(obsel, 1), 0.5 * peak_weight, 1e-12, "gaussian airy normalized high frequency");

  const auto low_weight  = obsel.weight_matrix()[1, 0];
  const auto high_weight = obsel.weight_matrix()[1, 1];
  ARTS_USER_ERROR_IF(low_weight[0] <= high_weight[0],
                     "Higher frequency Gaussian Airy sample must be narrower")
}

void test_gaussian_airy_matches_frequency_specific_gaussian_pattern() {
  const ZenGrid zen_grid{{0.1, 0.2}};
  const AziGrid azi_grid{{0.0, 0.2}};
  const Stokvec peak_weight{2.0, 0.0, 0.0, 0.0};
  const sensor::GaussianAiryAntenna ant{zen_grid, azi_grid, 1.0, peak_weight};
  const sensor::BoxChannel channel{AscendingGrid{100.0e9, 200.0e9}};

  const auto obsel = ant(channel, {600e3, 10.0, 20.0}, {45.0, 30.0});

  for (Size ifreq = 0; ifreq < channel.freq_grid().size(); ++ifreq) {
    const Numeric airy_std =
        gaussian_airy_expected_std(channel.freq_grid()[ifreq], 1.0);
    const sensor::GaussianAntenna gaussian{
        zen_grid, azi_grid, airy_std, airy_std, peak_weight};

    Numeric normalization = 0.0;
    for (Size izen = 0; izen < zen_grid.size(); ++izen) {
      for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
        normalization += gaussian.data[izen, iazi][0] / peak_weight[0];
      }
    }

    Size isample = 0;
    for (Size izen = 0; izen < zen_grid.size(); ++izen) {
      for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
        assert_stokvec(obsel.weight_matrix()[isample, ifreq],
                       (channel.weights()[ifreq] / normalization) *
                           gaussian.data[izen, iazi],
                       1e-12,
                       "gaussian airy per-frequency gaussian match");
        ++isample;
      }
    }

    assert_stokvec(sum_column_weights(obsel, ifreq),
                   channel.weights()[ifreq] * peak_weight,
                   1e-12,
                   "gaussian airy per-frequency column normalization");
  }
}

          void test_degenerate_azimuth_ring_is_collapsed() {
            const sensor::GaussianAntenna gaussian{
              ZenGrid{{0.0, 1.0}}, AziGrid{{0.0, 120.0, 240.0}}, 1.0, 1.0};
            const auto gaussian_obsel =
              gaussian(sensor::DiracChannel{100.0e9}, {600e3, 10.0, 20.0}, {45.0, 30.0});

            ARTS_USER_ERROR_IF(gaussian_obsel.poslos_grid().size() != 4,
                     "Expected zero-zenith azimuth collapse to leave 4 LOS samples, got {}",
                     gaussian_obsel.poslos_grid().size())
            assert_stokvec(gaussian_obsel.weight_matrix()[0, 0],
                   {3.0, 0.0, 0.0, 0.0},
                   1e-12,
                   "collapsed gaussian zero-zenith ring weight");

            const sensor::GaussianAiryAntenna airy{
              ZenGrid{{0.0, 0.2}}, AziGrid{{0.0, 120.0, 240.0}}, 1.0};
            const auto airy_obsel =
              airy(sensor::DiracChannel{100.0e9}, {600e3, 10.0, 20.0}, {45.0, 30.0});
            const Numeric gain = gaussian_airy_expected_gain(0.2, 100.0e9, 1.0);
            const Numeric retained_scale = 3.0 * (1.0 + gain) / (1.0 + 3.0 * gain);

            ARTS_USER_ERROR_IF(airy_obsel.poslos_grid().size() != 4,
                     "Expected Gaussian Airy zero-zenith azimuth collapse to leave 4 LOS samples, got {}",
                     airy_obsel.poslos_grid().size())
            assert_stokvec(airy_obsel.weight_matrix()[0, 0],
                           {retained_scale / (3.0 * (1.0 + gain)), 0.0, 0.0, 0.0},
                           1e-6,
                           "collapsed gaussian airy center is renormalized with retained points");
            assert_stokvec(airy_obsel.weight_matrix()[1, 0],
                           {retained_scale * gain / (3.0 * (1.0 + gain)), 0.0, 0.0, 0.0},
                           1e-6,
                           "collapsed gaussian airy off-axis weight is renormalized uniformly");
            assert_stokvec(sum_column_weights(airy_obsel, 0),
                           {1.0, 0.0, 0.0, 0.0},
                           1e-6,
                           "collapsed gaussian airy column renormalizes to unity");
          }

void test_gaussian_airy_rejects_nonpositive_frequencies() {
  const sensor::GaussianAiryAntenna ant{ZenGrid{{0.0}}, AziGrid{{0.0}}, 1.0};

  bool threw = false;
  try {
    static_cast<void>(
        ant(sensor::DiracChannel{}, {600e3, 10.0, 20.0}, {0.0, 0.0}));
  } catch (const std::runtime_error&) {
    threw = true;
  }

  ARTS_USER_ERROR_IF(not threw,
                     "Gaussian Airy antenna must reject nonpositive channel frequencies")
}
}  // namespace

int main() {
  test_gaussian_initialization_uses_antenna_frame_offsets();
  test_gaussian_airy_is_frequency_dependent();
  test_gaussian_airy_matches_frequency_specific_gaussian_pattern();
  test_degenerate_azimuth_ring_is_collapsed();
  test_gaussian_airy_rejects_nonpositive_frequencies();
  return 0;
}