#include <antenna_pattern.h>
#include <arts_constants.h>
#include <arts_conversions.h>
#include <planet_data.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <stdexcept>
#include <string_view>

namespace {
const Vector2 ell{Body::Earth::a, Body::Earth::b};
const Vector3 pos{600e3, 10.0, 20.0};
const Stokvec unpolarized{1.0, 0.0, 0.0, 0.0};

void assert_close(Numeric actual, Numeric expected, Numeric tol, std::string_view name) {
  if (std::abs(actual - expected) > tol)
    throw std::runtime_error(std::format("{} mismatch: actual {} expected {}", name, actual, expected));
}

void assert_stokvec(Stokvec actual, Stokvec expected, Numeric tol, std::string_view name) {
  for (Size i = 0; i < 4; ++i) {
    if (std::abs(actual[i] - expected[i]) > tol)
      throw std::runtime_error(std::format("{} mismatch at {}:\nactual {} expected {}", name, i, actual, expected));
  }
}

//! Solid angle [sr] of the spherical cap within a half angle [deg]
Numeric cap(Numeric half_angle) { return Constant::two_pi * (1.0 - Conversion::cosd(half_angle)); }

//! Unit vector of a line of sight in the local east-north-up frame
Vector3 enu(const Vector2& los) {
  using Conversion::cosd, Conversion::sind;
  return {sind(los[0]) * sind(los[1]), sind(los[0]) * cosd(los[1]), cosd(los[0])};
}

//! Angle [deg] between two lines of sight at the same position
Numeric offset(const Vector2& a, const Vector2& b) {
  return Conversion::acosd(std::clamp(dot(enu(a), enu(b)), -1.0, 1.0));
}

Numeric gaussian_airy_expected_gain(Numeric zenith_deg, Numeric frequency, Numeric aperture_diameter) {
  constexpr Numeric gaussian_airy_hwhm_factor = Constant::bessel_j_n1_k1_zero / Constant::pi;
  const Numeric     wavelength                = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg = Conversion::rad2deg(gaussian_airy_hwhm_factor * wavelength / aperture_diameter);
  const Numeric ratio    = zenith_deg / hwhm_deg;
  return std::exp(-std::numbers::ln2 * ratio * ratio);
}

void test_gaussian_rejects_nonpositive_std() try {
  const sensor::GaussianAntenna ant{ZenGrid{{0.0, 1.0}}, 0.0, 2};

  bool threw = false;
  try {
    static_cast<void>(ant(sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 30.0}, ell));
  } catch (const std::runtime_error&) { threw = true; }

  if (not threw) throw std::runtime_error("Gaussian antenna must reject a nonpositive std");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_rejects_nonpositive_std failed:\n{}", e.what()));
}

void test_gaussian_airy_is_frequency_dependent() try {
  const sensor::GaussianAiryAntenna ant{ZenGrid{{0.0, 0.2}}, 1.0, 1};
  const sensor::BoxChannel          channel{AscendingGrid{100.0e9, 200.0e9}};

  const auto obsel = ant(channel, unpolarized, pos, {45.0, 30.0}, ell);

  const auto low_weight  = obsel.weight_matrix()[1, 0];
  const auto high_weight = obsel.weight_matrix()[1, 1];

  if (low_weight[0] <= high_weight[0])
    throw std::runtime_error("Higher frequency Gaussian Airy sample must be narrower");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_airy_is_frequency_dependent failed:\n{}", e.what()));
}

//! The beam weight of each sample is its gain times the solid angle of its cell: the cap inside the first
//! midpoint for the boresight, and an equal share of the band between the midpoints for each ring sample
void test_cells_are_spherical_caps_and_bands() try {
  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0, 3};
  const auto gaussian_obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 30.0}, ell);

  if (gaussian_obsel.poslos_grid().size() != 4) {
    throw std::runtime_error(std::format("Expected zero-zenith azimuth collapse to leave 4 LOS samples, got {}",
                                         gaussian_obsel.poslos_grid().size()));
  }

  const Numeric center = cap(0.5);
  const Numeric band   = (cap(1.5) - cap(0.5)) / 3.0 * std::exp(-0.5);
  const Numeric total  = center + 3.0 * band;
  assert_close(gaussian_obsel.weight_matrix()[0, 0][0], center / total, 1e-12, "gaussian boresight cap");
  for (Size i = 1; i < 4; i++) {
    assert_close(gaussian_obsel.weight_matrix()[i, 0][0], band / total, 1e-12, "gaussian band sample");
  }

  constexpr Numeric                 f = 100.0e9;
  const sensor::GaussianAiryAntenna airy{ZenGrid{{0.0, 0.2}}, 1.0, 3};
  const auto                        airy_obsel = airy(sensor::DiracChannel{f}, unpolarized, pos, {45.0, 30.0}, ell);

  const Numeric airy_center = cap(0.1);
  const Numeric airy_band   = (cap(0.3) - cap(0.1)) / 3.0 * gaussian_airy_expected_gain(0.2, f, 1.0);
  const Numeric airy_total  = airy_center + 3.0 * airy_band;
  assert_close(airy_obsel.weight_matrix()[0, 0][0], airy_center / airy_total, 1e-6, "airy boresight cap");
  assert_close(airy_obsel.weight_matrix()[1, 0][0], airy_band / airy_total, 1e-6, "airy band sample");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_cells_are_spherical_caps_and_bands failed:\n{}", e.what()));
}

void test_two_azimuths_split_the_band() try {
  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0};
  const auto                    obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 0.0}, ell);

  if (obsel.poslos_grid().size() != 3) {
    throw std::runtime_error(
        std::format("Expected degenerate azimuth collapse to leave 3 LOS samples, got {}", obsel.poslos_grid().size()));
  }

  const Numeric center = cap(0.5);
  const Numeric wing   = (cap(1.5) - cap(0.5)) / 2.0 * std::exp(-0.5);
  const Numeric total  = center + 2.0 * wing;
  assert_close(obsel.weight_matrix()[0, 0][0], center / total, 1e-12, "gaussian center weight counted once");
  assert_close(obsel.weight_matrix()[1, 0][0], wing / total, 1e-12, "gaussian wing at azimuth 0");
  assert_close(obsel.weight_matrix()[2, 0][0], wing / total, 1e-12, "gaussian wing at azimuth 180");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_two_azimuths_split_the_band failed:\n{}", e.what()));
}

/** The beam reproduces the second moment of a 2D Gaussian
 *
 *  For the gain exp(-r^2 / 2 s^2) on the plane, truncated at r = R, the mean of r^2 is
 *  2 s^2 (1 - x exp(-x) / (1 - exp(-x))) with x = R^2 / 2 s^2.  The sphere differs from the plane by
 *  sin(r) / r - 1, about 1e-4 within the 2 degrees used here.
 */
void test_gaussian_second_moment() try {
  constexpr Numeric std  = 0.5;
  constexpr Numeric step = 0.025;
  const Vector2     bore{135.0, 20.0};

  const sensor::GaussianAntenna gaussian{ZenGrid{nlinspace(0.0, 2.0, 81)}, std, 72};
  const auto                    obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, bore, ell);

  Numeric m2 = 0.0, total = 0.0;
  for (Size i = 0; i < obsel.poslos_grid().size(); i++) {
    const Numeric w  = obsel.weight_matrix()[i, 0][0];
    const Numeric r  = offset(obsel.poslos_grid()[i].los, bore);
    m2              += w * r * r;
    total           += w;
  }

  const Numeric R        = 2.0 + 0.5 * step;
  const Numeric x        = R * R / (2.0 * std * std);
  const Numeric expected = 2.0 * std * std * (1.0 - x * std::exp(-x) / (1.0 - std::exp(-x)));

  assert_close(total, 1.0, 1e-12, "gaussian beam weights sum");
  assert_close(m2 / expected, 1.0, 2e-3, "gaussian beam second moment relative to the analytic");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_second_moment failed:\n{}", e.what()));
}

//! Each frequency of the channel keeps its channel weight; the beam is normalized per frequency.  The channel
//! weights of a box on 10 equidistant frequencies are the trapezoid rule's, 1/18 at the ends and 1/9 inside.
void test_gaussian_airy_keeps_the_channel_weights() try {
  const sensor::GaussianAiryAntenna ant{ZenGrid{nlinspace(0.0, 3.0, 13)}, 0.3, 24};
  const sensor::BoxChannel          channel{50e9, 500e9, 10};

  const auto obsel = ant(channel, unpolarized, pos, {135.0, 0.0}, ell);

  for (Size j = 0; j < obsel.f_grid().size(); j++) {
    Numeric marginal = 0.0;
    for (Size i = 0; i < obsel.poslos_grid().size(); i++) marginal += obsel.weight_matrix()[i, j][0];
    const Numeric expected = (j == 0 or j + 1 == obsel.f_grid().size()) ? 1.0 / 18.0 : 1.0 / 9.0;
    assert_close(marginal, expected, 1e-12, std::format("airy frequency marginal at {} Hz", obsel.f_grid()[j]));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_airy_keeps_the_channel_weights failed:\n{}", e.what()));
}

//! A cross-track scanner turns its linear polarization by scan_rotation times the signed scan angle
void test_scan_rotation() try {
  using Conversion::cosd, Conversion::sind;

  const sensor::Polarization vertical{.state = {1.0, 1.0, 0.0, 0.0}, .scan_rotation = 1.0};
  assert_stokvec(vertical.at(30.0), {1.0, cosd(60.0), sind(60.0), 0.0}, 1e-14, "vertical at nadir, +30 degrees");
  assert_stokvec(vertical.at(-30.0), {1.0, cosd(60.0), -sind(60.0), 0.0}, 1e-14, "vertical at nadir, -30 degrees");

  const sensor::Polarization opposite{.state = {1.0, 1.0, 0.0, 0.0}, .scan_rotation = -1.0};
  assert_stokvec(opposite.at(30.0), {1.0, cosd(60.0), -sind(60.0), 0.0}, 1e-14, "opposite optics, +30 degrees");

  const sensor::Polarization horizontal{.state = {1.0, -1.0, 0.0, 0.0}, .scan_rotation = 1.0};
  assert_stokvec(horizontal.at(30.0), {1.0, -cosd(60.0), -sind(60.0), 0.0}, 1e-14, "horizontal at nadir, +30 degrees");

  const sensor::Polarization fixed{.state = {1.0, 1.0, 0.0, 0.0}};
  assert_stokvec(fixed.at(30.0), {1.0, 1.0, 0.0, 0.0}, 0.0, "no scan rotation");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_scan_rotation failed:\n{}", e.what()));
}

/** A nadir-looking antenna receives one field direction in all its lines of sight
 *
 *  The boresight [180, 0] receives radiation going up, with the basis h = (-1, 0, 0) of its azimuth and
 *  v = h x k = (0, 1, 0), north.  A line of sight at a small angle r from nadir towards azimuth a
 *  receives along k = (-sin r sin a, -sin r cos a, cos r), with h = (-cos a, sin a, 0) and v = h x k,
 *  which is (sin a, cos a, 0) to first order in r.  The northward field has E_v = cos a and E_h = sin a
 *  there, so its Stokes weights are [1, cos 2a, sin 2a, 0], to second order in r.
 */
void test_polarization_is_carried_around_nadir() try {
  using Conversion::cosd, Conversion::sind;

  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0, 8};
  sensor::Channel               channel = sensor::DiracChannel{100.0e9};
  channel.polarization.state            = {1.0, 1.0, 0.0, 0.0};

  const auto obsel = gaussian(channel, channel.polarization.state, pos, {180.0, 0.0}, ell);

  for (Size i = 0; i < obsel.poslos_grid().size(); i++) {
    const Stokvec w = obsel.weight_matrix()[i, 0];
    const Numeric a = obsel.poslos_grid()[i].los[1];
    assert_stokvec(w / w.I(),
                   {1.0, cosd(2.0 * a), sind(2.0 * a), 0.0},
                   1e-3,
                   std::format("northward field in the line of sight towards azimuth {}", a));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_polarization_is_carried_around_nadir failed:\n{}", e.what()));
}

void test_normalize_keeps_the_polarization() try {
  sensor::SparseStokvecMatrix weights(1, 2);
  weights[0, 0] = Stokvec{1.0, 0.6, 0.8, 0.0};
  weights[0, 1] = Stokvec{1.0, 0.6, 0.8, 0.0};
  sensor::Obsel obsel{AscendingGrid{100e9, 101e9}, SensorPosLosVector(Size{1}, sensor::PosLos{.pos = pos, .los = Vector2{150.0, 0.0}}), weights};

  obsel.normalize();
  assert_stokvec(obsel.weight_matrix()[0, 0], {0.5, 0.3, 0.4, 0.0}, 1e-15, "normalized rotated linear weights");

  sensor::SparseStokvecMatrix unpolarized(1, 1);
  unpolarized[0, 0] = Stokvec{1.0, 0.0, 0.0, 0.0};
  sensor::Obsel plain{AscendingGrid{100e9}, SensorPosLosVector(Size{1}, sensor::PosLos{.pos = pos, .los = Vector2{150.0, 0.0}}), unpolarized};

  bool threw = false;
  try {
    plain.normalize(Stokvec{1.0, 1.0, 0.0, 0.0});
  } catch (const std::runtime_error&) { threw = true; }
  if (not threw) throw std::runtime_error("Normalizing unpolarized weights to a polarization must throw");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_normalize_keeps_the_polarization failed:\n{}", e.what()));
}

void test_gaussian_airy_rejects_nonpositive_frequencies() try {
  const sensor::GaussianAiryAntenna ant{ZenGrid{{0.0}}, 1.0, 1};

  bool threw = false;
  try {
    static_cast<void>(ant(sensor::DiracChannel{}, unpolarized, pos, {0.0, 0.0}, ell));
  } catch (const std::runtime_error&) { threw = true; }

  if (not threw) throw std::runtime_error("Gaussian Airy antenna must reject nonpositive channel frequencies");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_airy_rejects_nonpositive_frequencies failed:\n{}", e.what()));
}
}  // namespace

int main() try {
  test_gaussian_rejects_nonpositive_std();
  test_gaussian_airy_is_frequency_dependent();
  test_cells_are_spherical_caps_and_bands();
  test_two_azimuths_split_the_band();
  test_gaussian_second_moment();
  test_gaussian_airy_keeps_the_channel_weights();
  test_scan_rotation();
  test_polarization_is_carried_around_nadir();
  test_normalize_keeps_the_polarization();
  test_gaussian_airy_rejects_nonpositive_frequencies();
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Test failed:\n{}", e.what());
  return EXIT_FAILURE;
}
