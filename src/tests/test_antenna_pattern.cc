#include <antenna_pattern.h>
#include <arts_constants.h>
#include <arts_conversions.h>
#include <planet_data.h>

#include <algorithm>
#include <cmath>
#include <format>
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

//! The Airy pattern (2 J1(x) / x)^2 of a uniformly illuminated circular aperture
Numeric airy_pattern(Numeric x) {
  const Numeric a = 2.0 * std::cyl_bessel_j(1.0, x) / x;
  return a * a;
}

//! Each current antenna kind is a single beam: one pixel at no offset, and one obsel per channel
void test_single_beams_have_one_pixel() try {
  sensor::GriddedAntennaPattern gridded;
  gridded.data = sensor::AntennaPatternField{
      .data_name  = std::string{"gain"},
      .data       = Matrix(2, 2, 1.0),
      .grid_names = {std::string{"zenith"}, std::string{"azimuth"}},
      .grids      = {ZenGrid{{0.0, 1.0}}, AziGrid{{0.0, 180.0}}},
  };

  const std::vector<sensor::AntennaPattern> antennas{sensor::PencilBeamAntenna{},
                                                     sensor::GaussianAntenna{ZenGrid{{0.0, 1.0}}, 1.0, 2},
                                                     sensor::GaussianAiryAntenna{ZenGrid{{0.0, 0.2}}, 1.0, 2},
                                                     gridded};

  for (const auto& antenna : antennas) {
    const auto pixels = sensor::antenna_pixels(antenna);
    if (pixels.size() != 1 or pixels.rows.size() != 1 or pixels.cols.size() != 1)
      throw std::runtime_error(std::format("Antenna {} must have one pixel", antenna.index()));
    assert_close(pixels.rows[0], 0.0, 0.0, "pixel row offset");
    assert_close(pixels.cols[0], 0.0, 0.0, "pixel column offset");

    const auto obsels =
        sensor::antenna_obsels(antenna, sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 30.0}, ell);
    if (obsels.size() != pixels.size())
      throw std::runtime_error(std::format(
          "Antenna {} must give one obsel per pixel, got {} for {}", antenna.index(), obsels.size(), pixels.size()));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_single_beams_have_one_pixel failed:\n{}", e.what()));
}

/** The pixels of a thin-lens camera look through the centre of the lens, inverted, in CCD coordinates
 *
 *  For a camera looking north at the horizon the local east-north-up frame has the optical axis (0, 1, 0), and a
 *  pixel centred at (dx, dy) on a CCD at the image distance d_i behind the lens looks along (-dx, d_i, -dy): the
 *  ray from the pixel through the centre of the lens.  Its zenith angle is acos(-dy / |v|) and its azimuth
 *  atan2(-dx, d_i).  The pixels are ordered by row and then by column.
 */
void test_camera_pixel_directions() try {
  constexpr Size    n_rows = 3, n_cols = 4;
  constexpr Numeric ccd_h = 0.03, ccd_w = 0.02, f = 0.05, d_o = 100.0;
  const Numeric     d_i = f * d_o / (d_o - f);

  const auto camera = sensor::CameraAntenna::thin_lens(n_rows, n_cols, ccd_h, ccd_w, f, d_o);
  const auto obsels = camera(sensor::DiracChannel{100.0e9}, unpolarized, pos, {90.0, 0.0}, ell);

  if (camera.pixels().size() != n_rows * n_cols or obsels.size() != n_rows * n_cols)
    throw std::runtime_error(std::format("Expected {} pixels and obsels", n_rows * n_cols));

  for (Size i = 0; i < n_rows; i++) {
    for (Size j = 0; j < n_cols; j++) {
      const Numeric dy = (static_cast<Numeric>(i) - 0.5 * (n_rows - 1)) * ccd_h / n_rows;
      const Numeric dx = (static_cast<Numeric>(j) - 0.5 * (n_cols - 1)) * ccd_w / n_cols;
      const Numeric v  = std::hypot(dx, d_i, dy);

      const auto& obsel = obsels[i * n_cols + j];
      if (obsel.poslos_grid().size() != 1) throw std::runtime_error("A pencil-beam pixel samples one line of sight");
      const Vector2 los = obsel.poslos_grid()[0].los;
      assert_close(los[0], Conversion::acosd(-dy / v), 1e-12, std::format("zenith of pixel ({}, {})", i, j));
      assert_close(los[1], Conversion::atan2d(-dx, d_i), 1e-12, std::format("azimuth of pixel ({}, {})", i, j));
      assert_close(camera.grid.rows[i], Conversion::atan2d(dy, d_i), 1e-12, "row offset");
      assert_close(camera.grid.cols[j], Conversion::atan2d(dx, d_i), 1e-12, "column offset");
    }
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_camera_pixel_directions failed:\n{}", e.what()));
}

//! A nadir camera receives one field direction in all its pixels, as in test_polarization_is_carried_around_nadir
void test_camera_carries_polarization_at_nadir() try {
  using Conversion::cosd, Conversion::sind;

  const sensor::CameraAntenna camera{.grid = {.rows = AscendingGrid{-0.5, 0.0, 0.5}, .cols = AscendingGrid{-0.5, 0.0, 0.5}}};
  const auto obsels = camera(sensor::DiracChannel{100.0e9}, Stokvec{1.0, 1.0, 0.0, 0.0}, pos, {180.0, 0.0}, ell);

  for (const auto& obsel : obsels) {
    const Stokvec w = obsel.weight_matrix()[0, 0];
    const Numeric a = obsel.poslos_grid()[0].los[1];
    assert_stokvec(w / w.I(),
                   {1.0, cosd(2.0 * a), sind(2.0 * a), 0.0},
                   1e-3,
                   std::format("northward field in the pixel towards azimuth {}", a));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_camera_carries_polarization_at_nadir failed:\n{}", e.what()));
}

//! Each pixel of a camera has the camera's beam, centred on the pixel
void test_camera_pixels_have_the_beam() try {
  const sensor::CameraAntenna pencil{.grid = {.rows = AscendingGrid{Vector{0.0}}, .cols = AscendingGrid{-1.0, 1.0}}};
  sensor::CameraAntenna       gaussian = pencil;
  gaussian.beam                        = sensor::GaussianAntenna{ZenGrid{{0.0, 0.5}}, 0.5, 4};

  const auto centres = pencil(sensor::DiracChannel{100.0e9}, unpolarized, pos, {120.0, 30.0}, ell);
  const auto beams   = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, {120.0, 30.0}, ell);

  for (Size p = 0; p < 2; p++) {
    if (beams[p].poslos_grid().size() != 5)
      throw std::runtime_error(std::format("Pixel {} must sample its beam's 5 lines of sight", p));

    // The zero zenith offset of the beam is the pixel's own line of sight
    assert_close(beams[p].poslos_grid()[0].los[0], centres[p].poslos_grid()[0].los[0], 1e-12, "beam centre zenith");
    assert_close(beams[p].poslos_grid()[0].los[1], centres[p].poslos_grid()[0].los[1], 1e-12, "beam centre azimuth");

    Numeric total = 0.0;
    for (Size i = 0; i < 5; i++) total += beams[p].weight_matrix()[i, 0][0];
    assert_close(total, 1.0, 1e-12, "beam weights of a pixel");
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_camera_pixels_have_the_beam failed:\n{}", e.what()));
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

  const auto obsel = ant(channel, unpolarized, pos, {45.0, 30.0}, ell).front();

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
  const auto gaussian_obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 30.0}, ell).front();

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
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_cells_are_spherical_caps_and_bands failed:\n{}", e.what()));
}

/** The Gaussian Airy beam has the half-power width of the Airy pattern of its aperture
 *
 *  For small angles the Airy pattern of an aperture of diameter D is (2 J1(x) / x)^2 with x = pi D offset /
 *  wavelength.  The standard library's Bessel function confirms the argument where it is one half, and a sample
 *  at that offset must have half the gain of the boresight: its weight is that gain times the solid angle of its
 *  cell, as for the Gaussian above.
 */
void test_gaussian_airy_half_power_width() try {
  constexpr Numeric x = Constant::airy_half_power_argument;
  assert_close(airy_pattern(x), 0.5, 1e-15, "Airy pattern at the half-power argument");

  constexpr Numeric f          = 100.0e9;
  constexpr Numeric D          = 0.1;
  const Numeric     wavelength = Constant::speed_of_light / f;
  const Numeric     half_power = Conversion::rad2deg(x * wavelength / (Constant::pi * D));

  const sensor::GaussianAiryAntenna airy{ZenGrid{{0.0, half_power}}, D, 3};
  const auto                        obsel = airy(sensor::DiracChannel{f}, unpolarized, pos, {45.0, 30.0}, ell).front();

  const Numeric center = cap(0.5 * half_power);
  const Numeric band   = (cap(1.5 * half_power) - cap(0.5 * half_power)) / 3.0 * 0.5;
  const Numeric total  = center + 3.0 * band;
  assert_close(obsel.weight_matrix()[0, 0][0], center / total, 1e-12, "airy boresight cap");
  for (Size i = 1; i < 4; i++) {
    assert_close(obsel.weight_matrix()[i, 0][0], band / total, 1e-12, "airy sample at the half-power offset");
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_gaussian_airy_half_power_width failed:\n{}", e.what()));
}

void test_two_azimuths_split_the_band() try {
  const sensor::GaussianAntenna gaussian{ZenGrid{{0.0, 1.0}}, 1.0};
  const auto                    obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, {45.0, 0.0}, ell).front();

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
  const auto                    obsel = gaussian(sensor::DiracChannel{100.0e9}, unpolarized, pos, bore, ell).front();

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

  const auto obsel = ant(channel, unpolarized, pos, {135.0, 0.0}, ell).front();

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

  const auto obsel = gaussian(channel, channel.polarization.state, pos, {180.0, 0.0}, ell).front();

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
  test_single_beams_have_one_pixel();
  test_camera_pixel_directions();
  test_camera_carries_polarization_at_nadir();
  test_camera_pixels_have_the_beam();
  test_gaussian_rejects_nonpositive_std();
  test_gaussian_airy_is_frequency_dependent();
  test_cells_are_spherical_caps_and_bands();
  test_gaussian_airy_half_power_width();
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
