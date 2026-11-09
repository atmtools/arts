#pragma once

#include <span>
#include <vector>

#include "antenna_pattern.h"
#include "frequency_channel_selection.h"
#include "sensor.h"

namespace sensor {
struct FrequencyRange;

/** Combines channels with an antenna pattern and bore geometries into a sensor
 *
 *  The sensor simulates the frequencies of all channels.  Each placement of the
 *  instrument, a position with a bore line of sight, adds the lines of sight
 *  its antenna samples.  The elements are ordered by placement, then by pixel
 *  of the antenna, then by channel.  The meta has the same axes: placement by
 *  channel for an antenna with one pixel, and placement by rows by columns of
 *  pixels by channel otherwise.
 */
struct Builder {
  std::vector<Channel> channels;
  AntennaPattern       antenna{};

  //! The channels of a spectrometer
  [[nodiscard]] static Builder from(const Spectrometer& spectrometer, AntennaPattern antenna = {});

  //! The channels of a spectrometer behind a backend
  [[nodiscard]] static Builder from(const Spectrometer&   spectrometer,
                                    const FrequencyRange& backend,
                                    AntennaPattern        antenna = {});

  /** The sensor of paired positions and bore lines of sight, one placement each
   *
   *  The scan side of a line of sight is not known here, so every channel must
   *  have a polarization that does not depend on the scan (scan_rotation ==
   *  0); use cross_track() otherwise.  The sensor is named name, or after its
   *  number of channels if name is empty.
   */
  [[nodiscard]] Sensor operator()(std::span<const Vector3> pos,
                                  std::span<const Vector2> los,
                                  const Vector2&           ell,
                                  const String&            name = "") const;

  /** The sensor of cross-track scans
   *
   *  Every position scans the same signed scan angles [deg] in the plane
   *  across its flight direction, heading [deg] clockwise from north.  A scan
   *  angle is the nadir angle of the boresight, positive to the right of the
   *  flight direction: the bore line of sight is [180 - |s|, heading + 90] for
   *  s >= 0 and [180 - |s|, heading - 90] for s < 0.  Each channel's
   *  polarization is taken at s, see Polarization::at.  The placements are
   *  ordered by position and then scan angle.  The antenna must have one
   *  pixel.
   *
   *  @param pos         Platform positions
   *  @param heading     Flight direction at each position
   *  @param scan_angles Signed scan angles, the same for all positions
   *  @param ell         Reference ellipsoid
   *  @param name        The name of the sensor, or empty to name it after its number of channels
   */
  [[nodiscard]] Sensor cross_track(std::span<const Vector3> pos,
                                   std::span<const Numeric> heading,
                                   const Vector&            scan_angles,
                                   const Vector2&           ell,
                                   const String&            name = "") const;
};
}  // namespace sensor

using SensorBuilder = sensor::Builder;

template <> struct format_tag_aggregate<sensor::Builder> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<sensor::Builder> {
  static constexpr std::string_view name = "SensorBuilder";
};

template <> struct xml_io_stream_aggregate<sensor::Builder> {
  static constexpr bool value = true;
};
