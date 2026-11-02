#pragma once

#include <concepts>
#include <span>
#include <utility>
#include <vector>

#include "antenna_pattern.h"
#include "frequency_channel_selection.h"
#include "obsel.h"
#include "sensor_meta_info.h"

namespace sensor {
//! Combines channels with an antenna pattern and bore geometries into obsels.
struct Builder;
struct FrequencyRange;

//! Concept for selecting sensor builders.
template <typename T>
concept BuilderSelection = std::derived_from<T, Builder>;

struct Builder {
  std::vector<Channel> channels;
  AntennaPattern       antenna{};
  bool                 preserve_common_frequency_grid{false};

  //! The channels of a spectrometer, kept on one common frequency grid
  [[nodiscard]] static Builder from(const Spectrometer& spectrometer, AntennaPattern antenna = {});

  //! The channels of a spectrometer behind a backend, kept on one common frequency grid
  [[nodiscard]] static Builder from(const Spectrometer&   spectrometer,
                                    const FrequencyRange& backend,
                                    AntennaPattern        antenna = {});

  /** Obsels for paired positions and bore lines of sight
   *
   *  The scan side of a line of sight is not known here, so every channel
   *  must have a polarization that does not depend on the scan
   *  (scan_rotation == 0); use cross_track() otherwise.
   */
  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> operator()(std::span<const Vector3> pos,
                                                                                std::span<const Vector2> los,
                                                                                const Vector2&           ell) const;
  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> operator()(std::span<const Vector3> pos,
                                                                                std::span<const Vector2> los,
                                                                                std::span<const Time>    time,
                                                                                const Vector2&           ell) const;

  /** Obsels for cross-track scans
   *
   *  Every position scans the same signed scan angles [deg] in the plane
   *  across its flight direction, heading [deg] clockwise from north.  A scan
   *  angle is the nadir angle of the boresight, positive to the right of the
   *  flight direction: the bore line of sight is [180 - |s|, heading + 90] for
   *  s >= 0 and [180 - |s|, heading - 90] for s < 0.  Each channel's
   *  polarization is taken at s, see Polarization::at.  The obsels are ordered
   *  by position, then scan angle, then channel, with one metadata block per
   *  position and scan angle.
   *
   *  @param pos         Platform positions
   *  @param heading     Flight direction at each position
   *  @param scan_angles Signed scan angles, the same for all positions
   *  @param time        Empty or one time per position
   *  @param ell         Reference ellipsoid
   */
  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> cross_track(std::span<const Vector3> pos,
                                                                                 std::span<const Numeric> heading,
                                                                                 const Vector&            scan_angles,
                                                                                 std::span<const Time>    time,
                                                                                 const Vector2&           ell) const;
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
