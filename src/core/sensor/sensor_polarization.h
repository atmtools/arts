#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml_io_stream.h>

namespace sensor {
/** The polarization response of a receiver channel
 *
 *  The response is a Stokes vector of weights in the polarization basis of
 *  the boresight ray (see the polarization geometry concept): the basis of the
 *  received radiation's propagation direction, which is opposite to the line
 *  of sight.  The measurement is the dot product of the weights and the
 *  Stokes vector, so [1, 1, 0, 0] receives vertical and [1, -1, 0, 0]
 *  horizontal polarization.
 *
 *  Cross-track scanners turn the linear part of the response about the ray as
 *  they scan.  How much is a property of the instrument's optics, so it is
 *  given, not computed: the linear polarization turns by scan_rotation times
 *  the signed scan angle.  The scan angle is the nadir angle of the boresight,
 *  positive to the right of the flight direction and negative to the left
 *  (see Builder::cross_track), so scan_rotation is 0 for a response that does
 *  not depend on the scan and +1 or -1 for a cross-track scanner with a
 *  rotating reflector, the sign being that of its optics.  The turn is from v
 *  towards h in the basis of each side of the scan.
 *
 *  An antenna pattern carries the response from the boresight to its other
 *  lines of sight along the great circles to them.
 */
struct Polarization {
  Stokvec state{1.0, 0.0, 0.0, 0.0};
  Numeric scan_rotation{0.0};

  /** The response at a scan angle
   *
   *  @param scan_angle The signed scan angle [deg]
   *  @return The state with its linear part turned by scan_rotation times scan_angle
   */
  [[nodiscard]] Stokvec at(Numeric scan_angle) const;
};
}  // namespace sensor

template <> struct format_tag_aggregate<sensor::Polarization> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<sensor::Polarization> {
  static constexpr std::string_view name = "SensorPolarization";
};

template <> struct xml_io_stream_aggregate<sensor::Polarization> {
  static constexpr bool value = true;
};
