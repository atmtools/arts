#include "sensor_polarization.h"

#include <arts_conversions.h>

namespace sensor {
Stokvec Polarization::at(Numeric scan_angle) const {
  if (scan_rotation == 0.0) return state;

  // Turning the receiver by chi from v towards h is the change of basis by -chi
  const Numeric chi = scan_rotation * scan_angle;
  return rtepack::stokes_rotation(Conversion::cosd(2.0 * chi), -Conversion::sind(2.0 * chi)) * state;
}
}  // namespace sensor
