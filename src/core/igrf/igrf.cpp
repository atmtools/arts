#include "igrf.h"

#include <arts_conversions.h>
#include <debug.h>
#include <geodetic.h>
#include <legendre.h>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace IGRF {
Numeric decimal_year(const Time& time) {
  using namespace std::chrono;

  Time utc       = time;
  utc.utc_offset = minutes{0};

  const year_month_day date  = utc.date();
  const Numeric        start = static_cast<Numeric>((sys_days{date} - sys_days{date.year() / January / 1}).count());
  const Numeric        day   = start + duration<Numeric, days::period>{utc.time_of_day().to_duration()}.count();
  return static_cast<int>(date.year()) + day / (date.year().is_leap() ? 366.0 : 365.0);
}

std::pair<Coefficients, Coefficients> coefficients(const Time& time) {
  using Data::epochs, Data::nepoch, Data::valid_range;

  const Numeric year = decimal_year(time);

  ARTS_USER_ERROR_IF(
      year < valid_range.first or year > valid_range.second,
      "The IGRF coefficients are valid from the decimal year {} to {}, but the time {} is the decimal year {}",
      valid_range.first,
      valid_range.second,
      time,
      year)

  std::pair<Coefficients, Coefficients> out;
  auto& [g, h] = out;

  if (year >= epochs.back()) {
    const Numeric dt     = year - epochs.back();
    const auto& [g0, h0] = Data::coefficients(nepoch - 1);
    const auto& [dg, dh] = Data::secular_variation();

    for (Size k = 0; k < Coefficients::ndata; k++) {
      g.data[k] = g0.data[k] + dt * dg.data[k];
      h.data[k] = h0.data[k] + dt * dh.data[k];
    }
  } else {
    const Size    e      = static_cast<Size>(std::ranges::upper_bound(epochs, year) - epochs.begin()) - 1;
    const Numeric w      = (year - epochs[e]) / (epochs[e + 1] - epochs[e]);
    const auto& [g0, h0] = Data::coefficients(e);
    const auto& [g1, h1] = Data::coefficients(e + 1);

    for (Size k = 0; k < Coefficients::ndata; k++) {
      g.data[k] = (1.0 - w) * g0.data[k] + w * g1.data[k];
      h.data[k] = (1.0 - w) * h0.data[k] + w * h1.data[k];
    }
  }

  return out;
}

Vector3 igrf(const Vector3 pos, const Vector2 ell, const Time& time) {
  using Conversion::deg2rad;

  const auto [g, h] = coefficients(time);

  const Vector3 geoc = geodetic2geocentric(pos, ell);
  const Vector3 mag  = Legendre::schmidt_fieldcalc(g, h, r0, geoc);

  // Rotate from the geocentric to the geodetic vertical
  const Numeric ang = deg2rad(pos[1] - geoc[1]);
  const Numeric ca  = std::cos(ang);
  const Numeric sa  = std::sin(ang);

  return {1e-9 * mag[2], 1e-9 * (-ca * mag[1] - sa * mag[0]), 1e-9 * (-sa * mag[1] + ca * mag[0])};
}
}  // namespace IGRF
