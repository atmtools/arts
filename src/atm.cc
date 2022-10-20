#include "atm.h"

#include <algorithm>
#include <ostream>
#include <variant>

#include "debug.h"
#include "gridded_fields.h"
#include "interpolation_lagrange.h"
#include "matpackIV.h"

namespace Atm {
std::ostream& operator<<(std::ostream& os, const Point& atm) {
  os << "Temperature: " << atm.temperature << " K,\n";
  os << "Pressure: " << atm.pressure << " Pa,\n";
  os << "Wind Field: [u: " << atm.wind[0] << ", v: " << atm.wind[1]
     << ", w: " << atm.wind[2] << "] m/s,\n";
  os << "Magnetic Field: [u: " << atm.mag[0] << ", v: " << atm.mag[1]
     << ", w: " << atm.mag[2] << "] T";

  for (auto& spec : atm.species_content) {
    os << ",\n" << spec.first << ": " << spec.second;
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const Field& atm) {
  const auto printer = [&](auto& data) {
    if constexpr (isFunctionalDataType<decltype(data)>)
      os << "Functional\n";
    else
      os << data;
  };
  os << "Regularized Field: " << (atm.regularized ? "true " : "false");
  if (atm.regularized) {
    os << ",\nTime: [" << atm.time << "]";
    os << ",\nAltitude: [" << atm.alt << "]";
    os << ",\nLatitude: [" << atm.lat << "]";
    os << ",\nLongitude: [" << atm.lon << "]";
  }

  for (auto& vals : atm.other) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second);
  }

  for (auto& vals : atm.specs) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second);
  }

  return os;
}

Numeric field_interp(Numeric x, Time, Numeric, Numeric, Numeric) noexcept {
  return x;
}

Numeric field_interp(const Tensor4&, Time, Numeric, Numeric, Numeric){
    ARTS_USER_ERROR("Cannot field interp a tensor4")}

Numeric field_interp(const FunctionalData& f,
                     Time time_point,
                     Numeric alt_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  return f(time_point, alt_point, lat_point, lon_point);
}

Numeric field_interp(const GriddedField4& x,
                     Time time_point,
                     Numeric alt_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  const auto time_lagr =
      Interpolation::FixedLagrange<0>(0, time_point, x.get_time_grid(0));
  const auto alt_lagr =
      Interpolation::FixedLagrange<0>(0, alt_point, x.get_numeric_grid(1));
  const auto lat_lagr =
      Interpolation::FixedLagrange<0>(0, lat_point, x.get_numeric_grid(2));
  const auto lon_lagr =
      Interpolation::FixedLagrange<0>(0, lon_point, x.get_numeric_grid(3));

  const auto iww =
      Interpolation::interpweights(time_lagr, alt_lagr, lat_lagr, lon_lagr);

  return Interpolation::interp(
      x.data, iww, time_lagr, alt_lagr, lat_lagr, lon_lagr);
}

Point Field::at(Time time_point,
                Numeric alt_point,
                Numeric lat_point,
                Numeric lon_point) const {
  Point atm;

  if (not regularized) {
    const auto get_value = [&](auto& xs) {
      return std::visit(
          [&](auto& x) {
            return field_interp(x, time_point, alt_point, lat_point, lon_point);
          },
          xs);
    };

    for (auto& vals : specs) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : other) atm.set(vals.first, get_value(vals.second));
  } else {
    const auto time_lagr = Interpolation::FixedLagrange<0>(0, time_point, time);
    const auto alt_lagr = Interpolation::FixedLagrange<0>(0, alt_point, alt);
    const auto lat_lagr = Interpolation::FixedLagrange<0>(0, lat_point, lat);
    const auto lon_lagr = Interpolation::FixedLagrange<0>(0, lon_point, lon);
    const auto iww =
        Interpolation::interpweights(time_lagr, alt_lagr, lat_lagr, lon_lagr);

    const auto get_value = [&](auto& x) {
      return Interpolation::interp(*std::get_if<Tensor4>(&x),
                                   iww,
                                   time_lagr,
                                   alt_lagr,
                                   lat_lagr,
                                   lon_lagr);
    };

    for (auto& vals : specs) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : other) atm.set(vals.first, get_value(vals.second));
  }

  return atm;
}

//! Regularizes the calculations so that all data is on a single grid
Field& Field::regularize(const ArrayOfTime& times,
                         const Vector& altitudes,
                         const Vector& latitudes,
                         const Vector& longitudes) {
  ARTS_USER_ERROR_IF(regularized, "Cannot re-regularize a regularized grid")
  ArrayOfTensor4 specs_data(
      specs.size(),
      Tensor4(
          times.size(), altitudes.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor4 other_data(
      other.size(),
      Tensor4(
          times.size(), altitudes.size(), latitudes.size(), longitudes.size()));

#pragma omp parallel for collapse(4) if (!arts_omp_in_parallel())
  for (std::size_t i1 = 0; i1 < times.size(); i1++) {
    for (Index i2 = 0; i2 < altitudes.size(); i2++) {
      for (Index i3 = 0; i3 < latitudes.size(); i3++) {
        for (Index i4 = 0; i4 < longitudes.size(); i4++) {
          const auto pnt =
              at(times[i1], altitudes[i2], latitudes[i3], longitudes[i4]);

          for (Index i0{0}; auto& vals : specs)
            specs_data[i0++](i1, i2, i3, i4) = pnt[vals.first];
          for (Index i0{0}; auto& vals : other)
            other_data[i0++](i1, i2, i3, i4) = pnt[vals.first];
        }
      }
    }
  }

  regularized = true;
  time = times;
  alt = altitudes;
  lat = latitudes;
  lon = longitudes;

  for (Index i0{0}; auto& vals : specs)
    set(vals.first, std::move(specs_data[i0++]));
  for (Index i0{0}; auto& vals : other)
    set(vals.first, std::move(other_data[i0++]));

  return *this;
}
}  // namespace Atm
