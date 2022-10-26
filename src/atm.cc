#include "atm.h"

#include <algorithm>
#include <numeric>
#include <ostream>
#include <variant>

#include "arts_omp.h"
#include "compare.h"
#include "debug.h"
#include "gridded_fields.h"
#include "interpolation_lagrange.h"
#include "matpack.h"
#include "matpackIV.h"
#include "matpack_concepts.h"

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

  for (auto& vals : atm.nlte) {
    if (atm.nlte_energy and atm.nlte_energy->size())
      os << ",\n"
         << vals.first << ": " << vals.second << " K; "
         << atm.nlte_energy->at(vals.first) << " J";
    else
      os << ",\n" << vals.first << ": " << vals.second;
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
    os << ",\nPressure: [" << atm.pre << "]";
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

  for (auto& vals : atm.nlte) {
    os << ",\n" << vals.first << ":";
    if (atm.nlte_energy and atm.nlte_energy->size())
      os << " " << atm.nlte_energy->at(vals.first) << " J\n";
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
                     Numeric pre_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  return f(time_point, pre_point, lat_point, lon_point);
}

Numeric field_interp(const GriddedField4& x,
                     Time time_point,
                     Numeric pre_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  const auto time_lagr =
      Interpolation::FixedLagrange<0>(0, time_point, x.get_time_grid(0));
  const auto pre_lagr =
      Interpolation::FixedLagrange<0>(0, pre_point, x.get_numeric_grid(1), false, Interpolation::GridType::Log);
  const auto lat_lagr =
      Interpolation::FixedLagrange<0>(0, lat_point, x.get_numeric_grid(2));
  const auto lon_lagr =
      Interpolation::FixedLagrange<0>(0, lon_point, x.get_numeric_grid(3));

  const auto iww =
      Interpolation::interpweights(time_lagr, pre_lagr, lat_lagr, lon_lagr);

  return Interpolation::interp(
      x.data, iww, time_lagr, pre_lagr, lat_lagr, lon_lagr);
}

Point Field::at(Time time_point,
                Numeric pre_point,
                Numeric lat_point,
                Numeric lon_point) const {
  Point atm{Key::pressure, pre_point};

  if (not regularized) {
    const auto get_value = [&](auto& xs) {
      return std::visit(
          [&](auto& x) {
            return field_interp(x, time_point, pre_point, lat_point, lon_point);
          },
          xs);
    };

    for (auto& vals : specs) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : other) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : nlte) atm.set(vals.first, get_value(vals.second));
  } else {
    const auto time_lagr = Interpolation::FixedLagrange<0>(0, time_point, time);
    const auto pre_lagr = Interpolation::FixedLagrange<0>(0, pre_point, pre, false, Interpolation::GridType::Log);
    const auto lat_lagr = Interpolation::FixedLagrange<0>(0, lat_point, lat);
    const auto lon_lagr = Interpolation::FixedLagrange<0>(0, lon_point, lon);
    const auto iww =
        Interpolation::interpweights(time_lagr, pre_lagr, lat_lagr, lon_lagr);

    const auto get_value = [&](auto& x) {
      return Interpolation::interp(*std::get_if<Tensor4>(&x),
                                   iww,
                                   time_lagr,
                                   pre_lagr,
                                   lat_lagr,
                                   lon_lagr);
    };

    for (auto& vals : specs) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : other) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : nlte) atm.set(vals.first, get_value(vals.second));
  }

  atm.set_energy_level(nlte_energy);
  return atm;
}

//! Regularizes the calculations so that all data is on a single grid
Field& Field::regularize(const ArrayOfTime& times,
                         const Vector& pressures,
                         const Vector& latitudes,
                         const Vector& longitudes) {
  ARTS_USER_ERROR_IF(regularized, "Cannot re-regularize a regularized grid")
  ArrayOfTensor4 specs_data(
      specs.size(),
      Tensor4(
          times.size(), pressures.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor4 other_data(
      other.size(),
      Tensor4(
          times.size(), pressures.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor4 nlte_data(
      nlte.size(),
      Tensor4(
          times.size(), pressures.size(), latitudes.size(), longitudes.size()));

#pragma omp parallel for collapse(4) if (!arts_omp_in_parallel())
  for (std::size_t i1 = 0; i1 < times.size(); i1++) {
    for (Index i2 = 0; i2 < pressures.size(); i2++) {
      for (Index i3 = 0; i3 < latitudes.size(); i3++) {
        for (Index i4 = 0; i4 < longitudes.size(); i4++) {
          const auto pnt =
              at(times[i1], pressures[i2], latitudes[i3], longitudes[i4]);

          for (Index i0{0}; auto& vals : specs)
            specs_data[i0++](i1, i2, i3, i4) = pnt[vals.first];
          for (Index i0{0}; auto& vals : other)
            other_data[i0++](i1, i2, i3, i4) = pnt[vals.first];
          for (Index i0{0}; auto& vals : nlte)
            nlte_data[i0++](i1, i2, i3, i4) = pnt[vals.first];
        }
      }
    }
  }

  regularized = true;
  time = times;
  pre = pressures;
  lat = latitudes;
  lon = longitudes;

  for (Index i0{0}; auto& vals : specs)
    set(vals.first, std::move(specs_data[i0++]));
  for (Index i0{0}; auto& vals : other)
    set(vals.first, std::move(other_data[i0++]));
  for (Index i0{0}; auto& vals : nlte)
    set(vals.first, std::move(nlte_data[i0++]));

  return *this;
}

namespace internal {
  using namespace Cmp;

std::pair<std::array<Index, 4>, std::array<Index, 4>> shape(
    const GriddedField& gf) {
  std::array<Index, 4> size{1, 1, 1, 1};
  std::array<Index, 4> pos{-1, -1, -1, -1};

  for (Index i = 0; i < gf.get_dim(); i++) {
    if (auto& name = gf.get_grid_name(i); name == "Time") {
      size[0] = gf.get_grid_size(i);
      pos[0] = i;
    } else if (name == "Pressure") {
      size[1] = gf.get_grid_size(i);
      pos[1] = i;
    } else if (name == "Latitude") {
      size[2] = gf.get_grid_size(i);
      pos[2] = i;
    } else if (name == "Longitude") {
      size[3] = gf.get_grid_size(i);
      pos[3] = i;
    } else {
      ARTS_USER_ERROR("Bad grid name: ", name)
    }
  }

  ARTS_USER_ERROR_IF(std::count(pos.begin(), pos.end(), -1) not_eq gf.get_dim(), "Bad griddedfield:\n", gf)

  return {size, pos};
}

auto get_value(const auto& in_data,
               const std::array<Index, 4>& ind,
               const std::array<Index, 4>& pos) {
  using T = decltype(in_data);
  if constexpr (matpack::has_dim<T, 4>) {
    return in_data(ind[pos[0]], ind[pos[1]], ind[pos[2]], ind[pos[3]]);
  } else if constexpr (matpack::has_dim<T, 3>) {
    if (pos[0] == -1) return in_data(ind[pos[1]], ind[pos[2]], ind[pos[3]]);
    if (pos[1] == -1) return in_data(ind[pos[0]], ind[pos[2]], ind[pos[3]]);
    if (pos[2] == -1) return in_data(ind[pos[0]], ind[pos[1]], ind[pos[3]]);
    return in_data(ind[pos[0]], ind[pos[1]], ind[pos[2]]);
  } else if constexpr (matpack::has_dim<T, 2>) {
    if (pos[0] == -1) {
      if (pos[1] == -1) return in_data(ind[pos[2]], ind[pos[3]]);
      if (pos[2] == -1) return in_data(ind[pos[1]], ind[pos[3]]);
      return in_data(ind[pos[1]], ind[pos[2]]);
    }
    if (pos[1] == -1) {
      if (pos[2] == -1) return in_data(ind[pos[0]], ind[pos[3]]);
      return in_data(ind[pos[0]], ind[pos[2]]);
    }
    return in_data(ind[pos[0]], ind[pos[1]]);
  } else {
    return in_data[ind[*std::find_if_not(pos.begin(), pos.end(), ne(-1))]];
  }
}

void adapt_data(Tensor4& out_data,
                const auto& in_data,
                const GriddedField& gf) {
  auto [size, pos] = shape(gf);

  out_data = Tensor4(size[0], size[1], size[2], size[3]);
  for (Index i = 0; i < size[0]; i++) {
    for (Index j = 0; j < size[1]; j++) {
      for (Index k = 0; k < size[2]; k++) {
        for (Index m = 0; m < size[3]; m++) {
          out_data(i, j, k, m) = get_value(in_data, {i, j, k, m}, pos);
        }
      }
    }
  }
}

void adapt_grid(GriddedField& out, const GriddedField& in) {
  out.set_grid_name(0, "Time");
  out.set_grid_name(1, "Pressure");
  out.set_grid_name(2, "Latitude");
  out.set_grid_name(3, "Longitude");

  out.set_grid(0, ArrayOfTime(1, Time{}));
  out.set_grid(1, Vector(1, 0));
  out.set_grid(2, Vector(1, 0));
  out.set_grid(3, Vector(1, 0));

  for (Index i = 0; i < in.get_dim(); i++) {
    if (auto& name = in.get_grid_name(i); name == "Time") {
      out.set_grid(0, in.get_time_grid(i));
    } else if (name == "Pressure") {
      out.set_grid(0, in.get_numeric_grid(i));
    } else if (name == "Latitude") {
      out.set_grid(0, in.get_numeric_grid(i));
    } else if (name == "Longitude") {
      out.set_grid(0, in.get_numeric_grid(i));
    } else {
      ARTS_USER_ERROR("Bad grid name: ", name)
    }
  }
}

GriddedField4 adapt(const auto& in_data, const GriddedField& gf) {
  GriddedField4 out(gf.get_name());
  adapt_data(out.data, in_data, gf);
  adapt_grid(out, gf);
  return out;
}
}  // namespace internal

GriddedField4 fix(const GriddedField4& gf) {
  return internal::adapt(gf.data, gf);
}

GriddedField4 fix(const GriddedField3& gf) {
  return internal::adapt(gf.data, gf);
}

GriddedField4 fix(const GriddedField2& gf) {
  return internal::adapt(gf.data, gf);
}

GriddedField4 fix(const GriddedField1& gf) {
  return internal::adapt(gf.data, gf);
}
}  // namespace Atm
