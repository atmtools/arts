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

Numeric field_interp(Numeric x, Numeric, Numeric, Numeric) noexcept {
  return x;
}

Numeric field_interp(const Tensor3&, Numeric, Numeric, Numeric){
    ARTS_USER_ERROR("Cannot field interp a tensor4")}

Numeric field_interp(const FunctionalData& f,
                     Numeric pre_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  return f(pre_point, lat_point, lon_point);
}

Numeric field_interp(const GriddedField3& x,
                     Numeric pre_point,
                     Numeric lat_point,
                     Numeric lon_point) {
  const auto pre_lagr =
      Interpolation::FixedLagrange<0>(0, pre_point, x.get_numeric_grid(0), false, Interpolation::GridType::Log);
  const auto lat_lagr =
      Interpolation::FixedLagrange<0>(0, lat_point, x.get_numeric_grid(1));
  const auto lon_lagr =
      Interpolation::FixedLagrange<0>(0, lon_point, x.get_numeric_grid(2));

  const auto iww =
      Interpolation::interpweights(pre_lagr, lat_lagr, lon_lagr);

  return Interpolation::interp(
      x.data, iww, pre_lagr, lat_lagr, lon_lagr);
}

Point Field::at(Numeric pre_point,
                Numeric lat_point,
                Numeric lon_point) const {
  Point atm{Key::pressure, pre_point};

  if (not regularized) {
    const auto get_value = [&](auto& xs) {
      return std::visit(
          [&](auto& x) {
            return field_interp(x, pre_point, lat_point, lon_point);
          },
          xs);
    };

    for (auto& vals : specs) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : other) atm.set(vals.first, get_value(vals.second));
    for (auto& vals : nlte) atm.set(vals.first, get_value(vals.second));
  } else {
    const auto pre_lagr = Interpolation::FixedLagrange<0>(0, pre_point, pre, false, Interpolation::GridType::Log);
    const auto lat_lagr = Interpolation::FixedLagrange<0>(0, lat_point, lat);
    const auto lon_lagr = Interpolation::FixedLagrange<0>(0, lon_point, lon);
    const auto iww =
        Interpolation::interpweights(pre_lagr, lat_lagr, lon_lagr);

    const auto get_value = [&](auto& x) {
      return Interpolation::interp(*std::get_if<Tensor3>(&x),
                                   iww,
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
Field& Field::regularize(const Vector& pressures,
                         const Vector& latitudes,
                         const Vector& longitudes) {
  ARTS_USER_ERROR_IF(regularized, "Cannot re-regularize a regularized grid")
  ArrayOfTensor3 specs_data(
      specs.size(),
      Tensor3(pressures.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 other_data(
      other.size(),
      Tensor3(pressures.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 nlte_data(
      nlte.size(),
      Tensor3(pressures.size(), latitudes.size(), longitudes.size()));

#pragma omp parallel for collapse(3) if (!arts_omp_in_parallel())
  for (Index i2 = 0; i2 < pressures.size(); i2++) {
    for (Index i3 = 0; i3 < latitudes.size(); i3++) {
      for (Index i4 = 0; i4 < longitudes.size(); i4++) {
        const auto pnt = at(pressures[i2], latitudes[i3], longitudes[i4]);

        for (Index i0{0}; auto &vals : specs)
          specs_data[i0++](i2, i3, i4) = pnt[vals.first];
        for (Index i0{0}; auto &vals : other)
          other_data[i0++](i2, i3, i4) = pnt[vals.first];
        for (Index i0{0}; auto &vals : nlte)
          nlte_data[i0++](i2, i3, i4) = pnt[vals.first];
      }
    }
  }

  regularized = true;
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

std::pair<std::array<Index, 3>, std::array<Index, 3>> shape(
    const GriddedField& gf) {
  std::array<Index, 3> size{1, 1, 1};
  std::array<Index, 3> pos{ -1, -1, -1};

  for (Index i = 0; i < gf.get_dim(); i++) {
    if (const auto& name = gf.get_grid_name(i); name == "Pressure") {
      size[0] = gf.get_grid_size(i);
      pos[0] = i;
    } else if (name == "Latitude") {
      size[1] = gf.get_grid_size(i);
      pos[1] = i;
    } else if (name == "Longitude") {
      size[2] = gf.get_grid_size(i);
      pos[2] = i;
    } else {
      ARTS_USER_ERROR("Bad grid name: ", name)
    }
  }

  ARTS_USER_ERROR_IF(std::count(pos.begin(), pos.end(), -1) not_eq gf.get_dim(), "Bad griddedfield:\n", gf)

  return {size, pos};
}

auto get_value(const auto& in_data,
               const std::array<Index, 3>& ind,
               const std::array<Index, 3>& pos) {
  using T = decltype(in_data);
  if constexpr (matpack::has_dim<T, 3>) {
    return in_data(ind[pos[0]], ind[pos[1]], ind[pos[2]]);
  } else if constexpr (matpack::has_dim<T, 2>) {
    if (pos[0] == -1) return in_data(ind[pos[1]], ind[pos[2]]);
    if (pos[1] == -1) return in_data(ind[pos[0]], ind[pos[2]]);
    return in_data(ind[pos[0]], ind[pos[1]]);
  } else {
    return in_data[ind[*std::find_if_not(pos.begin(), pos.end(), ne(-1))]];
  }
}

void adapt_data(Tensor3& out_data,
                const auto& in_data,
                const GriddedField& gf) {
  auto [size, pos] = shape(gf);

  out_data = Tensor3(size[0], size[1], size[2]);
  for (Index i = 0; i < size[0]; i++) {
    for (Index j = 0; j < size[1]; j++) {
      for (Index k = 0; k < size[2]; k++) {
        out_data(i, j, k) = get_value(in_data, {i, j, k}, pos);
      }
    }
  }
}

void adapt_grid(GriddedField& out, const GriddedField& in) {
  out.set_grid_name(0, "Pressure");
  out.set_grid_name(1, "Latitude");
  out.set_grid_name(2, "Longitude");

  out.set_grid(0, Vector(1, 0));
  out.set_grid(1, Vector(1, 0));
  out.set_grid(2, Vector(1, 0));

  for (Index i = 0; i < in.get_dim(); i++) {
    if (auto& name = in.get_grid_name(i); name == "Pressure") {
      out.set_grid(0, in.get_numeric_grid(i));
    } else if (name == "Latitude") {
      out.set_grid(1, in.get_numeric_grid(i));
    } else if (name == "Longitude") {
      out.set_grid(2, in.get_numeric_grid(i));
    } else {
      ARTS_USER_ERROR("Bad grid name: ", name)
    }
  }
}

GriddedField3 adapt(const auto& in_data, const GriddedField& gf) {
  GriddedField3 out(gf.get_name());
  adapt_data(out.data, in_data, gf);
  adapt_grid(out, gf);
  return out;
}
}  // namespace internal

GriddedField3 fix(const GriddedField3& gf) {
  return internal::adapt(gf.data, gf);
}

GriddedField3 fix(const GriddedField2& gf) {
  return internal::adapt(gf.data, gf);
}

GriddedField3 fix(const GriddedField1& gf) {
  return internal::adapt(gf.data, gf);
}
}  // namespace Atm
