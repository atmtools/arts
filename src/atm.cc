#include "atm.h"

#include <algorithm>
#include <exception>
#include <iomanip>
#include <limits>
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

  for (auto& spec : atm.specs) {
    os << ",\n" << spec.first << ": " << spec.second;
  }

  for (auto& vals : atm.nlte) {
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
    os << ",\nAltitude: [" << atm.grid[0] << "]";
    os << ",\nLatitude: [" << atm.grid[1] << "]";
    os << ",\nLongitude: [" << atm.grid[2] << "]";
  }

  for (auto& vals : atm.other) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second.data);
  }

  for (auto& vals : atm.specs) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second.data);
  }

  for (auto& vals : atm.nlte) {
    os << ",\n" << vals.first << ":";
    std::visit(printer, vals.second.data);
  }

  return os;
}

Numeric Point::mean_mass() const {
  Numeric out = 0;
  Numeric total_vmr = 0;
  for (auto& x: specs) {
    Numeric this_mass_sum = 0.0;

    for (auto& spec: x.first) {
      this_mass_sum += spec.Mass();
    }

    if (x.first.nelem()) this_mass_sum *= x.second / static_cast<Numeric>(x.first.nelem());
    out += this_mass_sum;
    total_vmr += x.second;
  }

  if (total_vmr not_eq 0) out /= total_vmr;
  return out;
}

std::vector<KeyVal> Point::keys() const {
  std::vector<KeyVal> out;
  out.reserve(nelem());
  for (auto &a : enumtyps::KeyTypes)
    out.emplace_back(a);
  for (auto &a : specs)
    out.emplace_back(a.first);
  for (auto &a : nlte)
    out.emplace_back(a.first);
  return out;
}

std::vector<KeyVal> Field::keys() const {
  std::vector<KeyVal> out;
  out.reserve(nelem());
  for (auto &a : other)
    out.emplace_back(a.first);
  for (auto &a : specs)
    out.emplace_back(a.first);
  for (auto &a : nlte)
    out.emplace_back(a.first);
  return out;
}

Index Point::nspec() const { return static_cast<Index>(specs.size()); }
Index Point::nnlte() const { return static_cast<Index>(nlte.size()); }
Index Point::nelem() const { return nspec() + nnlte() + nother(); }

Index Field::nspec() const { return static_cast<Index>(specs.size()); }
Index Field::nnlte() const { return static_cast<Index>(nlte.size()); }
Index Field::nother() const { return static_cast<Index>(other.size()); }
Index Field::nelem() const { return nspec() + nnlte() + nother(); }

String Data::data_type() const {
  if (std::holds_alternative<GriddedField3>(data)) return "GriddedField3";
  if (std::holds_alternative<Tensor3>(data)) return "Tensor3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";
  ARTS_ASSERT(false, "Cannot be reached, you have added a new type but not doen the plumbing...")
  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

void Field::throwing_check() const {
  if (regularized) {
    for (auto &key_val : other) {
      std::visit(
          [&](auto &&x) {
            if constexpr (not isTensor3<decltype(x)>) {
              ARTS_USER_ERROR(
                  "The data for ", key_val.first,
                  " is not a Tensor3 even though the data is regularized")
            } else {
              ARTS_USER_ERROR_IF(x.npages() not_eq grid[0].nelem() or
                                     x.nrows() not_eq grid[1].nelem() or
                                     x.ncols() not_eq grid[2].nelem(),
                                 "Mismatch dimensions.  Expects shape:\n[",
                                 grid[0].nelem(), ", ", grid[1].nelem(), ", ",
                                 grid[2].nelem(), "]\nGot shape:\n[", x.npages(),
                                 ", ", grid[1].nelem(), ", ", grid[2].nelem(),
                                 "]\nFor field: ", key_val.first)
            }
          },
          key_val.second.data);
    }

    for (auto &key_val : specs) {
      std::visit(
          [&](auto &&x) {
            if constexpr (not isTensor3<decltype(x)>) {
              ARTS_USER_ERROR(
                  "The data for ", key_val.first,
                  " is not a Tensor3 even though the data is regularized")
            } else {
              ARTS_USER_ERROR_IF(x.npages() not_eq grid[0].nelem() or
                                     x.nrows() not_eq grid[1].nelem() or
                                     x.ncols() not_eq grid[2].nelem(),
                                 "Mismatch dimensions.  Expects shape:\n[",
                                 grid[0].nelem(), ", ", grid[1].nelem(), ", ",
                                 grid[2].nelem(), "]\nGot shape:\n[", x.npages(),
                                 ", ", grid[1].nelem(), ", ", grid[2].nelem(),
                                 "]\nFor field: ", key_val.first)
            }
          },
          key_val.second.data);
    }

    for (auto &key_val : nlte) {
      std::visit(
          [&](auto &&x) {
            if constexpr (not isTensor3<decltype(x)>) {
              ARTS_USER_ERROR(
                  "The data for ", key_val.first,
                  " is not a Tensor3 even though the data is regularized")
            } else {
              ARTS_USER_ERROR_IF(x.npages() not_eq grid[0].nelem() or
                                     x.nrows() not_eq grid[1].nelem() or
                                     x.ncols() not_eq grid[2].nelem(),
                                 "Mismatch dimensions.  Expects shape:\n[",
                                 grid[0].nelem(), ", ", grid[1].nelem(), ", ",
                                 grid[2].nelem(), "]\nGot shape:\n[", x.npages(),
                                 ", ", grid[1].nelem(), ", ", grid[2].nelem(),
                                 "]\nFor field: ", key_val.first)
            }
          },
          key_val.second.data);
    }
  }
}

namespace detail {
constexpr Numeric field_interp(Numeric x, Numeric, Numeric, Numeric) noexcept {
  return x;
}

Numeric field_interp(const Tensor3 &, Numeric, Numeric, Numeric) {
  ARTS_USER_ERROR("Cannot field interp a Tensor3")
}

Numeric field_interp(const FunctionalData &f, Numeric alt_point,
                     Numeric lat_point, Numeric lon_point) {
  return f(alt_point, lat_point, lon_point);
}

template <std::size_t NALT, std::size_t NLAT, std::size_t NLON>
struct InterpCapture {
  Interpolation::FixedLagrange<NALT> alt;
  Interpolation::FixedLagrange<NLAT> lat;
  Interpolation::FixedLagrange<NLON> lon;
  FixedGrid<Numeric, NALT + 1, NLAT + 1, NLON + 1> iw;

  InterpCapture(const Vector &alt_grid, const Vector &lat_grid,
                const Vector &lon_grid, Numeric alt_point, Numeric lat_point,
                Numeric lon_point)
      : alt(0, alt_point, alt_grid), lat(0, lat_point, lat_grid),
        lon(0, lon_point, lon_grid, false, Interpolation::GridType::Cyclic,
            {-180, 180}),
        iw(Interpolation::interpweights(alt, lat, lon)) {}
};

using LinearInterpolation = detail::InterpCapture<1, 1, 1>;

Numeric field_interp(const GriddedField3 &x, Numeric alt_point,
                     Numeric lat_point, Numeric lon_point) {
  auto [alt, lat, lon, iw] = LinearInterpolation(
      x.get_numeric_grid(0), x.get_numeric_grid(1), x.get_numeric_grid(2),
      alt_point, lat_point, lon_point);

  return Interpolation::interp(x.data, iw, alt, lat, lon);
}

bool limits(Numeric &x, const Vector &r, Extrapolation low, Extrapolation upp, const char* type) {
  ARTS_USER_ERROR_IF(r.size() == 0, "The ", std::quoted(type), " grid is empty")

  using enum Extrapolation;
  switch (low) {
  case Zero:
    if (x < r.front()) {
      x = 0;
      return true;
    }
    break;
  case Hydrostatic: [[fallthrough]];
  case Nearest:
    x = std::max(x, r.front());
    break;
  case None:
    ARTS_USER_ERROR_IF(x < r.front(), "The ", std::quoted(type), " value ", x,
                       " is below the range of ", r)
    break;
  case Linear:
    break;
  case FINAL: {
  }
  }

  switch (upp) {
  case Zero:
    if (x > r.back()) {
      x = 0;
      return true;
    }
    break;
  case Hydrostatic: [[fallthrough]];
  case Nearest:
    x = std::min(x, r.back());
    break;
  case None:
    ARTS_USER_ERROR_IF(x > r.back(), "The ", std::quoted(type), " value ", x, " is below the range of ",
                       r)
    break;
  case Linear:
    break;
  case FINAL: {
  }
  }

  return false;
}

Numeric compute_value(const Data &data, Numeric alt, Numeric lat, Numeric lon) {
  auto compute = [alt, lat, lon, alt_low = data.alt_low, alt_upp = data.alt_upp,
                  lat_low = data.lat_low, lat_upp = data.lat_upp,
                  lon_low = data.lon_low,
                  lon_upp = data.lon_upp](auto &&x) -> Numeric {
    using T = decltype(x);

    if constexpr (isGriddedField3<T>) {
      const Vector &alts = x.get_numeric_grid(0);
      const Vector &lats = x.get_numeric_grid(1);
      const Vector &lons = x.get_numeric_grid(2);

      Numeric test;
      if (limits(test=alt, alts, alt_low, alt_upp, "altitude")) return test;
      if (limits(test=lat, lats, lat_low, lat_upp, "latitude")) return test;
      if (limits(test=lon, lons, lon_low, lon_upp, "longitude")) return test;
    }

    return field_interp(x, alt, lat, lon);
  };

  return std::visit(compute, data.data);
}

Numeric fix_hydrostatic(Numeric x0, const Point& atm, const Data& data, Numeric g, Numeric alt) {
  const auto need_fixing = [&](auto& x) {
    using T = decltype(x);

    if constexpr (isGriddedField3<T>) {
      const Vector& alts =  x.get_numeric_grid(0);
      if (data.alt_low == Extrapolation::Hydrostatic and alt < alts.front())
      return SimpleHydrostaticExpansion{x0, alts.front(), atm.temperature, atm.mean_mass(), g}(alt);

      if (data.alt_upp == Extrapolation::Hydrostatic and alt > alts.back())
      return SimpleHydrostaticExpansion{x0, alts.front(), atm.temperature, atm.mean_mass(), g}(alt);
    }

    return x0;
  };

  return std::visit(need_fixing, data.data);
}
} // namespace detail

std::tuple<Numeric, Numeric, Numeric>
Data::hydrostatic_coord(Numeric alt, Numeric lat, Numeric lon) const {
  if (const auto *ptr = std::get_if<GriddedField3>(&data); ptr) {
    detail::limits(alt, ptr->get_numeric_grid(0), alt_low, alt_upp, "altitude");
    detail::limits(lat, ptr->get_numeric_grid(1), lat_low, lat_upp, "latitude");
    detail::limits(lon, ptr->get_numeric_grid(2), lon_low, lon_upp,
                   "longitude");
  }
  return {alt, lat, lon};
}

Point Field::internal_fitting(Numeric alt_point, Numeric lat_point, Numeric lon_point) const {
  Point atm;

  if (alt_point > top_of_atmosphere)
    return atm;

  if (not regularized) {
    const auto get_value = [alt = alt_point, lat = lat_point,
                            lon = lon_point](auto &&x) {
      return detail::compute_value(x, alt, lat, lon);
    };

    for (auto &vals : specs)
      atm[vals.first] = get_value(vals.second);
    for (auto &vals : other)
      atm[vals.first] = get_value(vals.second);
    for (auto &vals : nlte)
      atm[vals.first] = get_value(vals.second);
  } else {
    const auto get_value = [v = detail::LinearInterpolation(
                                grid[0], grid[1], grid[2], alt_point,
                                lat_point, lon_point)](auto &x) {
      return Interpolation::interp(*std::get_if<Tensor3>(&x.data), v.iw, v.alt,
                                   v.lat, v.lon);
    };

    for (auto &vals : specs)
      atm[vals.first] = get_value(vals.second);
    for (auto &vals : other)
      atm[vals.first] = get_value(vals.second);
    for (auto &vals : nlte)
      atm[vals.first] = get_value(vals.second);
  }
  
  return atm;
}

Point Field::at(Numeric alt_point, Numeric lat_point, Numeric lon_point,
                const FunctionalData &g) const {
  throwing_check();

  Point atm = internal_fitting(alt_point, lat_point, lon_point);

  // Fix for hydrostatic equilibrium
  for (auto &vals : specs) {
    if (vals.second.need_hydrostatic()) {
      const auto [base_alt, base_lat, base_lon] =
          vals.second.hydrostatic_coord(alt_point, lat_point, lon_point);
      if (base_alt not_eq alt_point) {
        auto base_atm = internal_fitting(base_alt, base_lat, base_lon);
        auto base_gra = g(base_alt, base_lat, base_lon);
        atm[vals.first] =
                detail::fix_hydrostatic(base_atm[vals.first], base_atm,
                                        vals.second, base_gra, alt_point);
      }
    }
  }

  for (auto &vals : other) {
    if (vals.second.need_hydrostatic()) {
      const auto [base_alt, base_lat, base_lon] =
          vals.second.hydrostatic_coord(alt_point, lat_point, lon_point);
      if (base_alt not_eq alt_point) {
        auto base_atm = internal_fitting(base_alt, base_lat, base_lon);
        auto base_gra = g(base_alt, base_lat, base_lon);
        atm[vals.first] =
                detail::fix_hydrostatic(base_atm[vals.first], base_atm,
                                        vals.second, base_gra, alt_point);
      }
    }
  }

  for (auto &vals : nlte) {
    if (vals.second.need_hydrostatic()) {
      const auto [base_alt, base_lat, base_lon] =
          vals.second.hydrostatic_coord(alt_point, lat_point, lon_point);
      if (base_alt not_eq alt_point) {
        auto base_atm = internal_fitting(base_alt, base_lat, base_lon);
        auto base_gra = g(base_alt, base_lat, base_lon);
        atm[vals.first] =
                detail::fix_hydrostatic(base_atm[vals.first], base_atm,
                                        vals.second, base_gra, alt_point);
      }
    }
  }

  return atm;
}

//! Regularizes the calculations so that all data is on a single grid
Field& Field::regularize(const Vector& altitudes,
                         const Vector& latitudes,
                         const Vector& longitudes) {
  ARTS_USER_ERROR_IF(regularized, "Cannot re-regularize a regularized grid")

  ArrayOfTensor3 specs_data(
      specs.size(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 other_data(
      other.size(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 nlte_data(
      nlte.size(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));

  for (Index i2 = 0; i2 < altitudes.size(); i2++) {
    for (Index i3 = 0; i3 < latitudes.size(); i3++) {
      for (Index i4 = 0; i4 < longitudes.size(); i4++) {
        const auto pnt = at(altitudes[i2], latitudes[i3], longitudes[i4]);

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
  grid = {altitudes, latitudes, longitudes};

  for (Index i0{0}; auto& vals : specs)
    specs[vals.first] = std::move(specs_data[i0++]);
  for (Index i0{0}; auto& vals : other)
    other[vals.first] = std::move(other_data[i0++]);
  for (Index i0{0}; auto& vals : nlte)
    nlte[vals.first] = std::move(nlte_data[i0++]);

  return *this;
}

std::array<Index, 3> Field::regularized_shape() const {
  return {grid[0].nelem(), grid[1].nelem(), grid[2].nelem()};
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
