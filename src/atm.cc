#include "atm.h"

#include <algorithm>
#include <exception>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ostream>
#include <type_traits>
#include <variant>
#include <vector>

#include "arts_omp.h"
#include "compare.h"
#include "debug.h"
#include "gridded_fields.h"
#include "grids.h"
#include "interp.h"
#include "interpolation.h"
#include "matpack_algo.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_iter.h"

namespace Atm {
const std::unordered_map<QuantumIdentifier, Data> &Field::nlte() const {
  return map<QuantumIdentifier>();
}
const std::unordered_map<ArrayOfSpeciesTag, Data> &Field::specs() const {
  return map<ArrayOfSpeciesTag>();
}
const std::unordered_map<Key, Data> &Field::other() const {
  return map<Key>();
}
const std::unordered_map<ParticulatePropertyTag, Data> &Field::partp() const {
  return map<ParticulatePropertyTag>();
}

std::unordered_map<QuantumIdentifier, Data> &Field::nlte() {
  return map<QuantumIdentifier>();
}
std::unordered_map<ArrayOfSpeciesTag, Data> &Field::specs() {
  return map<ArrayOfSpeciesTag>();
}
std::unordered_map<Key, Data> &Field::other() { return map<Key>(); }
std::unordered_map<ParticulatePropertyTag, Data> &Field::partp() {
  return map<ParticulatePropertyTag>();
}

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
  

  for (auto& vals : atm.other()) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second.data);
  }

  for (auto& vals : atm.specs()) {
    os << ",\n" << vals.first << ":\n";
    std::visit(printer, vals.second.data);
  }

  for (auto& vals : atm.nlte()) {
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
  for (auto &a : partp)
    out.emplace_back(a.first);
  return out;
}

Index Point::nspec() const { return static_cast<Index>(specs.size()); }
Index Point::npart() const { return static_cast<Index>(partp.size()); }
Index Point::nnlte() const { return static_cast<Index>(nlte.size()); }
Index Point::nelem() const { return nspec() + nnlte() + nother() + npart(); }

Index Field::nspec() const { return static_cast<Index>(specs().size()); }
Index Field::npart() const { return static_cast<Index>(partp().size()); }
Index Field::nnlte() const { return static_cast<Index>(nlte().size()); }
Index Field::nother() const { return static_cast<Index>(other().size()); }

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
    for (auto &key_val : other()) {
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

    for (auto &key_val : specs()) {
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

    for (auto &key_val : nlte()) {
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
struct Limits {
  Numeric alt_low{std::numeric_limits<Numeric>::lowest()};
  Numeric alt_upp{std::numeric_limits<Numeric>::max()};
  Numeric lat_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lat_upp{std::numeric_limits<Numeric>::max()};
  Numeric lon_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lon_upp{std::numeric_limits<Numeric>::max()};
};

struct ComputeLimit {
  Extrapolation type{Extrapolation::Linear};
  Numeric alt, lat, lon;
};

Limits find_limits(const Numeric&) {
  return {};
}

Limits find_limits(const FunctionalData&) {
  return {};
}

Limits find_limits(const GriddedField3 &gf3) {
  return {gf3.get_numeric_grid(0).front(), gf3.get_numeric_grid(0).back(),
          gf3.get_numeric_grid(1).front(), gf3.get_numeric_grid(1).back(),
          gf3.get_numeric_grid(2).front(), gf3.get_numeric_grid(2).back()};
}

Limits find_limits(const Tensor3 &) {
  ARTS_ASSERT(false, "This must be dealt with earlier");
  return {};
}

Vector vec_interp(const Numeric& v, const Vector& alt, const Vector&, const Vector&) {
  Vector out(alt.nelem(), v);
  return out;
}

Vector vec_interp(const FunctionalData& v, const Vector& alt, const Vector& lat, const Vector& lon) {
  const Index n=alt.nelem();
  Vector out(n);
  for (Index i=0; i<n; i++) out[i] = v(alt[i], lat[i], lon[i]);
  return out;
}

template <Index poly_alt, Index poly_lat, Index poly_lon,
          bool precompute = false>
struct interp_helper {
  using AltLag = my_interp::Lagrange<poly_alt>;
  using LatLag = my_interp::Lagrange<poly_lat>;
  using LonLag = std::conditional_t<
      poly_lon == 0, my_interp::Lagrange<0>,
      my_interp::Lagrange<poly_lon, false, my_interp::GridType::Cyclic,
                          my_interp::cycle_m180_p180>>;

  Array<AltLag> lags_alt;
  Array<LatLag> lags_lat;
  Array<LonLag> lags_lon;
  [[no_unique_address]] std::conditional_t<
      precompute,
      decltype(my_interp::flat_interpweights(lags_alt, lags_lat, lags_lon)),
      my_interp::Empty>
      iws{};

  constexpr interp_helper(const Vector &alt_grid, const Vector &lat_grid,
                          const Vector &lon_grid, const Vector &alt,
                          const Vector &lat, const Vector &lon)
    requires(not precompute)
      : lags_alt(
            my_interp::lagrange_interpolation_list<AltLag>(alt, alt_grid, -1)),
        lags_lat(
            my_interp::lagrange_interpolation_list<LatLag>(lat, lat_grid, -1)),
        lags_lon(my_interp::lagrange_interpolation_list<LonLag>(lon, lon_grid,
                                                                -1)) {}

  constexpr interp_helper(const Vector &alt_grid, const Vector &lat_grid,
                          const Vector &lon_grid, const Vector &alt,
                          const Vector &lat, const Vector &lon)
    requires(precompute)
      : lags_alt(
            my_interp::lagrange_interpolation_list<AltLag>(alt, alt_grid, -1)),
        lags_lat(
            my_interp::lagrange_interpolation_list<LatLag>(lat, lat_grid, -1)),
        lags_lon(
            my_interp::lagrange_interpolation_list<LonLag>(lon, lon_grid, -1)),
        iws(flat_interpweights(lags_alt, lags_lat, lags_lon)) {}

  Vector operator()(const Tensor3 &data) const {
    if constexpr (precompute)
      return flat_interp(data, iws, lags_alt, lags_lat, lags_lon);
    else
      return flat_interp(data, lags_alt, lags_lat, lags_lon);
  }
};

template <Index poly_alt, Index poly_lat, Index poly_lon>
Vector tvec_interp(const Tensor3 &v, const Vector &alt_grid,
                   const Vector &lat_grid, const Vector &lon_grid,
                   const Vector &alt, const Vector &lat, const Vector &lon) {
  const interp_helper<poly_alt, poly_lat, poly_lon> interpolater(
      alt_grid, lat_grid, lon_grid, alt, lat, lon);
  return interpolater(v);
}

Vector vec_interp(const GriddedField3& v, const Vector& alt, const Vector& lat, const Vector& lon) {

  ARTS_ASSERT(v.get_grid_size(0) > 0)
  ARTS_ASSERT(v.get_grid_size(1) > 0)
  ARTS_ASSERT(v.get_grid_size(2) > 0)

  const bool d1 = v.get_grid_size(0) == 1;
  const bool d2 = v.get_grid_size(1) == 1;
  const bool d3 = v.get_grid_size(2) == 1;

  const Index n=alt.nelem();
  
  if (d1 and d2 and d3) return Vector(n, v.data(0, 0, 0));
  if (d1 and d2)        return tvec_interp<0, 0, 1>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  if (d1 and d3)        return tvec_interp<0, 1, 0>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  if (d2 and d3)        return tvec_interp<1, 0, 0>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  if (d1)               return tvec_interp<0, 1, 1>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  if (d2)               return tvec_interp<1, 0, 1>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  if (d3)               return tvec_interp<1, 1, 0>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
  return tvec_interp<1, 1, 1>(v.data, v.get_numeric_grid(0), v.get_numeric_grid(1), v.get_numeric_grid(2), alt, lat, lon);
}

Vector vec_interp(const Tensor3 &, const Vector &, const Vector &,
                  const Vector &) {
  ARTS_ASSERT(false, "This must be dealt with earlier")
  return {};
}

Numeric limit(const Data &data, ComputeLimit lim, Numeric orig) {
  ARTS_ASSERT(lim.type not_eq Extrapolation::FINAL)

  ARTS_USER_ERROR_IF(lim.type == Extrapolation::None,
                     "Altitude limit breaced.  Position (", lim.alt, ", ",
                     lim.lat, ", ", lim.lon,
                     ") is out-of-bounds when no extrapolation is wanted")

  if (lim.type == Extrapolation::Zero)
    return 0;

  if (lim.type == Extrapolation::Nearest)
    return std::visit(
        [&](auto &d) {
          return vec_interp(d, {lim.alt}, {lim.lat}, {lim.lon})[0];
        },
        data.data);

  return orig;
}

constexpr Extrapolation combine(Extrapolation a, Extrapolation b) {
  using enum Extrapolation;
  switch(a) {
    case None: return None;
    case Zero: {
      switch(b) {
        case None: return None;
        case Zero: return Zero;
        case Nearest: return Zero;
        case Linear: return Zero;
        case FINAL: ARTS_ASSERT(false);
      }
    } ARTS_ASSERT(false); return FINAL;
    case Nearest: {
      switch (b) {
        case None: return None;
        case Zero: return Zero;
        case Nearest: return Nearest;
        case Linear: return Nearest;
        case FINAL: ARTS_ASSERT(false);
      }
    } ARTS_ASSERT(false); return FINAL;
    case Linear: return b;
    case FINAL: ARTS_ASSERT(false);
  }

  return FINAL;
}

constexpr Extrapolation combine(Extrapolation a, Extrapolation b, Extrapolation c) {
  return combine(combine(a, b), c);
}

void select(Extrapolation lowt, Extrapolation uppt, Numeric lowv, Numeric uppv, Numeric v, Numeric& outv, Extrapolation& outt) {
  if (v < lowv) {
    outt = lowt;
    if (outt == Extrapolation::Nearest) v = lowv;
  } else if (uppv < v) {
    outt = uppt;
    if (outt == Extrapolation::Nearest) v = uppv;
  }

  outv = v;
}

ComputeLimit find_limit(const Data& data, const Limits& lim, Numeric alt, Numeric lat, Numeric lon) {
  ComputeLimit out;
  Extrapolation a{Extrapolation::Linear}, b{Extrapolation::Linear}, c{Extrapolation::Linear};
  
  select(data.alt_low, data.alt_upp, lim.alt_low, lim.alt_upp, alt, out.alt, a);
  select(data.lat_low, data.lat_upp, lim.lat_low, lim.lat_upp, lat, out.lat, b);
  select(data.lon_low, data.lon_upp, lim.lon_low, lim.lon_upp, lon, out.lon, c);
  out.type = combine(a, b, c);

  return out;
}

Vector vec_interp(const Data& data, const Vector& alt, const Vector& lat, const Vector& lon) {
  const auto compute = [&](auto& d) {
    return vec_interp(d, alt, lat, lon);
  };

  // Perform the interpolation
  Vector out = std::visit(compute, data.data);

  // Fix the extrapolations for ZERO and NONE and NEAREST
  const auto lim =
      std::visit([](auto &d) { return find_limits(d); }, data.data);
  const Index n = alt.nelem();
  for (Index i=0; i<n; i++) {
    out[i] = limit(data, find_limit(data, lim, alt[i], lat[i], lon[i]), out[i]);
  }

  return out;
}

template <Index poly_alt, Index poly_lat, Index poly_lon>
void tensor_interpolator_fun(std::vector<Point>& out, const Field& f, const Vector& alt, const Vector& lat, const Vector& lon) {
  const auto n = out.size();
  const auto interpolator = interp_helper<poly_alt, poly_lat, poly_lon, true>(f.grid[0], f.grid[1], f.grid[2], alt, lat, lon);
  for (const auto& key: f.keys()) {
    const Vector f_vec = interpolator(f[key].get<const Tensor3&>());
    for (std::size_t i=0; i<n; i++) {
      out[i][key] = f_vec[i];
    }
  }
}

void tensor_interpolator(std::vector<Point>& out, const Field& f, const Vector& alt, const Vector& lat, const Vector& lon) {
  const bool d1 = f.grid[0].nelem() == 1;
  const bool d2 = f.grid[1].nelem() == 1;
  const bool d3 = f.grid[2].nelem() == 1;

  if (d1 and d2 and d3) tensor_interpolator_fun<0, 0, 0>(out, f, alt, lat, lon);
  else if (d1 and d2)   tensor_interpolator_fun<0, 0, 1>(out, f, alt, lat, lon);
  else if (d1 and d3)   tensor_interpolator_fun<0, 1, 0>(out, f, alt, lat, lon);
  else if (d2 and d3)   tensor_interpolator_fun<1, 0, 0>(out, f, alt, lat, lon);
  else if (d1)          tensor_interpolator_fun<0, 1, 1>(out, f, alt, lat, lon);
  else if (d2)          tensor_interpolator_fun<1, 0, 1>(out, f, alt, lat, lon);
  else if (d3)          tensor_interpolator_fun<1, 1, 0>(out, f, alt, lat, lon);
  else                  tensor_interpolator_fun<1, 1, 1>(out, f, alt, lat, lon);
}
}  // namespace detail

void Field::at(std::vector<Point>& out, const Vector& alt, const Vector& lat, const Vector& lon) const {
  throwing_check();
  ARTS_USER_ERROR_IF(
      std::any_of(alt.begin(), alt.end(), Cmp::gt(top_of_atmosphere)),
      "Cannot get values above the top of the atmosphere, which is at: ",
      top_of_atmosphere, " m.\nYour max input altitude is: ", max(alt), " m.")
  
  const auto n = static_cast<Index>(out.size());
  ARTS_ASSERT(n == alt.nelem() and n == lat.nelem() and n == lon.nelem())
  
  if (not regularized) {
    const auto compute = [&](const auto& key, const Data& data) {
      const auto interpolate = [&](const Vector& alts, const Vector& lats, const Vector& lons) -> Vector {
        return detail::vec_interp(data, alts, lats, lons);
      };

      const Vector field_val = interpolate(alt, lat, lon);
      for (Index i=0; i<n; i++) out[i][key] = field_val[i];
    };

    for (auto&& key: keys()) compute(key, operator[](key));
  } else {
    detail::tensor_interpolator(out, *this, alt, lat, lon);
  }
}

std::vector<Point> Field::at(const Vector& alt, const Vector& lat, const Vector& lon) const {
  std::vector<Point> out(alt.nelem());
  at(out, alt, lat, lon);
  return out;
}

//! Regularizes the calculations so that all data is on a single grid
Field& Field::regularize(const Vector& altitudes,
                         const Vector& latitudes,
                         const Vector& longitudes) {
  ARTS_USER_ERROR_IF(regularized, "Cannot re-regularize a regularized grid")

  ArrayOfTensor3 specs_data(
      nspec(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 other_data(
      nother(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));
  ArrayOfTensor3 nlte_data(
      nnlte(),
      Tensor3(altitudes.size(), latitudes.size(), longitudes.size()));

  const auto grids = matpack::repeat(altitudes, latitudes, longitudes);
  const Index nalt = altitudes.nelem(), nlat = latitudes.nelem(),
              nlon = longitudes.nelem();
  const auto pnts =
      matpack::matpack_data<Point, 1>{at(grids[0], grids[1], grids[2])}.reshape(
          nalt, nlat, nlon);

  for (Index i2 = 0; i2 < altitudes.size(); i2++) {
    for (Index i3 = 0; i3 < latitudes.size(); i3++) {
      for (Index i4 = 0; i4 < longitudes.size(); i4++) {
        const auto &pnt = pnts(i2, i3, i4);

        for (Index i0{0}; auto &vals : specs())
          specs_data[i0++](i2, i3, i4) = pnt[vals.first];
        for (Index i0{0}; auto &vals : other())
          other_data[i0++](i2, i3, i4) = pnt[vals.first];
        for (Index i0{0}; auto &vals : nlte())
          nlte_data[i0++](i2, i3, i4) = pnt[vals.first];
      }
    }
  }

  regularized = true;
  grid = {altitudes, latitudes, longitudes};
  top_of_atmosphere = max(altitudes);

  for (Index i0{0}; auto& vals : specs())
    specs()[vals.first] = std::move(specs_data[i0++]);
  for (Index i0{0}; auto& vals : other())
    other()[vals.first] = std::move(other_data[i0++]);
  for (Index i0{0}; auto& vals : nlte())
    nlte()[vals.first] = std::move(nlte_data[i0++]);

  return *this;
}

std::array<Index, 3> Field::regularized_shape() const {
  ARTS_USER_ERROR_IF(not regularized, "Error in atmospheric field:\nCalling a regularized function with an irregular field")
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
  if constexpr (matpack::rank<T>() == 3) {
    return in_data(ind[pos[0]], ind[pos[1]], ind[pos[2]]);
  } else if constexpr (matpack::rank<T>() == 2) {
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

Numeric Point::operator[](Species::Species x) const noexcept {
  for (auto &spec : specs) {
    if (spec.first.Species() == x)
      return spec.second;
  }
  return 0.0;
}

bool Point::is_lte() const noexcept { return nlte.empty(); }

Tensor4 extract_specs_content(const Field &atm,
                              const ArrayOfArrayOfSpeciesTag &specs) {
  Tensor4 out(atm.nspec(), atm.regularized_shape()[0],
              atm.regularized_shape()[1], atm.regularized_shape()[2]);
  std::transform(specs.begin(), specs.end(), out.begin(),
                 [&](auto &spec) { return atm[spec].template get<Tensor3>(); });
  return out;
}

Tensor4 extract_partp_content(const Field &atm, const ArrayOfString &specs) {
  Tensor4 out(atm.npart(), atm.regularized_shape()[0],
              atm.regularized_shape()[1], atm.regularized_shape()[2]);
  std::transform(specs.begin(), specs.end(), out.begin(),
                 [&](auto &spec) { return atm[ParticulatePropertyTag{spec}].template get<Tensor3>(); });
  return out;
}

template <class Key, class T, class Hash, class KeyEqual, class Allocator>
std::vector<Key> get_keys(const std::unordered_map<Key, T, Hash, KeyEqual, Allocator>& map) {
  std::vector<Key> out(map.size());
  std::transform(map.begin(), map.end(), out.begin(), [](auto& v){return v.first;});
  return out;
}

ArrayOfQuantumIdentifier Field::nlte_keys() const {
  return keys<QuantumIdentifier>();
}

void Data::rescale(Numeric x) {
  std::visit(
      [x](auto &v) {
        using T = decltype(v);
        if constexpr (isFunctionalDataType<T>) {
          v = [x, f = v](Numeric alt, Numeric lat, Numeric lon) -> Numeric {
            return x * f(alt, lat, lon);
          };
        } else if constexpr (isGriddedField3<T>) {
          v.data *= x;
        } else {
          v *= x;
        }
      },
      data);
}

void Point::setZero() {
  pressure=0;
  temperature=0;
  wind={0., 0., 0.};
  mag = {0., 0., 0.};
  for (auto& v: specs) v.second = 0;
  for (auto& v: nlte) v.second = 0;
}

Numeric Point::operator[](const KeyVal &k) const {
  return std::visit(
      [this](auto &key) { return this->template operator[](key); }, k);
}

Numeric &Point::operator[](const KeyVal &k) {
  return std::visit(
      [this](auto &key) -> Numeric & { return this->template operator[](key); },
      k);
}
} // namespace Atm

std::ostream &operator<<(std::ostream &os, const ParticulatePropertyTag &ppt) {
  return os << ppt.name;
}
