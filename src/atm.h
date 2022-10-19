#pragma once

#include <limits>
#include <map>
#include <ostream>
#include <type_traits>
#include <utility>
#include <variant>

#include "artstime.h"
#include "debug.h"
#include "enums.h"
#include "gridded_fields.h"
#include "matpack.h"
#include "matpackI.h"
#include "matpackIV.h"
#include "species.h"
#include "species_tags.h"

namespace Atmosphere {
ENUMCLASS(Key,
          char,
          temperature,
          pressure,
          wind_u,
          wind_v,
          wind_w,
          mag_u,
          mag_v,
          mag_w)

class Point {
  Numeric pressure{std::numeric_limits<Numeric>::min()};
  Numeric temperature{std::numeric_limits<Numeric>::min()};
  std::array<Numeric,3> wind{0. ,0. ,0.};
  std::array<Numeric,3> mag{0., 0., 0.};
  std::map<ArrayOfSpeciesTag, Numeric> species_content{};
};

//! All the field data; if these types grow too much we might want to reconsider...
using FieldData = std::variant<GriddedField4, Tensor4, Numeric>;

class Field {
  std::map<Key, FieldData> other;
  std::map<ArrayOfSpeciesTag, FieldData> specs;

  //! The below only exist if regularized is true
  bool regularized{false};
  Vector alt, lat, lon;
  ArrayOfTime time;

  template <typename... Ts>
  void internal_set(const ArrayOfSpeciesTag& lhs,
                    const FieldData& rhs,
                    Ts&&... ts) {
    specs[lhs] = rhs;
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

  template <typename... Ts>
  void internal_set(Key lhs, const FieldData& rhs, Ts&&... ts) {
    other[lhs] = rhs;
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

 public:
  template <typename... Ts>
  Field(Ts&&... ts) {
    static_assert((sizeof...(Ts) % 2) == 0, "Uneven number of inputs");
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

  [[nodiscard]] const FieldData& get(Key x) const { return other.at(x); }

  [[nodiscard]] const FieldData& get(const ArrayOfSpeciesTag& x) const {
    return specs.at(x);
  }

  [[nodiscard]] Shape<4> regularized_shape() const {
    return {time.nelem(), alt.nelem(), lat.nelem(), lon.nelem()};
  }

  FieldData& set(Key x, const FieldData& y) {
    const auto* const t4 = std::get_if<Tensor4>(&y);
    ARTS_USER_ERROR_IF(t4 and t4->shape() not_eq regularized_shape(),
                       "The shape is wrong.  The input field is shape ",
                       t4->shape(),
                       " but we have a regularized shape of ",
                       regularized_shape())
    ARTS_USER_ERROR_IF(
        regularized and not t4,
        "Expects a regularized field (i.e., Tensor4-style gridded)")
    return other[x] = y;
  }

  FieldData& set(const ArrayOfSpeciesTag& x, const FieldData& y) {
    const auto* const t4 = std::get_if<Tensor4>(&y);
    ARTS_USER_ERROR_IF(t4 and t4->shape() not_eq regularized_shape(),
                       "The shape is wrong.  The input field is shape ",
                       t4->shape(),
                       " but we have a regularized shape of ",
                       regularized_shape())
    ARTS_USER_ERROR_IF(
        regularized and not t4,
        "Expects a regularized field (i.e., Tensor4-style gridded)")
    return specs[x] = y;
  }

  //! Regularizes the calculations so that all data is one alt-lat-lon grids, and alt-lat-lon is in the grid map
  void regularize(const ArrayOfTime& times,
                  const Vector& altitudes,
                  const Vector& latitudes,
                  const Vector& longitudes);

  //! Compute the values at a single point
  [[nodiscard]] Point at(Time time_point,
                         Numeric alt_point,
                         Numeric lat_point,
                         Numeric lon_point) const;

  friend std::ostream& operator<<(std::ostream&, const Field&);
};
}  // namespace Atmosphere
