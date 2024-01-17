#pragma once

#include <array.h>
#include <mystring.h>

#include <algorithm>

#include "matpack_data.h"

namespace matpack {
template <typename T, typename... Grids>
struct gridded_data {
  static constexpr Size dim = sizeof...(Grids);

  using data_t = matpack::matpack_data<T, dim>;

  using grids_t = std::tuple<Grids...>;

  template <Size Grid>
  using grid_t = std::tuple_element_t<Grid, grids_t>;

  data_t data;
  grids_t grids;
  String data_name{};
  std::array<String, dim> grid_names{};

  template <Size Grid>
  [[nodiscard]] grid_t<Grid>& grid()
    requires(Grid < dim)
  {
    return std::get<Grid>(grids);
  }

  template <Size Grid>
  [[nodiscard]] const grid_t<Grid>& grid() const
    requires(Grid < dim)
  {
    return std::get<Grid>(grids);
  }

  template <Size Grid>
  [[nodiscard]] String& gridname()
    requires(Grid < dim)
  {
    return std::get<Grid>(grid_names);
  }

  template <Size Grid>
  [[nodiscard]] const String& gridname() const
    requires(Grid < dim)
  {
    return std::get<Grid>(grid_names);
  }

 private:
  template <Size... Ints>
  [[nodiscard]] std::array<Size, dim> shape(std::index_sequence<Ints...>) const
    requires(sizeof...(Ints) == dim)
  {
    return {static_cast<Size>(grid<Ints>().size())...};
  }

 public:
  [[nodiscard]] std::array<Size, dim> shape() const {
    return shape(std::make_index_sequence<dim>());
  }

  [[nodiscard]] bool check() const {
    const auto a = shape();
    const auto b = data.shape();
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  template <access_operator... Access>
  [[nodiscard]] decltype(auto) operator()(Access&&... access) const {
    return data(std::forward<Access>(access)...);
  }

  template <access_operator... Access>
  [[nodiscard]] decltype(auto) operator()(Access&&... access) {
    return data(std::forward<Access>(access)...);
  }

  template <typename... sz>
  void resize(sz&&... access) {
    data.resize(std::forward<sz>(access)...);
  }

  friend std::ostream& operator<<(std::ostream& os, const gridded_data& gd) {
    if (gd.data_name.size()) {
      os << gd.data_name << " (";
    } else {
      os << "gridded_data (";
    }

    for (Size i: gd.shape()) {
      os << i << ", ";
    }
    os << ") {\n";

    if constexpr (constexpr Size N = 0; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 1; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 2; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 3; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 4; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 5; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 6; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 7; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 8; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    if constexpr (constexpr Size N = 9; dim > N) {
      os << "  " << gd.grid_names[N] << ": " << gd.grid<N>() << '\n';
    }

    static_assert(dim <= 10, "Too many dimensions");

    os << "  data:\n" << gd.data << '\n';

    return os << '}';
  }

  auto operator<=>(const gridded_data&) const = default;
};
}  // namespace matpack

using GriddedField1 = matpack::gridded_data<Numeric, Vector>;
using GriddedField2 = matpack::gridded_data<Numeric, Vector, Vector>;
using GriddedField3 = matpack::gridded_data<Numeric, Vector, Vector, Vector>;
using GriddedField4 =
    matpack::gridded_data<Numeric, Vector, Vector, Vector, Vector>;
using GriddedField5 =
    matpack::gridded_data<Numeric, Vector, Vector, Vector, Vector, Vector>;
using GriddedField6 = matpack::
    gridded_data<Numeric, Vector, Vector, Vector, Vector, Vector, Vector>;

using ComplexGriddedField2 = matpack::gridded_data<Complex, Vector, Vector>;

using NamedGriddedField2 =
    matpack::gridded_data<Numeric, Array<String>, Vector, Vector>;
using NamedGriddedField3 =
    matpack::gridded_data<Numeric, Array<String>, Vector, Vector, Vector>;

using GriddedField1Named =
    matpack::gridded_data<Numeric, Vector, Array<String>>;

using ArrayOfGriddedField1 = Array<GriddedField1>;
using ArrayOfGriddedField2 = Array<GriddedField2>;
using ArrayOfGriddedField3 = Array<GriddedField3>;
using ArrayOfGriddedField4 = Array<GriddedField4>;
using ArrayOfGriddedField5 = Array<GriddedField5>;

using ArrayOfNamedGriddedField2 = Array<NamedGriddedField2>;
using ArrayOfNamedGriddedField3 = Array<NamedGriddedField3>;

using ArrayOfGriddedField1Named = Array<GriddedField1Named>;

using ArrayOfArrayOfGriddedField1 = Array<ArrayOfGriddedField1>;
using ArrayOfArrayOfGriddedField2 = Array<ArrayOfGriddedField2>;
using ArrayOfArrayOfGriddedField3 = Array<ArrayOfGriddedField3>;
using ArrayOfArrayOfGriddedField4 = Array<ArrayOfGriddedField4>;

namespace matpack {
template <typename T, typename... Grids>
std::ostream& operator<<(std::ostream& os,
                         const Array<gridded_data<T, Grids...>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

template <typename T, typename... Grids>
std::ostream& operator<<(std::ostream& os,
                         const Array<Array<gridded_data<T, Grids...>>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

template <typename T, typename... Grids>
std::ostream& operator<<(
    std::ostream& os, const Array<Array<Array<gridded_data<T, Grids...>>>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
}  // namespace matpack
