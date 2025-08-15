#pragma once

#include <array.h>
#include <configtypes.h>
#include <mystring.h>

#include <algorithm>

#include "lagrange_interp.h"
#include "matpack_mdspan.h"
#include "matpack_mdspan_helpers_grid_t.h"

namespace matpack {
template <typename T, typename Grid0, typename... Grids>
struct gridded_data_t {
  static constexpr Size dim       = 1 + sizeof...(Grids);
  static constexpr bool same_grid = (std::same_as<Grid0, Grids> and ...);

  using data_t = matpack::data_t<T, dim>;

  using grids_t = std::conditional_t<same_grid,
                                     std::array<Grid0, dim>,
                                     std::tuple<Grid0, Grids...>>;

  template <Size Grid>
  using grid_t = std::tuple_element_t<Grid, grids_t>;

  template <Size Grid>
  using grid_value_t = typename grid_t<Grid>::value_type;

  String data_name{};
  data_t data;
  std::array<String, dim> grid_names{};
  grids_t grids;

 private:
  template <Size... sz>
  void grid_copy(auto& out, std::index_sequence<sz...>) const {
    ((out.template grid<sz>() = grid<sz>()), ...);
  }

 public:
  // Allow copy conversion from another gridded_data_t with the same
  // number of grids, but possibly different types.  E.g., sorted to unsorted or vice versa.
  template <typename OtherT, typename OtherGrid0, typename... OtherGrids>
  operator gridded_data_t<OtherT, OtherGrid0, OtherGrids...>() const
    requires(std::is_convertible_v<T, OtherT> and
             std::is_convertible_v<Grid0, OtherGrid0> and
             (std::is_convertible_v<Grids, OtherGrids> and ...) and
             sizeof...(OtherGrids) == dim - 1)
  {
    gridded_data_t<OtherT, OtherGrid0, OtherGrids...> out;
    out.data_name  = data_name;
    out.data       = data;
    out.grid_names = grid_names;

    grid_copy(out, std::make_index_sequence<dim>{});

    return out;
  }

  template <Size Grid>
  [[nodiscard]] constexpr grid_t<Grid>& grid()
    requires(Grid < dim)
  {
    return std::get<Grid>(grids);
  }

  template <Size Grid>
  [[nodiscard]] constexpr const grid_t<Grid>& grid() const
    requires(Grid < dim)
  {
    return std::get<Grid>(grids);
  }

  template <Size Grid>
  [[nodiscard]] constexpr String& gridname()
    requires(Grid < dim)
  {
    return std::get<Grid>(grid_names);
  }

  template <Size Grid>
  [[nodiscard]] constexpr const String& gridname() const
    requires(Grid < dim)
  {
    return std::get<Grid>(grid_names);
  }

 private:
  template <Size... Ints>
  [[nodiscard]] constexpr std::array<Size, dim> shape(
      std::index_sequence<Ints...>) const
    requires(sizeof...(Ints) == dim)
  {
    return {static_cast<Size>(grid<Ints>().size())...};
  }

 public:
  [[nodiscard]] constexpr std::array<Size, dim> shape() const {
    return shape(std::make_index_sequence<dim>());
  }

  [[nodiscard]] constexpr bool ok() const {
    const auto a = shape();
    const auto b = data.shape();
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  template <typename Self, access_operator... Access>
  [[nodiscard]] constexpr decltype(auto) operator[](this Self&& self,
                                                    Access&&... access) {
    return std::forward<Self>(self).data[std::forward<Access>(access)...];
  }

  template <typename... sz>
  constexpr void resize(sz&&... access) {
    data.resize(std::forward<sz>(access)...);
  }

  constexpr auto operator<=>(const gridded_data_t&) const = default;

  template <Size Grid,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(const grid_t<Grid>& other,
                         Index order,
                         Numeric extrapol) const {
    return lagrange_interp::make_lags<transform>(
        grid<Grid>(), other, order, extrapol, gridname<Grid>().c_str());
  }

  template <Size Grid,
            Index N,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(const grid_t<Grid>& other, Numeric extrapol) const {
    return lagrange_interp::make_lags<N, transform>(
        grid<Grid>(), other, extrapol, gridname<Grid>().c_str());
  }

  template <Size Grid,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(const grid_value_t<Grid>& other, Index order) const {
    if constexpr (std::same_as<grid_t<Grid>, AscendingGrid>) {
      return lagrange_interp::lag_t<-1, transform>(
          grid<Grid>(), other, order, lagrange_interp::ascending_grid_t{});
    } else if constexpr (std::same_as<grid_t<Grid>, DescendingGrid>) {
      return lagrange_interp::lag_t<-1, transform>(
          grid<Grid>(), other, order, lagrange_interp::descending_grid_t{});
    } else {
      if (grid<Grid>().size() <= 1 or grid<Grid>()[0] < grid<Grid>()[1]) {
        return lagrange_interp::lag_t<-1, transform>(
            grid<Grid>(), other, order, lagrange_interp::ascending_grid_t{});
      }

      return lagrange_interp::lag_t<-1, transform>(
          grid<Grid>(), other, order, lagrange_interp::descending_grid_t{});
    }
  }

  template <Size Grid,
            Index N,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(const grid_value_t<Grid>& other) const {
    if constexpr (std::same_as<grid_t<Grid>, AscendingGrid>) {
      return lagrange_interp::lag_t<N, transform>(
          grid<Grid>(), other, lagrange_interp::ascending_grid_t{});
    } else if constexpr (std::same_as<grid_t<Grid>, DescendingGrid>) {
      return lagrange_interp::lag_t<N, transform>(
          grid<Grid>(), other, lagrange_interp::descending_grid_t{});
    } else {
      if (grid<Grid>().size() <= 1 or grid<Grid>()[0] < grid<Grid>()[1]) {
        return lagrange_interp::lag_t<N, transform>(
            grid<Grid>(), other, lagrange_interp::ascending_grid_t{});
      }

      return lagrange_interp::lag_t<N, transform>(
          grid<Grid>(), other, lagrange_interp::descending_grid_t{});
    }
  }
};
}  // namespace matpack

using GriddedField1 = matpack::gridded_data_t<Numeric, Vector>;
using GriddedField2 = matpack::gridded_data_t<Numeric, Vector, Vector>;
using GriddedField3 = matpack::gridded_data_t<Numeric, Vector, Vector, Vector>;
using GriddedField4 =
    matpack::gridded_data_t<Numeric, Vector, Vector, Vector, Vector>;
using GriddedField5 =
    matpack::gridded_data_t<Numeric, Vector, Vector, Vector, Vector, Vector>;
using GriddedField6 = matpack::
    gridded_data_t<Numeric, Vector, Vector, Vector, Vector, Vector, Vector>;

using SortedGriddedField1 = matpack::gridded_data_t<Numeric, AscendingGrid>;
using SortedGriddedField2 =
    matpack::gridded_data_t<Numeric, AscendingGrid, AscendingGrid>;
using SortedGriddedField3 = matpack::
    gridded_data_t<Numeric, AscendingGrid, AscendingGrid, AscendingGrid>;
using SortedGriddedField4 = matpack::gridded_data_t<Numeric,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid>;
using SortedGriddedField5 = matpack::gridded_data_t<Numeric,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid>;
using SortedGriddedField6 = matpack::gridded_data_t<Numeric,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid,
                                                    AscendingGrid>;

using CartesianSubsurfaceGriddedField3 = matpack::
    gridded_data_t<Numeric, DescendingGrid, AscendingGrid, AscendingGrid>;

using ComplexGriddedField2 = matpack::gridded_data_t<Complex, Vector, Vector>;

using NamedGriddedField2 =
    matpack::gridded_data_t<Numeric, Array<String>, Vector, Vector>;
using NamedGriddedField3 =
    matpack::gridded_data_t<Numeric, Array<String>, Vector, Vector, Vector>;

using GriddedField1Named =
    matpack::gridded_data_t<Numeric, Vector, Array<String>>;

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

template <typename T, typename... Grids>
struct std::formatter<matpack::gridded_data_t<T, Grids...>> {
  static constexpr Size n = sizeof...(Grids);

  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  using fmt_t = matpack::gridded_data_t<T, Grids...>;

  template <Size N>
  using grid_t = typename fmt_t::template grid_t<N>;

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

 private:
  template <Size N, class FmtContext>
  void grid(const fmt_t& v, FmtContext& ctx) const {
    if (tags.names) {
      const std::string_view quote = tags.quote();
      tags.format(ctx, quote, v.grid_names[N], quote, ": "sv);
    }

    tags.format(ctx, std::get<N>(v.grids), tags.sep(true));
  }

 public:
  template <class FmtContext>
  FmtContext::iterator format(const fmt_t& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '{');
    tags.add_if_bracket(ctx, '\n');

    if constexpr (constexpr Size N = 0; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 1; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 2; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 3; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 4; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 5; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 6; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 7; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 8; n > N) grid<N>(v, ctx);
    if constexpr (constexpr Size N = 9; n > N) grid<N>(v, ctx);

    if (tags.names) {
      const std::string_view quote = tags.quote();
      tags.format(ctx, quote, v.data_name, quote, ": "sv);
    }

    tags.format(ctx, v.data);

    tags.add_if_bracket(ctx, '\n');
    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
