#pragma once

#include <array.h>
#include <mystring.h>

#include <algorithm>

#include "interp.h"
#include "matpack_data.h"

namespace matpack {
template <typename T, typename... Grids>
struct gridded_data {
  static constexpr Size dim = sizeof...(Grids);

  using data_t = matpack::matpack_data<T, dim>;

  using grids_t = std::tuple<Grids...>;

  template <Size Grid>
  using grid_t = std::tuple_element_t<Grid, grids_t>;

  template <Size Grid>
  using grid_value_t = typename grid_t<Grid>::value_type;

  String data_name{};
  data_t data;
  std::array<String, dim> grid_names{};
  grids_t grids;

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

  [[nodiscard]] bool ok() const {
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

    for (Size i : gd.shape()) {
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

  template <Size Grid, my_interp::lagrange_type lag_t>
  [[nodiscard]] Array<lag_t> lag(const grid_t<Grid>& other,
                                 Index order,
                                 Numeric extrapol = 0.5) const
    requires(Grid < dim and std::remove_cvref_t<lag_t>::runtime_polyorder())
  {
    return my_interp::lagrange_interpolation_list<lag_t>(
        other, grid<Grid>(), order, extrapol);
  }

  template <Size Grid, my_interp::lagrange_type lag_t>
  [[nodiscard]] Array<lag_t> lag(const grid_t<Grid>& other,
                                 Numeric extrapol = 0.5) const
    requires(Grid < dim and not std::remove_cvref_t<lag_t>::runtime_polyorder())
  {
    return my_interp::lagrange_interpolation_list<lag_t>(
        other, grid<Grid>(), extrapol);
  }

 private:
  template <my_interp::lagrange_type... lag_ts, Size... sz>
  [[nodiscard]] data_t reinterp(const Grids&... other,
                                Index order,
                                Numeric extrapol,
                                std::integer_sequence<Size, sz...>) const {
    if (not ok()) throw std::runtime_error("bad field");
    return my_interp::reinterp(data,
                               lag<sz, lag_ts>(other, order, extrapol)...);
  }

  template <my_interp::lagrange_type... lag_ts, Size... sz>
  [[nodiscard]] data_t reinterp(const Grids&... other,
                                Numeric extrapol,
                                std::integer_sequence<Size, sz...>) const {
    if (not ok()) throw std::runtime_error("bad field");
    return my_interp::reinterp(data, lag<sz, lag_ts>(other, extrapol)...);
  }

 public:
  template <my_interp::lagrange_type... lag_ts>
  [[nodiscard]] data_t reinterp(const Grids&... other,
                                Index order,
                                Numeric extrapol = 0.5) const
    requires(0 < dim and
             (std::remove_cvref_t<lag_ts>::runtime_polyorder() and ...) and
             sizeof...(lag_ts) == dim)
  {
    if (not ok()) throw std::runtime_error("bad field");
    return reinterp<lag_ts...>(
        other..., order, extrapol, std::make_integer_sequence<Size, dim>{});
  }

  template <my_interp::lagrange_type... lag_ts>
  [[nodiscard]] data_t reinterp(const Grids&... other,
                                Numeric extrapol = 0.5) const
    requires(0 < dim and
             ((not std::remove_cvref_t<lag_ts>::runtime_polyorder()) and
              ...) and
             sizeof...(lag_ts) == dim)
  {
    if (not ok()) throw std::runtime_error("bad field");
    return reinterp<lag_ts...>(
        other..., extrapol, std::make_integer_sequence<Size, dim>{});
  }

  template <Size Grid, my_interp::lagrange_type lag_t>
  [[nodiscard]] lag_t lag(const grid_value_t<Grid>& other, Index order) const
    requires(Grid < dim and std::remove_cvref_t<lag_t>::runtime_polyorder())
  {
    return lag_t(
        std::ranges::lower_bound(grid<Grid>(), other) - grid<Grid>().begin(),
        other,
        grid<Grid>(),
        order);
  }

  template <Size Grid, my_interp::lagrange_type lag_t>
  [[nodiscard]] lag_t lag(const grid_value_t<Grid>& other) const
    requires(Grid < dim and not std::remove_cvref_t<lag_t>::runtime_polyorder())
  {
    return lag_t(
        std::ranges::lower_bound(grid<Grid>(), other) - grid<Grid>().begin(),
        other,
        grid<Grid>());
  }

 private:
  template <my_interp::lagrange_type... lag_ts, Size... sz>
  [[nodiscard]] T interp(const typename Grids::value_type&... other,
                         Index order,
                         std::integer_sequence<Size, sz...>) const {
    if (not ok()) throw std::runtime_error("bad field");
    return my_interp::interp(data, lag<sz, lag_ts>(other, order)...);
  }

  template <my_interp::lagrange_type... lag_ts, Size... sz>
  [[nodiscard]] T interp(const typename Grids::value_type&... other,
                         std::integer_sequence<Size, sz...>) const {
    if (not ok()) throw std::runtime_error("bad field");
    return my_interp::interp(data, lag<sz, lag_ts>(other)...);
  }

 public:
  template <my_interp::lagrange_type... lag_ts>
  [[nodiscard]] T interp(const typename Grids::value_type&... other,
                         Index order) const
    requires(0 < dim and
             (std::remove_cvref_t<lag_ts>::runtime_polyorder() and ...) and
             sizeof...(lag_ts) == dim)
  {
    if (not ok()) throw std::runtime_error("bad field");
    return interp<lag_ts...>(
        other..., order, std::make_integer_sequence<Size, dim>{});
  }

  template <my_interp::lagrange_type... lag_ts>
  [[nodiscard]] T interp(const typename Grids::value_type&... other) const
    requires(0 < dim and
             ((not std::remove_cvref_t<lag_ts>::runtime_polyorder()) and
              ...) and
             sizeof...(lag_ts) == dim)
  {
    if (not ok()) throw std::runtime_error("bad field");
    return interp<lag_ts...>(other..., std::make_integer_sequence<Size, dim>{});
  }
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

template <typename T, typename... Grids>
struct std::formatter<matpack::gridded_data<T, Grids...>> {
  static constexpr Size n = sizeof...(Grids);

  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  using fmt_t = matpack::gridded_data<T, Grids...>;

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
      std::ranges::copy(R"(")"sv, ctx.out());
      std::ranges::copy(v.grid_names[N], ctx.out());
      std::ranges::copy(R"(")"sv, ctx.out());
      std::ranges::copy(": "sv, ctx.out());
    }

    std::formatter<grid_t<N>> gridfmt{};
    tags.compat(gridfmt);
    gridfmt.format(std::get<N>(v.grids), ctx);
    std::ranges::copy(",\n"sv, ctx.out());
  }

 public:
  template <class FmtContext>
  FmtContext::iterator format(const fmt_t& v, FmtContext& ctx) const {
    if (tags.bracket) {
      std::ranges::copy("{\n"sv, ctx.out());
    }

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
      std::ranges::copy(R"(")"sv, ctx.out());
      std::ranges::copy(v.data_name, ctx.out());
      std::ranges::copy(R"(")"sv, ctx.out());
      std::ranges::copy(": "sv, ctx.out());
    }

    std::formatter<matpack::matpack_data<T, n>> datafmt{};
    datafmt.inner_fmt().tags = tags;
    datafmt.format(v.data, ctx);

    if (tags.bracket) {
      std::ranges::copy("\n}"sv, ctx.out());
    }

    return ctx.out();
  }
};
