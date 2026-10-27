#pragma once

#include <debug.h>
#include <enumsInterpolationExtrapolation.h>
#include <matpack.h>

#include <array>
#include <functional>
#include <utility>
#include <variant>
#include <vector>

#include "functional_numeric_ternary.h"

/** Interpolation of geodetic field data
 *
 * The data is either gridded with the last two grids being latitude
 * and longitude, a constant, or a function of the position.  The
 * atmospheric, surface, and subsurface fields are all of this kind,
 * whatever their leading grids (altitude, frequency, ...) are.
 */
namespace field_interp {
//! Longitude grids are cyclic, all other grids are not
template <typename Grid> using grid_transform =
    std::conditional_t<std::same_as<Grid, LonGrid>, lagrange_interp::loncross, lagrange_interp::grid_identity>;

//! Linear interpolation, or constant for single-point grids
template <typename Grid> using lags_t =
    std::variant<lagrange_interp::lag_t<0, grid_transform<Grid>>, lagrange_interp::lag_t<1, grid_transform<Grid>>>;

template <typename Grid> lags_t<Grid> lag(const Grid& grid, Numeric x) {
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid");
  if (grid.size() == 1) return grid.template lag<0, grid_transform<Grid>>(x);
  return grid.template lag<1, grid_transform<Grid>>(x);
}

//! Interpolation weights for reuse with interp below
template <typename... Lags> matpack::data_t<Numeric, sizeof...(Lags)> interpweights(const Lags&... lags) {
  return std::visit(
      [](const auto&... l) {
        matpack::data_t<Numeric, sizeof...(Lags)> out(l.size()...);
        lagrange_interp::interpweights(out, l...);
        return out;
      },
      lags...);
}

//! Interpolate data with precomputed weights
template <lagrange_interp::value_transformer transform, typename... Lags>
Numeric interp(const auto& data, const auto& iw, const Lags&... lags) requires(sizeof...(Lags) > 0) {
  return std::visit([&](const auto&... l) { return lagrange_interp::interp<transform>(data, iw, l...); }, lags...);
}

//! Interpolate gridded data at a position
template <lagrange_interp::value_transformer transform = lagrange_interp::value_identity, typename... Grids>
Numeric interp(const matpack::gridded_data_t<Numeric, Grids...>& field,
               const std::array<Numeric, sizeof...(Grids)>&      pos) {
  if (not field.ok())
    arts::user_error(
        "The field data shape {:B,} does not match its grid sizes {:B,}", field.data.shape(), field.shape());

  return [&]<Size... I>(std::index_sequence<I...>) {
    return std::visit([&data = field.data](const auto&... l) { return lagrange_interp::interp<transform>(data, l...); },
                      lag(field.template grid<I>(), pos[I])...);
  }(std::make_index_sequence<sizeof...(Grids)>{});
}

//! Flat indices and weights of the gridded data points contributing to a position
template <typename... Grids>
std::vector<std::pair<Index, Numeric>> flat_weight(const matpack::gridded_data_t<Numeric, Grids...>& field,
                                                   const std::array<Numeric, sizeof...(Grids)>&      pos) {
  if (not field.ok())
    arts::user_error(
        "The field data shape {:B,} does not match its grid sizes {:B,}", field.data.shape(), field.shape());

  constexpr Size N     = sizeof...(Grids);
  const auto     shape = field.data.shape();

  return [&]<Size... I>(std::index_sequence<I...>) {
    return std::visit(
        [&shape](const auto&... l) {
          const std::array<Index, N>             n{l.size()...};
          std::vector<std::pair<Index, Numeric>> out;
          out.reserve((l.size() * ...));

          for (std::array<Index, N> s{}; s.front() < n.front();) {
            std::apply(
                [&](auto... i) {
                  Index idx = 0;
                  Size  d   = 0;
                  ((idx = idx * shape[d++] + l.indx[i]), ...);
                  out.emplace_back(idx, (... * l.data[i]));
                },
                s);

            for (Size d = N; d-- > 0;) {
              if (++s[d] < n[d] or d == 0) break;
              s[d] = 0;
            }
          }

          return out;
        },
        lag(field.template grid<I>(), pos[I])...);
  }(std::make_index_sequence<N>{});
}

/** Retrieval hooks of functional data
 *
 * Only registered function structs that expose x() and w(...) have
 * retrievable values.  See xml_io_stream_functional.
 */
inline ConstVectorView functional_x(const NumericTernary& f) { return ternary::get_xc(f); }

inline VectorView functional_x(NumericTernary& f) { return ternary::get_xm(f); }

inline std::vector<std::pair<Index, Numeric>> functional_w(const NumericTernary& f, const std::array<Numeric, 3>& pos) {
  return ternary::get_w(f, pos[0], pos[1], pos[2]);
}

template <typename R, typename... Args>
concept unregistered_functional =
    std::same_as<typename xml_io_stream_functional<std::function<R(Args...)>>::structs_t, void>;

template <typename R, typename... Args> ConstVectorView functional_x(const std::function<R(Args...)>&) {
  static_assert(unregistered_functional<R, Args...>,
                "This function signature has registered structs, add its retrieval hooks in functional_field_interp.h");
  return {};
}

template <typename R, typename... Args> VectorView functional_x(std::function<R(Args...)>&) {
  static_assert(unregistered_functional<R, Args...>,
                "This function signature has registered structs, add its retrieval hooks in functional_field_interp.h");
  return {};
}

template <typename R, typename... Args>
std::vector<std::pair<Index, Numeric>> functional_w(const std::function<R(Args...)>&,
                                                    const std::array<Numeric, sizeof...(Args)>&) {
  static_assert(unregistered_functional<R, Args...>,
                "This function signature has registered structs, add its retrieval hooks in functional_field_interp.h");
  return {};
}

//! Extrapolation below and above the grid of one dimension
struct extrapolation {
  InterpolationExtrapolation low;
  InterpolationExtrapolation upp;
};

//! The extrapolation type of a position along a grid, clamping the position for nearest-neighbor extrapolation
InterpolationExtrapolation limit(const extrapolation& extrap, const auto& grid, Numeric& x) {
  using enum InterpolationExtrapolation;

  if (x < grid.front()) {
    if (extrap.low == Nearest) x = grid.front();
    return extrap.low;
  }

  if (grid.back() < x) {
    if (extrap.upp == Nearest) x = grid.back();
    return extrap.upp;
  }

  return Linear;
}

//! The strictest extrapolation type: None, then Zero, then Nearest, then Linear
constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a, InterpolationExtrapolation b) {
  using enum InterpolationExtrapolation;
  for (auto type : {None, Zero, Nearest}) {
    if (a == type or b == type) return type;
  }
  return Linear;
}

//! The value of the data at a position, respecting the extrapolation outside of gridded data
template <typename G, typename F, Size N = G::dim> Numeric at(const std::variant<G, Numeric, F>&  data,
                                                              const std::array<extrapolation, N>& extrap,
                                                              const std::array<Numeric, N>&       pos,
                                                              bool                                log = false) {
  if (const auto* field = std::get_if<G>(&data)) {
    using enum InterpolationExtrapolation;

    std::array<Numeric, N>     x    = pos;
    InterpolationExtrapolation type = Linear;
    [&]<Size... I>(std::index_sequence<I...>) {
      ((type = combine(type, limit(extrap[I], field->template grid<I>(), x[I]))), ...);
    }(std::make_index_sequence<N>{});

    if (type == None)
      arts::user_error("Limit breached.  Position {:B,} is out-of-bounds when no extrapolation is wanted", pos);

    if (type == Zero) return 0.0;

    return log ? interp<lagrange_interp::log_transform>(*field, x) : interp(*field, x);
  }

  if (const auto* x = std::get_if<Numeric>(&data)) return *x;

  return std::apply(std::get<F>(data), pos);
}

//! Flat indices and weights of the data points contributing to a position
template <typename G, typename F, Size N = G::dim>
std::vector<std::pair<Index, Numeric>> flat_weight(const std::variant<G, Numeric, F>& data,
                                                   const std::array<Numeric, N>&      pos) {
  return std::visit(
      [&pos]<typename T>(const T& x) -> std::vector<std::pair<Index, Numeric>> {
        if constexpr (std::same_as<T, G>)
          return flat_weight(x, pos);
        else if constexpr (std::same_as<T, Numeric>)
          return {{0, 1.0}};
        else
          return functional_w(x.f, pos);
      },
      data);
}

//! The data as a flat vector
template <typename G, typename F> ConstVectorView flat_view(const std::variant<G, Numeric, F>& data) {
  return std::visit(
      []<typename T>(const T& x) -> ConstVectorView {
        if constexpr (std::same_as<T, G>)
          return x.data.view_as(x.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return ConstVectorView{x};
        else
          return functional_x(x.f);
      },
      data);
}

//! The data as a flat vector
template <typename G, typename F> VectorView flat_view(std::variant<G, Numeric, F>& data) {
  return std::visit(
      []<typename T>(T& x) -> VectorView {
        if constexpr (std::same_as<T, G>)
          return x.data.view_as(x.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{x};
        else
          return functional_x(x.f);
      },
      data);
}

//! Dimensions with a single grid point; all dimensions of constant and functional data
template <typename G, typename F, Size N = G::dim>
std::array<bool, N> single_point_dims(const std::variant<G, Numeric, F>& data) {
  std::array<bool, N> out;
  out.fill(true);

  if (const auto* field = std::get_if<G>(&data)) {
    [&]<Size... I>(std::index_sequence<I...>) {
      ((out[I] = field->template grid<I>().size() == 1), ...);
    }(std::make_index_sequence<N>{});
  }

  return out;
}

//! Gridded data must have valid latitudes and longitudes
template <typename G, typename F> bool ok(const std::variant<G, Numeric, F>& data) {
  const auto* field = std::get_if<G>(&data);
  if (field == nullptr) return true;

  const auto grid_ok = []<typename Grid>(const Grid& grid) {
    if constexpr (std::same_as<Grid, LonGrid>)
      return lagrange_interp::loncross::cycle(grid.front()) == grid.front() and
             lagrange_interp::loncross::cycle(grid.back()) == grid.back();
    else if constexpr (std::same_as<Grid, LatGrid>)
      return grid.front() >= -90 and grid.back() <= 90;
    else
      return true;
  };

  return field->ok() and [&]<Size... I>(std::index_sequence<I...>) {
    return (grid_ok(field->template grid<I>()) and ...);
  }
  (std::make_index_sequence<G::dim>{});
}
}  // namespace field_interp
