#pragma once

#include <nonstd.h>
#include <xml.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "matpack_mdspan_common_types.h"

namespace lagrange_interp {
namespace {
/*! Defines a transformer concept
 *
 * A transformer is a type that can be used to transform a numeric value
 * to another numeric value.  It is used to transform values in the Lagrange
 * interpolation mechanism.
 *
 * The core expression of a Lagrange interpolation weight is:
 *
 *   wj(x) = (x - x0) * (x - x1) * ... * (x - xm) / ((xj - x0) * (xj - x1) * ... * (xj - xm)) for all j != m.
 *
 * A transformer redefines this to:
 *
 *   wj(x) = (f(x) - f(x0)) * (f(x) - f(x1)) * ... * (f(x) - f(xm)) / ((f(xj) - f(x0)) * (f(xj) - f(x1)) * ... * (f(xj) - f(xm))) for all j != m,
 *
 * where f is the transformation function defined as:
 *
 *   [static] [constexpr] Numeric operator()(Numeric x) [noexcept].
 *
 * The transformer is also responsible for cyclicity of a grid.  The transformation is considered
 * cyclic if the transformer defines the member methods:
 *
 *  static [constexpr] Numeric cycle(Numeric), and
 *  static [consteval] Numeric cycle(), and
 *  static [consteval] Numeric midpoint().
 *
 * The first method cycles the value x to the range defined by the transformer.
 * The second returns the full cycle of the transformer.
 * The third returns the midpoint of the cycle.
 *
 * Note that cyclicity is only dealt with by the cycle method and not by the operator().
 */
template <typename T>
concept transformer = std::same_as<decltype(T{}(Numeric{})), Numeric>;

template <typename T>
concept cyclic =
    transformer<T> and std::same_as<decltype(T::cycle(Numeric{})), Numeric> and
    std::same_as<decltype(T::cycle()), Numeric> and T::cycle() > 0.0 and
    std::same_as<decltype(T::midpoint()), Numeric>;

template <typename T>
concept lagrange_type = requires(T a) {
  a.indx;
  a.data;
  a.size();
};

template <typename T>
concept constant_lagrange_type = lagrange_type<T> and not T::runtime;

template <typename T>
concept runtime_lagrange_type =
    lagrange_type<T> and not constant_lagrange_type<T>;

template <typename T>
concept lagrange_type_list =
    lagrange_type<typename T::value_type> and requires(T a) { a.size(); };

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
inline constexpr bool is_std_array_v = is_std_array<T>::value;

template <typename T>
concept constant_lagrange_type_list =
    lagrange_type_list<T> and is_std_array_v<T>;

template <typename T>
concept fully_constant_lagrange_type_list =
    constant_lagrange_type<typename T::value_type> and
    constant_lagrange_type_list<T>;

template <typename T>
concept runtime_lagrange_type_list =
    lagrange_type_list<T> and not constant_lagrange_type_list<T>;
}  // namespace

//!  A transformer that does not transform the value
struct identity {
  static constexpr Numeric operator()(Numeric x) noexcept { return x; }
};

/*! A transformer that cycles the value to the range [lower, upper)
 *
 * Note that only values [2 * lower - upper, 2 * upper - lower) are
 * cycled.  This is to ensure that the transformation is quick.  If
 * you need a greater range, please implement another transformer
 * and keep this fast.
 */
template <Numeric lower, Numeric upper>
  requires(lower < upper)
struct cycler {
  static consteval Numeric midpoint() noexcept {
    return std::midpoint(upper, lower);
  }

  static consteval Numeric cycle() noexcept { return upper - lower; }

  static constexpr Numeric cycle(Numeric x) noexcept {
    return x + ((x < lower) - (x >= upper)) * cycle();
  }

  static constexpr Numeric operator()(Numeric x) noexcept { return x; }
};

//! [-180, 180) cycler
using loncross = cycler<-180.0, 180.0>;

/******************************************************************
 * Knowledge about the order of the interpolation is very important
 * for the interpolation to work correctly.  The order is defined
 * by the transformer and the size of the input array.
 ******************************************************************/

struct ascending_grid_t {};

struct descending_grid_t {};

namespace {
using order_t = std::variant<ascending_grid_t, descending_grid_t>;

template <typename T>
concept grid_order =
    std::same_as<T, ascending_grid_t> or std::same_as<T, descending_grid_t>;

template <grid_order T>
consteval bool ascending() {
  return std::same_as<T, ascending_grid_t>;
}

/******************************************************************
 * Update the position in the input array so that the closest value
 * for a given polynomial order is in index 0 or 1
 * Polynomial order  |        x <- value here
 * ------------------|--------------------------
 * 0                 |       0 
 * 1                 |       0 1
 * 2                 |    -1 0 1
 * 3                 |    -1 0 1 2
 * 4                 | -2 -1 0 1 2
 * ------------------|--------------------------
 *
 * The update is done by cycling through the input array.  Thus this
 * works best if the input indx list is close to the real value and
 * if this method is called for consecutive values (e.g., as is done
 * in the make_lags methods), it is best if the next x is close to
 * the previous x.
 ******************************************************************/

template <transformer transform, Size X, grid_order grid>
constexpr void update_pos(std::span<Index, X> indx,
                          std::span<const Numeric> xi,
                          Numeric x) {
  const Size n = xi.size();
  const Size P = X == std::dynamic_extent ? indx.size() : X;

  if (n <= P) return;

  if constexpr (cyclic<transform>) x = transform::cycle(x);

  const Size N  = P - 1;
  const Size Of = N / 2;
  std::span<const Numeric>::iterator xp =
      xi.begin() + indx[Of];
  const std::span<const Numeric>::iterator xf =
      xi.begin() + Of;
  const std::span<const Numeric>::iterator xe =
      xi.end() - P / 2 - 1;

  if constexpr (cyclic<transform>) {
    if constexpr (ascending<grid>()) {
      Numeric xl = xi.front(), xu = xi.back();
      if (x < xl or x >= xu) {
        if (N == 0) {
          xu      = xu - (x < xl) * transform::cycle();
          xl      = xl + (x >= xu) * transform::cycle();
          indx[0] = (nonstd::abs(x - xl) > nonstd::abs(x - xu)) * (n - 1);
        } else {
          for (Size i = 0; i < P; i++) indx[i] = (n - Of + i) % n;
        }

        return;
      }
    } else {
      Numeric xu = xi.front(), xl = xi.back();
      if (x < xl or x >= xu) {
        if (N == 0) {
          xu      = xu - (x < xl) * transform::cycle();
          xl      = xl + (x >= xu) * transform::cycle();
          indx[0] = (nonstd::abs(x - xl) < nonstd::abs(x - xu)) * (n - 1);
        } else {
          for (Size i = 0; i < P; i++) indx[i] = (n - Of + i) % n;
        }

        return;
      }
    }
  }

  if constexpr (ascending<grid>()) {
    while (xp < xe and *(xp + 1) < x) ++xp;
    while (xp > xf and *xp > x) --xp;
  } else {
    while (xp < xe and *(xp + 1) > x) ++xp;
    while (xp > xf and *xp < x) --xp;
  }

  if (N == 0) {
    if constexpr (cyclic<transform>) {
      const auto xn = ((xp + 1) == xi.end()) ? xf : (xp + 1);
      xp            = (nonstd::abs(x - *xn) < nonstd::abs(x - *xp)) ? xn : xp;
    } else {
      const auto xn = xp + 1;
      xp = (xn == xi.end() or nonstd::abs(x - *xn) > nonstd::abs(x - *xp)) ? xp
                                                                           : xn;
    }
  }

  if constexpr (cyclic<transform>) {
    indx[0] = (std::clamp(xp, xf, xe) - xf) % n;
    for (Size i = 1; i < P; i++) indx[i] = (indx[i - 1] + 1) % n;
  } else {
    std::iota(indx.begin(), indx.end(), std::clamp(xp, xf, xe) - xf);
  }
}

/******************************************************************
 * Starter method to find a position in the input array, works best
 * with sorted evenly spaced arrays.  
 ******************************************************************/

template <transformer transform, Size X, grid_order grid>
constexpr void find_pos(std::span<Index, X> indx,
                        std::span<const Numeric> xi,
                        Numeric x) {
  const Size n = xi.size();
  const Size P = X == std::dynamic_extent ? indx.size() : X;

  if (n <= P) {
    if constexpr (cyclic<transform>)
      for (Size i = 0; i < P; ++i) indx[i] = i % n;
    else
      for (Size i = 0; i < P; ++i) indx[i] = i;
    return;
  }

  if constexpr (cyclic<transform>) x = transform::cycle(x);

  constexpr auto fractional_index =
      [](Numeric x, Numeric x0, Numeric x1, Index n) static {
        const Numeric frac = (x - x0) / (x1 - x0);
        const Index p0     = static_cast<Index>(frac * (Numeric)(n));
        return std::clamp<Index>(p0, 0, n);
      };

  const Index p0 = fractional_index(
      x, xi.front(), xi.back(), n - (cyclic<transform> ? 0 : P));

  std::iota(indx.begin(), indx.end(), p0);

  if constexpr (cyclic<transform>)
    for (auto& i : indx) i = (i + n) % n;

  update_pos<transform, X, grid>(indx, xi, x);
}

/******************************************************************
 * The Lagrange interpolation weights
 ******************************************************************/

template <transformer transform, Size M, grid_order grid>
constexpr void set_weights(std::span<Numeric, M> data,
                           std::span<const Index, M> indx,
                           std::span<const Numeric> xi,
                           Numeric x) {
  const Size N = (M == std::dynamic_extent) ? (indx.size() - 1) : (M - 1);

  if constexpr (cyclic<transform>) x = transform::cycle(x);
  x = transform{}(x);

  // Last value must make the sum equal to 1
  // And i != j in the internal loop

  if constexpr (cyclic<transform>) {
    if constexpr (ascending<grid>()) {
      if (x < xi.front()) {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);
          if (xj > transform::midpoint()) xj -= transform::cycle();

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m     = indx[k + (k >= j)];  // i != j
            Numeric xm = transform{}(xi[m]);
            if (xm > transform::midpoint()) xm -= transform::cycle();

            numer *= x - xm;
            denom *= xj - xm;
          }

          data[j] = numer / denom;
        }
      } else if (x > xi.back()) {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);
          if (xj < transform::midpoint()) xj += transform::cycle();

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m     = indx[k + (k >= j)];  // i != j
            Numeric xm = transform{}(xi[m]);
            if (xm < transform::midpoint()) xm += transform::cycle();
            numer *= x - xm;
            denom *= xj - xm;
          }
          data[j] = numer / denom;
        }
      } else {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m      = indx[k + (k >= j)];  // i != j
            Numeric xm  = transform{}(xi[m]);
            numer      *= x - xm;
            denom      *= xj - xm;
          }

          data[j] = numer / denom;
        }
      }
    } else {
      if (x < xi.back()) {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);
          if (xj > transform::midpoint()) xj -= transform::cycle();

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m     = indx[k + (k >= j)];  // i != j
            Numeric xm = transform{}(xi[m]);
            if (xm > transform::midpoint()) xm -= transform::cycle();

            numer *= x - xm;
            denom *= xj - xm;
          }

          data[j] = numer / denom;
        }
      } else if (x > xi.front()) {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);
          if (xj < transform::midpoint()) xj += transform::cycle();

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m     = indx[k + (k >= j)];  // i != j
            Numeric xm = transform{}(xi[m]);
            if (xm < transform::midpoint()) xm += transform::cycle();
            numer *= x - xm;
            denom *= xj - xm;
          }
          data[j] = numer / denom;
        }
      } else {
        for (Size j = 0; j < N; ++j) {
          Numeric xj = transform{}(xi[indx[j]]);

          Numeric numer = 1.0;
          Numeric denom = 1.0;

          for (Size k = 0; k < N; ++k) {
            Size m      = indx[k + (k >= j)];  // i != j
            Numeric xm  = transform{}(xi[m]);
            numer      *= x - xm;
            denom      *= xj - xm;
          }

          data[j] = numer / denom;
        }
      }
    }
  } else {
    for (Size j = 0; j < N; ++j) {
      Numeric xj = transform{}(xi[indx[j]]);

      Numeric numer = 1.0;
      Numeric denom = 1.0;

      for (Size k = 0; k < N; ++k) {
        Size m      = indx[k + (k >= j)];  // i != j
        Numeric xm  = transform{}(xi[m]);
        numer      *= x - xm;
        denom      *= xj - xm;
      }

      data[j] = numer / denom;
    }
  }

  data[N] = 1.0;
  for (Size j = 0; j < N; ++j) data[N] -= data[j];
}
}  // namespace

/******************************************************************
 * The core type for Lagrange interpolation
 *
 * Index -1 means that the polynomial order is not fixed
 * and the size of the data is dynamic.  Otherwise, the size is fixed
 * to N + 1, where N is the polynomial order.
 ******************************************************************/

template <Index N, transformer transform = identity>
struct lag_t {
  static constexpr bool cyclic     = lagrange_interp::cyclic<transform>;
  static constexpr bool runtime    = N == -1;
  static constexpr Index PolyOrder = N;
  using transform_t                = transform;

  static_assert(
      N >= -1,
      "N must be -1 or greater, -1 is only accepted for runtime polynomials");

  using data_t = std::
      conditional_t<runtime, std::vector<Numeric>, std::array<Numeric, N + 1>>;

  using indx_t =
      std::conditional_t<runtime, std::vector<Index>, std::array<Index, N + 1>>;

  data_t data{};
  indx_t indx{};

  constexpr lag_t()                            = default;
  constexpr lag_t(const lag_t&)                = default;
  constexpr lag_t(lag_t&&) noexcept            = default;
  constexpr lag_t& operator=(const lag_t&)     = default;
  constexpr lag_t& operator=(lag_t&&) noexcept = default;

  template <grid_order grid>
  constexpr lag_t(std::span<const Numeric> xi,
                  Numeric x,
                  grid)
    requires(not runtime)
  {
    find_pos<transform, N + 1, grid>(indx, xi, x);
    set_weights<transform, N + 1, grid>(data, indx, xi, x);
  }

  template <grid_order grid>
  constexpr lag_t(std::span<const Numeric> xi,
                  Numeric x,
                  Index M,
                  grid)
    requires(runtime)
      : data(M + 1), indx(M + 1) {
    find_pos<transform, std::dynamic_extent, grid>(indx, xi, x);
    set_weights<transform, std::dynamic_extent, grid>(data, indx, xi, x);
  }

  template <grid_order grid>
  constexpr lag_t(indx_t pos,
                  std::span<const Numeric> xi,
                  Numeric x,
                  grid)
      : indx(std::move(pos)) {
    if constexpr (runtime) data.resize(indx.size());
    set_weights<transform, N + 1, grid>(data, indx, xi, x);
  }

  [[nodiscard]] constexpr Index size() const
    requires(runtime)
  {
    return data.size();
  }

  [[nodiscard]] static constexpr Index size()
    requires(not runtime)
  {
    return N + 1;
  }
};

/******************************************************************
 * Common helper functions for linear Lagrange interpolation
 ******************************************************************/

namespace {
template <transformer transform, Size poly>
lag_t<poly, transform> poly_lag(
    std::span<const Numeric> xi, Numeric x) {
  assert(xi.size() > poly);
  if (xi[0] < xi[1]) return lag_t<poly, transform>(xi, x, ascending_grid_t{});
  return lag_t<poly, transform>(xi, x, descending_grid_t{});
}

template <transformer transform, typename T, Size N, Size... Ms>
T variant_lag_helper(std::span<const Numeric> xi,
                     Numeric x) {
  if constexpr (N == 0) {
    return lag_t<0, transform>(xi, x, ascending_grid_t{});
  } else {
    if (xi.size() > N) return poly_lag<transform, N>(xi, x);
    return variant_lag_helper<transform, T, Ms...>(xi, x);
  }
}

template <transformer transform, Size... poly>
auto variant_lag_helper(std::span<const Numeric> xi,
                        Numeric x,
                        std::index_sequence<poly...>) {
  return variant_lag_helper<transform,
                            std::variant<lag_t<poly, transform>...>,
                            sizeof...(poly) - 1 - poly...>(xi, x);
}
}  // namespace

/*!  Creates a variant object that have polynomial possibilities from 0 to N
 *
 * By default, either N=1 (linear interpolation) is used, which means nearest
 * neighbor or linear interpolation.  You can specify other values for N to
 * get higher order polynomial interpolation.
 */
template <transformer transform, Size N = 1>
auto variant_lag(std::span<const Numeric> xi, Numeric x) {
  assert(xi.size() > 0);
  return variant_lag_helper<transform>(
      xi, x, std::make_index_sequence<N + 1>{});
}

/******************************************************************
 * Check limits
 ******************************************************************/

template <transformer transform, Size Extent>
constexpr order_t check_limit(
    const std::span<const Numeric>& xi,
    const std::span<const Numeric, Extent>& xn,
    Numeric extrapolation_limit,
    const Index polyorder,
    const char* info) try {
  const Index n        = xi.size();
  const bool ascending = n <= 1 or xi[0] < xi[1];
  auto retval =
      ascending ? order_t{ascending_grid_t{}} : order_t{descending_grid_t{}};

  if (n == 0) return retval;

  if (polyorder >= n) {
    throw std::runtime_error(
        "Too few grid points for the given polynomial order");
  }

  if constexpr (not cyclic<transform>) {
    if (polyorder == 0 or extrapolation_limit <= 0.0) return retval;

    const auto [minptr, maxptr] = stdr::minmax_element(xn);
    const Numeric xmin          = *minptr;
    const Numeric xmax          = *maxptr;

    const Numeric xmax_lim =
        ascending ? xi.back() + extrapolation_limit * (xi.back() - xi[n - 2])
                  : xi.front() + extrapolation_limit * (xi.front() - xi[1]);

    const Numeric xmin_lim =
        ascending ? xi.front() - extrapolation_limit * (xi[1] - xi.front())
                  : xi.back() - extrapolation_limit * (xi[n - 2] - xi.back());

    if (xmax_lim < xmax or xmin_lim > xmin) {
      throw std::runtime_error(std::format(
          R"(Extrapolation limit ({0}) yields limits to the extrapolation of the grid that are outside the input grid.

These limits are computed from the input grid: [{9}, {10}, ... {11}, {12}]

The maximum value we can extrapolate to is {1} + {0} * ({1} - {2}) = {3}
The minimum value we can extrapolate to is {4} - {0} * ({5} - {4}) = {6}

The actual maximum value is {7}
The actual minimum value is {8}
)",
          extrapolation_limit,
          ascending ? xi.back() : xi.front(),
          ascending ? xi[n - 2] : xi[1],
          xmax_lim,
          ascending ? xi.front() : xi.back(),
          ascending ? xi[1] : xi[n - 2],
          xmin_lim,
          xmax,
          xmin,
          xi.front(),
          xi[1],
          xi[n - 2],
          xi.back()));
    }
  } else {
    if (transform::cycle(xi.front()) != xi.front() or
        transform::cycle(xi.back()) != xi.back()) {
      throw std::runtime_error(std::format(
          "The grid cycles.  This is not allowed.\n\n"
          "The grid covers [{}, {}] but the limits cycle to {} and {}, respectively.",
          xi.front(),
          xi.back(),
          transform::cycle(xi.front()),
          transform::cycle(xi.back())));
    }
  }

  return retval;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error in check_limit for {}:\n{}", info, e.what()));
}

template <transformer transform,
          matpack::ranked_md<1> Orig,
          matpack::ranked_md<1> New>
constexpr order_t check_limit(const Orig& xi,
                              const New& xn,
                              Numeric extrapolation_limit,
                              const Index polyorder,
                              const char* info) {
  if constexpr (matpack::any_cdata<New>) {
    return check_limit<transform, New::ndata>(
        xi, xn, extrapolation_limit, polyorder, info);
  } else {
    return check_limit<transform, std::dynamic_extent>(
        xi, xn, extrapolation_limit, polyorder, info);
  }
}

template <transformer transform, matpack::ranked_md<1> Orig, Size N>
constexpr order_t check_limit(const Orig& xi,
                              const std::array<Numeric, N>& xn,
                              Numeric extrapolation_limit,
                              const Index polyorder,
                              const char* info) {
  return check_limit<transform, N>(
      xi, xn, extrapolation_limit, polyorder, info);
}

/******************************************************************
 * Create vectors of interpolation coordinates and weights
 ******************************************************************/

//! Fixed version of make_lags
template <Size N,
          transformer transform,
          Size Extent,
          class FlagT = lag_t<N, transform>>
constexpr auto make_lags(std::span<const Numeric> xi,
                         std::span<const Numeric, Extent> xn,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  constexpr bool has_dynamic_extent = Extent == std::dynamic_extent;

  const order_t order =
      check_limit<transform, Extent>(xi, xn, extrapolation_limit, N, info);

  if constexpr (has_dynamic_extent) {
    std::vector<FlagT> lags;
    lags.reserve(xn.size());

    std::visit(
        [&]<typename grid>(const grid& ord) {
          if (not xn.empty()) lags.emplace_back(xi, xn.front(), ord);

          for (Size i = 1; i < xn.size(); ++i) {
            const Numeric x = xn[i];
            auto& f         = lags.emplace_back(lags[i - 1]);

            update_pos<transform, N + 1, grid>(f.indx, xi, x);
            set_weights<transform, N + 1, grid>(f.data, f.indx, xi, x);
          }
        },
        order);

    return lags;
  } else {
    std::array<FlagT, Extent> lags{};

    if constexpr (Extent > 0) {
      std::visit(
          [&]<typename grid>(const grid& ord) {
            if (not xn.empty()) lags[0] = FlagT{xi, xn.front(), ord};

            for (Size i = 1; i < xn.size(); ++i) {
              const Numeric x = xn[i];
              auto& f         = lags[i];
              f.indx          = lags[i - 1].indx;

              update_pos<transform, N + 1, grid>(f.indx, xi, x);
              set_weights<transform, N + 1, grid>(f.data, f.indx, xi, x);
            }
          },
          order);
    }

    return lags;
  }
}

template <Size N,
          transformer transform,
          matpack::ranked_md<1> Orig,
          matpack::ranked_md<1> New,
          class FlagT = lag_t<N, transform>>
constexpr auto make_lags(const Orig& xi,
                         const New& xn,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  if constexpr (matpack::any_cdata<New>) {
    return make_lags<N, transform, New::ndata>(
        xi, xn, extrapolation_limit, info);
  } else {
    return make_lags<N, transform, std::dynamic_extent>(
        xi, xn, extrapolation_limit, info);
  }
}

template <Size N,
          transformer transform,
          matpack::ranked_md<1> Orig,
          Size M,
          class FlagT = lag_t<N, transform>>
constexpr auto make_lags(const Orig& xi,
                         const std::array<Numeric, M>& xn,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  return make_lags<N, transform, M>(xi, xn, extrapolation_limit, info);
}

//! Dynamic version of make_lags
template <transformer transform,
          Size Extent,
          class FlagT = lag_t<-1, transform>>
constexpr auto make_lags(std::span<const Numeric> xi,
                         std::span<const Numeric, Extent> xn,
                         const Index polyorder       = 1,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  constexpr bool dynamic_extent = Extent == std::dynamic_extent;

  const order_t order = check_limit<transform, Extent>(
      xi, xn, extrapolation_limit, polyorder, info);

  if constexpr (dynamic_extent) {
    std::vector<FlagT> lags;
    lags.reserve(xn.size());

    std::visit(
        [&]<typename grid>(const grid& ord) {
          if (not xn.empty()) lags.emplace_back(xi, xn.front(), polyorder, ord);

          for (Size i = 1; i < xn.size(); ++i) {
            const Numeric x = xn[i];
            auto& f         = lags.emplace_back(lags[i - 1]);

            update_pos<transform, std::dynamic_extent, grid>(f.indx, xi, x);
            set_weights<transform, std::dynamic_extent, grid>(
                f.data, f.indx, xi, x);
          }
        },
        order);

    return lags;
  } else {
    std::array<FlagT, Extent> lags{};

    if constexpr (Extent > 0) {
      std::visit(
          [&]<typename grid>(const grid& ord) {
            if (not xn.empty()) lags[0] = FlagT{xi, xn.front(), polyorder, ord};

            for (Size i = 1; i < xn.size(); ++i) {
              const Numeric x = xn[i];
              auto& f         = lags[i];
              f.indx          = lags[i - 1].indx;

              update_pos<transform, std::dynamic_extent, grid>(f.indx, xi, x);
              set_weights<transform, std::dynamic_extent, grid>(
                  f.data, f.indx, xi, x);
            }
          },
          order);
    }

    return lags;
  }
}

template <transformer transform,
          matpack::ranked_md<1> Orig,
          matpack::ranked_md<1> New,
          class FlagT = lag_t<-1, transform>>
constexpr auto make_lags(const Orig& xi,
                         const New& xn,
                         const Index polyorder       = 1,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  if constexpr (matpack::any_cdata<New>) {
    return make_lags<transform, New::ndata>(
        xi, xn, polyorder, extrapolation_limit, info);
  } else {
    return make_lags<transform, std::dynamic_extent>(
        xi, xn, polyorder, extrapolation_limit, info);
  }
}

template <transformer transform,
          matpack::ranked_md<1> Orig,
          Size M,
          class FlagT = lag_t<-1, transform>>
constexpr auto make_lags(const Orig& xi,
                         const std::array<Numeric, M>& xn,
                         const Index polyorder       = 1,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") {
  return make_lags<transform, M>(xi, xn, polyorder, extrapolation_limit, info);
}

/******************************************************************
 * Index manipulation
 ******************************************************************/

namespace {
//! Does not support 0-size.  Return early if size is 0.
template <typename Indx, Size N>
constexpr void inc(std::array<Indx, N>& s, const std::array<Indx, N>& n) {
  s.back()++;
  for (Size i = s.size() - 1; i > 0; --i) {
    if (s[i] == n[i]) {
      s[i] = 0;
      ++s[i - 1];
    } else {
      break;
    }
  }
}

template <lagrange_type f0, typename... Ts>
constexpr Size size(const f0& f, const Ts&...) {
  return f.size();
}

template <lagrange_type_list f0, typename... Ts>
constexpr Size size(const f0& f, const Ts&...) {
  return f.size();
}

template <fully_constant_lagrange_type_list T>
consteval Size inner_size() {
  return T::value_type::size();
}
}  // namespace

/******************************************************************
 * Interpolation
 *
 * This reduces the input N-ranked field to a single value.
 *
 * There are 2 types
 *
 * Direct interpolation.
 * Reuse interpolation weights.
 *
 ******************************************************************/

//! Reuse interpolation weights.
template <lagrange_type... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto interp(const matpack::ranked_md<N> auto& field,
                      const matpack::ranked_md<N> auto& itw,
                      const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  T out{};

  const std::array<Index, N> n{itw.shape()};

  for (std::array<Index, N> s{}; s.front() < n.front(); inc(s, n)) {
    out += std::apply(
        [&](auto&&... i) { return field[lags.indx[i]...] * itw[i...]; }, s);
  }

  return out;
}

//! Direct interpolation.
template <lagrange_type... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto interp(const matpack::ranked_md<N> auto& field,
                      const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  T out{};

  const std::array<Index, N> n{lags.size()...};

  for (std::array<Index, N> s{}; s.front() < n.front(); inc(s, n)) {
    out += std::apply(
        [&](auto&&... i) {
          return (field[lags.indx[i]...] * ... * lags.data[i]);
        },
        s);
  }

  return out;
}

/******************************************************************
 * Re-interpolation
 *
 * This retains the rank of the input field on new coordinates.
 *
 * Tip: use reshape on the rvalue to reduce ranks if any are 1.
 *
 * There are 4 types
 *
 * Direct reinterpolation.
 * Reuse interpolation weights.
 * Reuse output.
 * Reuse output and interpolation weights.
 *
 ******************************************************************/

//! Reuse output and interpolation weights.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void reinterp(matpack::mut_ranked_md<N> auto&& out,
                        const matpack::ranked_md<N> auto& field,
                        const matpack::ranked_md<2 * N> auto& itw,
                        const FlagTs&... lags) {
  if (out.empty()) return;

  const std::array<Index, N> n{out.shape()};

  for (std::array<Index, N> s{}; s.front() < n.front(); inc(s, n)) {
    std::apply(
        [&](auto&&... i) { out[i...] = interp(field, itw[i...], lags[i]...); },
        s);
  }
}

//! Reuse output.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void reinterp(matpack::mut_ranked_md<N> auto&& out,
                        const matpack::ranked_md<N> auto& field,
                        const FlagTs&... lags) {
  if (out.empty()) return;

  const std::array<Index, N> n{out.shape()};
  for (std::array<Index, N> s{}; s.front() < n.front(); inc(s, n)) {
    std::apply([&](auto&&... i) { out[i...] = interp(field, lags[i]...); }, s);
  }
}

//! Reuse interpolation weights.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto reinterp(const matpack::ranked_md<N> auto& field,
                        const matpack::ranked_md<2 * N> auto& itw,
                        const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  if constexpr ((constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<T, FlagTs::size()...> out{};
    reinterp(out, field, itw, lags...);
    return out;
  } else {
    matpack::data_t<T, N> out(lags.size()...);
    reinterp(out, field, itw, lags...);
    return out;
  }
}

//! Direct reinterpolation.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto reinterp(const matpack::ranked_md<N> auto& field,
                        const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  if constexpr ((constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<T, FlagTs::size()...> out{};
    reinterp(out, field, lags...);
    return out;
  } else {
    matpack::data_t<T, N> out(lags.size()...);
    reinterp(out, field, lags...);
    return out;
  }
}

/******************************************************************
 * Flat interpolation
 *
 * This reduces the input N-ranked field to a vector of values.
 *
 * There are 4 types
 *
 * Direct interpolation.
 * Reuse output.
 * Reuse interpolation weights.
 * Reuse output and interpolation weights.
 *
 ******************************************************************/

//! Reuse output and interpolation weights.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void flat_interp(matpack::mut_ranked_md<1> auto&& out,
                           const matpack::ranked_md<N> auto& field,
                           const matpack::ranked_md<N + 1> auto& itw,
                           const FlagTs&... lags) {
  const Size n = out.size();
  for (Size i = 0; i < n; ++i) out[i] = interp(field, itw[i], lags[i]...);
}

//! Reuse output.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void flat_interp(matpack::mut_ranked_md<1> auto&& out,
                           const matpack::ranked_md<N> auto& field,
                           const FlagTs&... lags) {
  const Size n = out.size();
  for (Size i = 0; i < n; ++i) out[i] = interp(field, lags[i]...);
}

//! Reuse interpolation weights.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto flat_interp(const matpack::ranked_md<N> auto& field,
                           const matpack::ranked_md<N + 1> auto& itw,
                           const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  if constexpr ((constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<T, size(FlagTs{}...)> out{};
    flat_interp(out, field, itw, lags...);
    return out;
  } else {
    matpack::data_t<T, 1> out(size(lags...));
    flat_interp(out, field, itw, lags...);
    return out;
  }
}

//! Direct interpolation.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto flat_interp(const matpack::ranked_md<N> auto& field,
                           const FlagTs&... lags) {
  using T = std::remove_cvref_t<decltype(*field.elem_begin())>;

  if constexpr ((constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<T, size(FlagTs{}...)> out{};
    flat_interp(out, field, lags...);
    return out;
  } else {
    matpack::data_t<T, 1> out(size(lags...));
    flat_interp(out, field, lags...);
    return out;
  }
}

/******************************************************************
 * Interpolation weights 
 *
 * There are 6 types:
 *
 * Direct for interp.
 * Reuse memory for interp.
 * Direct for flat_interp.
 * Reuse memory for flat_interp.
 * Direct for reinterp.
 * Reuse memory for reinterp.
 *
 * The rank for interp is the number of lags.
 * The rank for flat_interp is the number of lags plus one.
 * The rank for reinterp is twice the number of lags.
 *
 ******************************************************************/

//! Reuse memory for interp.
template <lagrange_type... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void interpweights(matpack::mut_ranked_md<N> auto&& itw,
                             const FlagTs&... lags) {
  const std::array<Index, N> n{itw.shape()};
  for (std::array<Index, N> s{}; s.front() < n.front(); inc(s, n)) {
    std::apply([&](auto&&... i) { itw[i...] = (... * lags.data[i]); }, s);
  }
}

//! Reuse memory for flat_interp.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void flat_interpweights(matpack::mut_ranked_md<1 + N> auto&& itw,
                                  const FlagTs&... lags) {
  const Size n = itw.extent(0);
  for (Size i = 0; i < n; ++i) interpweights(itw[i], lags[i]...);
}

//! Reuse memory for reinterp.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr void reinterpweights(matpack::mut_ranked_md<2 * N> auto&& itw,
                               const FlagTs&... lags) {
  if (itw.empty()) return;

  const std::array<Size, N> n{lags.size()...};
  for (std::array<Size, N> s{}; s.front() < n.front(); inc(s, n)) {
    std::apply([&](auto&&... i) { interpweights(itw[i...], lags[i]...); }, s);
  }
}

//! Direct for interp.
template <lagrange_type... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto interpweights(const FlagTs&... lags) {
  if constexpr ((constant_lagrange_type<FlagTs> and ...)) {
    matpack::cdata_t<Numeric, FlagTs::size()...> out{};
    interpweights(out, lags...);
    return out;
  } else {
    matpack::data_t<Numeric, N> out(lags.size()...);
    interpweights(out, lags...);
    return out;
  }
}

//! Direct for flat_interp.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto flat_interpweights(const FlagTs&... lags) {
  if constexpr ((fully_constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<Numeric, size(FlagTs{}...), inner_size<FlagTs>()...> out{};
    flat_interpweights(out, lags...);
    return out;
  } else {
    matpack::data_t<Numeric, 1 + N> out(
        size(lags...), (lags.size() ? lags.front().size() : 0)...);
    flat_interpweights(out, lags...);
    return out;
  }
}

//! Direct for reinterp.
template <lagrange_type_list... FlagTs, Size N = sizeof...(FlagTs)>
constexpr auto reinterpweights(const FlagTs&... lags) {
  if constexpr ((fully_constant_lagrange_type_list<FlagTs> and ...)) {
    matpack::cdata_t<Numeric, FlagTs::size()..., inner_size<FlagTs>()...> out{};
    reinterpweights(out, lags...);
    return out;
  } else {
    matpack::data_t<Numeric, 2 * N> out(
        lags.size()..., (lags.size() ? lags.front().size() : 0)...);
    reinterpweights(out, lags...);
    return out;
  }
}
}  // namespace lagrange_interp

template <Index N, lagrange_interp::transformer transform>
struct std::formatter<lagrange_interp::lag_t<N, transform>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lagrange_interp::lag_t<N, transform>& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.data, tags.sep(), v.indx);
  }
};

template <Index N, lagrange_interp::transformer transform>
struct xml_io_stream_name<lagrange_interp::lag_t<N, transform>> {
  static constexpr std::string_view name = "lagrange_interp";
};

template <Index N, lagrange_interp::transformer transform>
struct xml_io_stream<lagrange_interp::lag_t<N, transform>> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<lagrange_interp::lag_t<N, transform>>;

  static void write(std::ostream& os,
                    const lagrange_interp::lag_t<N, transform>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    XMLTag tag(type_name, "name"sv, name, "N"sv, N);
    tag.write_to_stream(os);

    xml_write_to_stream(os, x.data, pbofs, "data");
    xml_write_to_stream(os, x.indx, pbofs, "indx");

    tag.write_to_end_stream(os);
  }

  static void read(std::istream& is,
                   lagrange_interp::lag_t<N, transform>& x,
                   bifstream* pbifs = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);

    tag.check_name(type_name);
    tag.check_attribute("N", N);

    xml_read_from_stream(is, x.data, pbifs);
    xml_read_from_stream(is, x.indx, pbifs);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};
