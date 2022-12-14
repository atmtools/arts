#include "debug.h"
#include "matpack.h"

#include <mdspan/include/experimental/mdspan>

#include <array>
#include <concepts>
#include <functional>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace matpack::md {
using Joker = std::experimental::full_extent_t;
inline constexpr Joker joker = std::experimental::full_extent;

template <typename T>
concept access_operator = std::integral<T> or std::is_same_v<std::remove_cvref_t<Joker>, T>;

namespace detail {
using size_t = Index;

template <std::size_t Rank>
using rank = std::experimental::dextents<size_t, Rank>;

using strided = std::experimental::layout_stride;

//! The type in mdspan for data that is continuous (two items are always
//! naturally close in memory).
template <typename T, detail::size_t N>
using exhaustive_mdspan = std::experimental::mdspan<T, detail::rank<N>>;

//! The type in mdspan for data that is strided (there could be more than the
//! natural step-length between two items).
template <typename T, detail::size_t N>
using strided_mdspan =
    std::experimental::mdspan<T, detail::rank<N>, detail::strided>;

//! Holds true for mdspan(s) that are always continuous in memory, i.e., exhaustive_mdspan and not strided_mdspan
template <typename T>
concept is_always_exhaustive_v = T::is_always_exhaustive();

//! Test the ideas above
static_assert(is_always_exhaustive_v<exhaustive_mdspan<int, 4>>);
static_assert(not is_always_exhaustive_v<strided_mdspan<int, 4>>);
} // namespace detail

// Necessary for different operators to predefine these
template <typename T, detail::size_t N, bool constant> class simple_view;
template <typename T, detail::size_t N, bool constant> class strided_view;

template <access_operator Access, typename... arguments>
consteval detail::size_t get_access_size() noexcept {
  constexpr bool index = std::integral<Access>;
  if constexpr (sizeof...(arguments)) {
    return index + get_access_size<arguments...>();
  } else {
    return index;
  }
}

template <access_operator T, access_operator ... Ts>
bool check_index_sizes(auto&& arr, detail::size_t r, T first, Ts... rest) {
  constexpr bool test = std::integral<T>;
  if constexpr (test) {
    if (first >= arr.extent(r)) return false;
  }

  if constexpr (sizeof...(Ts) == 0) {
    return true;
  } else {
    return check_index_sizes(arr, r+1, rest...);
  }
}

template <std::size_t N>
consteval std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  out.fill(joker);
  return out;
}

template <detail::size_t M, detail::size_t N>
[[nodiscard]] constexpr auto tup(detail::size_t i) {
  if constexpr (M > 0 and N - M - 1 > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{i}, jokers<N - M - 1>());
  } else if constexpr (M > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{i});
  } else if constexpr (N - M - 1 > 0) {
    return std::tuple_cat(std::tuple{i}, jokers<N - M - 1>());
  } else {
    return std::tuple{i};
  }
}

template <detail::size_t M, detail::size_t N, class mdspan_type>
[[nodiscard]] auto sub(mdspan_type &v, detail::size_t i) {
  using namespace std::experimental;
  return std::apply(
      [&v](auto &&...slices) {
        if constexpr (requires{v(slices...);}) {
          return v(slices...);
        } else {
          return submdspan(v, slices...);
        }
      },
      tup<M, N>(i));
}

#define ACCESS_OPS                                                             \
  template <access_operator... access,                                         \
            detail::size_t M = get_access_size<access...>()>                   \
  auto operator()(access... ind)                                               \
      ->std::conditional_t<                                                    \
          M == N, T &,                                                         \
          std::conditional_t<                                                  \
              detail::is_always_exhaustive_v<                                  \
                  decltype(std::experimental::submdspan(view, ind...))>,       \
              simple_view<T, N - M, false>, strided_view<T, N - M, false>>>    \
    requires(sizeof...(access) == N and not constant)                          \
  {                                                                            \
    ARTS_ASSERT(check_index_sizes(view, 0, ind...), "Out-of-bounds")           \
    if constexpr (M == N)                                                      \
      return view(ind...);                                                     \
    else                                                                       \
      return std::experimental::submdspan(view, ind...);                       \
  }                                                                            \
  template <access_operator... access,                                         \
            detail::size_t M = get_access_size<access...>()>                   \
  auto operator()(access... ind) const->std::conditional_t<                    \
      M == N, const T &,                                                       \
      std::conditional_t<                                                      \
          detail::is_always_exhaustive_v<                                      \
              decltype(std::experimental::submdspan(view, ind...))>,           \
          simple_view<T, N - M, true>, strided_view<T, N - M, true>>>          \
    requires(sizeof...(access) == N)                                           \
  {                                                                            \
    ARTS_ASSERT(check_index_sizes(view, 0, ind...), "Out-of-bounds")           \
    if constexpr (M == N)                                                      \
      return view(ind...);                                                     \
    else                                                                       \
      return std::experimental::submdspan(view, ind...);                       \
  }

#define ITERATION_OPS                                                          \
  auto operator[](detail::size_t i)                                            \
      ->std::conditional_t<                                                    \
          N == 1, T &,                                                         \
          std::conditional_t<                                                  \
              detail::is_always_exhaustive_v<decltype(std::apply(              \
                  [&v = view](auto &&...ext) {                                 \
                    return std::experimental::submdspan(v, 0, ext...);         \
                  },                                                           \
                  jokers<N - 1>()))>,                                          \
              simple_view<T, N - 1, false>, strided_view<T, N - 1, false>>>    \
    requires(not constant)                                                     \
  {                                                                            \
    ARTS_ASSERT(i < view.extent(0), "Out of bounds")                           \
    if constexpr (N == 1)                                                      \
      return view[i];                                                          \
    else {                                                                     \
      return sub<0, N>(view, i);                                               \
    }                                                                          \
  }                                                                            \
  auto operator[](detail::size_t i) const->std::conditional_t<                 \
      N == 1, const T &,                                                       \
      std::conditional_t<detail::is_always_exhaustive_v<decltype(std::apply(   \
                             [&v = view](auto &&...ext) {                      \
                               return std::experimental::submdspan(v, 0,       \
                                                                   ext...);    \
                             },                                                \
                             jokers<N - 1>()))>,                               \
                         simple_view<T, N - 1, true>,                          \
                         strided_view<T, N - 1, true>>> {                      \
    ARTS_ASSERT(i < view.extent(0), "Out of bounds")                           \
    if constexpr (N == 1)                                                      \
      return view[i];                                                          \
    else {                                                                     \
      return sub<0, N>(view, i);                                               \
    }                                                                          \
  }

#define ARTS_SIZES                                                             \
  auto size() const noexcept {return view.size();}                             \
  auto extent(detail::size_t i) const noexcept {return view.extent(i);}        \
  auto nelem() const noexcept requires(N==1) {return extent(0);}               \
  auto ncols() const noexcept requires(N>=1) {return extent(0);}               \
  auto nrows() const noexcept requires(N>=2) {return extent(1);}               \
  auto npages() const noexcept requires(N>=3) {return extent(2);}              \
  auto nbooks() const noexcept requires(N>=4) {return extent(3);}              \
  auto nshelves() const noexcept requires(N>=5) {return extent(4);}            \
  auto nvitrines() const noexcept requires(N>=6) {return extent(5);}           \
  auto nlibraries() const noexcept requires(N>=7) {return extent(6);}

template <detail::size_t M, bool constant, class mdspan_type, class submdspan_type> class mditer {
  static constexpr auto N = mdspan_type::rank();
  using T = std::conditional_t<constant, const typename mdspan_type::value_type, typename mdspan_type::value_type>;

  static_assert(N > 0, "Not for 0-rank mdspan");
  static_assert(M >= 0 and M < N, "Out-of-rank");

  detail::size_t pos{0};
  std::conditional_t<constant, const mdspan_type *, mdspan_type *> orig;
public:
  using difference_type = detail::size_t;
  using value_type = std::conditional_t<N == 1, T, submdspan_type>;

  mditer() = default;
  mditer(mditer&&) noexcept = default;
  mditer(const mditer&) = default;
  mditer& operator=(mditer&&) noexcept = default;
  mditer& operator=(const mditer&) = default;

  mditer(mdspan_type &x) requires(not constant) : orig(&x) {}
  mditer(const mdspan_type &x) requires(constant) : orig(&x) {}

  constexpr mditer& operator++() {pos++; return *this;}
  constexpr mditer& operator--() {pos--; return *this;}
  constexpr mditer operator++(int) const {mditer out(*this); ++out.pos; return out;}
  constexpr mditer operator--(int) const {mditer out(*this); --out.pos; return out;}

  constexpr mditer& operator+=(detail::size_t i) {pos+=i; return *this;}
  constexpr mditer& operator-=(detail::size_t i) {pos-=i; return *this;}
  constexpr mditer operator+(detail::size_t i) const {mditer out(*this); out.pos+=i; return out;}
  constexpr mditer operator-(detail::size_t i) const {mditer out(*this); out.pos-=i; return out;}
  constexpr friend mditer operator+(detail::size_t i, const mditer& m) {return m + i;}
  constexpr friend mditer operator-(detail::size_t i, const mditer& m) {return m - i;}

  constexpr difference_type operator-(const mditer& other) const {return pos-other.pos;}
  constexpr difference_type operator+(const mditer& other) const {return pos+other.pos;}

  constexpr auto operator<=>(const mditer&) const = default;

  auto operator*() const -> std::conditional_t<N == 1, std::add_lvalue_reference_t<T>, value_type> {
    if constexpr (N == 1) {
      //return orig->operator[](pos);
      static double x = 2.0;
      return x;
    } else {
      return sub<M, N>(*orig, pos);
    }
  }

  auto operator[](detail::size_t i) const -> std::conditional_t<N == 1, std::add_lvalue_reference_t<T>, value_type> {
    if constexpr (N == 1) {
      return orig->operator[](pos+i);
    } else {
      return sub<M, N>(*orig, pos+i);
    }
  }
};

namespace static_tests {
using T0 = mditer<0, false, detail::exhaustive_mdspan<double, 10>, detail::exhaustive_mdspan<double, 9>>;
using T1 = mditer<0, false, detail::exhaustive_mdspan<double, 2>, detail::exhaustive_mdspan<double, 1>>;
using T2 = mditer<0, false, detail::exhaustive_mdspan<double, 1>, detail::exhaustive_mdspan<double, 0>>;
using T3 = mditer<0, true, detail::exhaustive_mdspan<double, 10>, detail::exhaustive_mdspan<double, 9>>;
using T4 = mditer<0, true, detail::exhaustive_mdspan<double, 2>, detail::exhaustive_mdspan<double, 1>>;
using T5 = mditer<0, true, detail::exhaustive_mdspan<double, 1>, detail::exhaustive_mdspan<double, 0>>;

static_assert(std::random_access_iterator<T0>);
static_assert(std::random_access_iterator<T1>);
static_assert(std::random_access_iterator<T2>);
static_assert(std::random_access_iterator<T3>);
static_assert(std::random_access_iterator<T4>);
static_assert(std::random_access_iterator<T5>);
} // namespace static_tests

template <typename T, detail::size_t N> class simple_data {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  static constexpr bool constant = false;

  std::vector<T> data;
  detail::exhaustive_mdspan<T, N> view;

  template <typename... arguments>
  static constexpr T defdata(arguments&&... args) noexcept
    requires(sizeof...(arguments) == N + 1)
  {
    return static_cast<T>(std::get<sizeof...(arguments) - 1>(
        std::make_tuple(std::forward<arguments>(args)...)));
  }

  template <typename... arguments>
  static constexpr void
  sz_data_impl(std::array<detail::size_t, N> &out,
                         std::integral auto &&x,
                         arguments&&... args [[maybe_unused]])
    requires(sizeof...(arguments) <= N)
  {
    std::get<N - sizeof...(arguments)>(out) =
        static_cast<detail::size_t>(std::forward<decltype(x)>(x));
    if constexpr (sizeof...(arguments) > 1)
      sz_data_impl(out, std::forward<arguments>(args)...);
  }

  template <typename... arguments>
  static constexpr std::array<detail::size_t, N>
  sz_data(arguments&&... args) noexcept
    requires(sizeof...(arguments) == N or sizeof...(arguments) == N + 1)
  {
    if constexpr (sizeof...(arguments) == N) {
      return {
          static_cast<detail::size_t>(std::forward<arguments>(args))...};
    } else {
      std::array<detail::size_t, N> out;
      sz_data_impl(out, std::forward<arguments>(args)...);
      return out;
    }
  }

public:
  // Standard constructor
  simple_data(const std::array<detail::size_t, N> &sz={}, const T &x={})
      : data(std::accumulate(sz.begin(), sz.end(), detail::size_t{1},
                             std::multiplies<>()),
             x),
      view(data.data(), sz) {}

  // Create the type with no known starting value
  template <std::integral... inds>
  simple_data(inds&&... ind)
    requires(sizeof...(inds) == N)
      : simple_data(sz_data(std::forward<inds>(ind)...), 0.0) {
  }

  // Create the type with no known starting value
  template <typename... arguments>
  simple_data(arguments... args)
    requires(sizeof...(args) == N + 1)
      : simple_data(sz_data(std::forward<arguments>(args)...),
                    defdata(std::forward<arguments>(args)...)) {
  }

  using iterator = mditer<0, false, simple_data, simple_view<T, N-1, false>>;
  using const_iterator = mditer<0, true, const simple_data, const simple_view<T, N-1, true>>;
  using value_type = T;
  static constexpr auto rank() {return detail::exhaustive_mdspan<T, N>::rank();}
  auto begin() {if constexpr(N==1) return data.begin(); else return iterator{*this};}
  auto end() {if constexpr(N==1) return data.end(); else return iterator{*this} + extent(0);}
  auto begin() const {if constexpr(N==1) return data.begin(); else return const_iterator{*this};}
  auto end() const {if constexpr(N==1) return data.end(); else return const_iterator{*this} + extent(0);}
  auto cbegin() const {if constexpr(N==1) return data.cbegin(); else return const_iterator{*this};}
  auto cend() const {if constexpr(N==1) return data.cend(); else return const_iterator{*this} + extent(0);}

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES

  friend std::ostream &operator<<(std::ostream &os, const simple_data &sd) {
    return os << simple_view<T, N, true>{sd.view};
  }
};

template <typename T, detail::size_t N, bool constant> class simple_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::exhaustive_mdspan<T, N> view;

  friend class strided_view<T, N, true>;
  friend class strided_view<T, N, false>;

public:
  constexpr simple_view() = default;
  constexpr simple_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr simple_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {} 
  constexpr simple_view(const strided_view<T, N, true>& v) requires(constant) : view(v.view) {}
  constexpr simple_view(const strided_view<T, N, false>& v) : view(v.view) {}

  using iterator = mditer<0, false, simple_view, simple_view<T, N-1, false>>;
  using const_iterator = mditer<0, true, const simple_view, const simple_view<T, N-1, true>>;
  using value_type = T;
  static constexpr auto rank() {return detail::exhaustive_mdspan<T, N>::rank();}
  auto begin() requires(not constant) {return iterator{*this};}
  auto end() requires(not constant) {return iterator{*this} + extent(0);}
  auto begin() const {return const_iterator{*this};}
  auto end() const {return const_iterator{*this} + extent(0);}
  auto cbegin() const {return const_iterator{*this};}
  auto cend() const {return const_iterator{*this} + extent(0);}

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES

  friend std::ostream &operator<<(std::ostream &os, const simple_view &sv) {
    return os << strided_view<T, N, true>{sv};
  }
};

template<typename T, detail::size_t N, bool constant>
class strided_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::strided_mdspan<T, N> view;

  friend class simple_view<T, N, true>;
  friend class simple_view<T, N, false>;

public:
  constexpr strided_view() = default;
  constexpr strided_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {} 
  constexpr strided_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr strided_view(const simple_view<T, N, true>& sv) requires(constant) : view(sv.view) {}
  constexpr strided_view(const simple_view<T, N, false>& sv) : view(sv.view) {}

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES

  [[nodiscard]] constexpr bool is_exhaustive() const noexcept {return view.is_exhaustive();}

  //! A common output operator
  friend std::ostream &operator<<(std::ostream &os, const strided_view &sv) {
    constexpr char space = N == 1 ? ' ' : '\n';

    for (detail::size_t i = 0; i < sv.view.extent(0); i++) {
      if (i not_eq 0) os << space;
      os << sv[i];
    }
    return os;
  }
};

#undef ACCESS_OPS
#undef ITERATION_OPS
#undef ARTS_SIZES
}  // namespace matpack::md

using TMPVector = matpack::md::simple_data<Numeric, 1>;
using TMPMatrix = matpack::md::simple_data<Numeric, 2>;
using TMPTensor3 = matpack::md::simple_data<Numeric, 3>;
