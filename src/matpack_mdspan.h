#pragma once

#include "debug.h"
#include "matpack.h"

#include "matpack_mdspan_concepts.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace matpack::md {
namespace detail {
template <Index Rank>
using rank = std::experimental::dextents<Index, Rank>;

using strided = std::experimental::layout_stride;

//! The type in mdspan for data that is continuous (two items are always
//! naturally close in memory).
template <typename T, Index N>
using exhaustive_mdspan = std::experimental::mdspan<T, detail::rank<N>>;

//! The type in mdspan for data that is strided (there could be more than the
//! natural step-length between two items).
template <typename T, Index N>
using strided_mdspan =
    std::experimental::mdspan<T, detail::rank<N>, detail::strided>;

//! Holds true for mdspan(s) that are always continuous in memory, i.e., exhaustive_mdspan and not strided_mdspan
template <typename T>
concept is_always_exhaustive_v = T::is_always_exhaustive();

//! Test the ideas above
static_assert(is_always_exhaustive_v<exhaustive_mdspan<int, 4>>);
static_assert(not is_always_exhaustive_v<strided_mdspan<int, 4>>);
} // namespace detail

//! A strided accessor to get the range of a type
struct strided_access {
  Index offset{0};
  Index extent{-1};
  Index stride{1};
  constexpr strided_access(Index i0=0, Index len=-1, Index d=1) noexcept : offset(i0), extent(len), stride(d) {}
  constexpr strided_access(Index i0, Joker, Index d=1) noexcept : offset(i0), stride(d) {}
  constexpr strided_access(Joker, Index d=1) noexcept : stride(d) {}
  constexpr strided_access(Index max_size, const strided_access& r)  ARTS_NOEXCEPT
      : offset(r.offset),
        extent(r.extent),
        stride(r.stride) {
    // Start must be >= 0:
    ARTS_ASSERT(0 <= offset, "offset=", offset);
    // ... and < max_size:
    ARTS_ASSERT(offset < max_size or (max_size == 0 and offset == 0), "offset=", offset, "; max_size=", max_size);

    // Stride must be != 0:
    ARTS_ASSERT(0 != stride, "stride=", stride);

    // Convert negative extent (joker) to explicit extent
    if (extent < 0) {
      if (0 < stride)
        extent = 1 + (max_size - 1 - offset) / stride;
      else
        extent = 1 + (0 - offset) / stride;
    } else {
      // Check that extent is ok:
      ARTS_ASSERT(0 <= (offset + (extent - 1) * stride));
      ARTS_ASSERT((offset + (extent - 1) * stride) < max_size);
    }
  }
  constexpr strided_access(const strided_access& p, const strided_access& n) ARTS_NOEXCEPT
      : offset(p.offset + n.offset * p.stride),
        extent(n.extent),
        stride(p.stride* n.stride) {
    // We have to juggle here a bit with previous, new, and resulting
    // quantities. I.e.;
    // p.mstride: Previous stride
    // n.mstride: New stride (as specified)
    // mstride:   Resulting stride (old*new)

    // Get the previous final element:
    Index const prev_fin = p.offset + (p.extent - 1) * p.stride;

    // Resulting start must be >= previous start:
    ARTS_ASSERT(p.offset <= offset);
    // and <= prev_fin, except for Joker:
    ARTS_ASSERT(offset <= prev_fin || extent == -1);

    // Resulting stride must be != 0:
    ARTS_ASSERT(0 != stride,
                "Problem: stride==",
                stride,
                " for offset==",
                offset,
                " and extent==",
                extent);

    // Convert negative extent (joker) to explicit extent
    if (extent < 0) {
      if (0 < stride)
        extent = 1 + (prev_fin - offset) / stride;
      else
        extent = 1 + (p.offset - offset) / stride;
    } else {
      ARTS_ASSERT(p.offset <= (offset + (extent - 1) * stride));
      ARTS_ASSERT((offset + (extent - 1) * stride) <= prev_fin);
    }
  }

  constexpr operator Joker() const noexcept {return joker;}
  constexpr operator Joker() noexcept {return joker;}
};

template <access_operator Access, typename... arguments>
[[nodiscard]] consteval Index get_access_size() noexcept {
  constexpr bool index = integral<Access>;
  if constexpr (sizeof...(arguments)) {
    return index + get_access_size<arguments...>();
  } else {
    return index;
  }
}

template <access_operator first, access_operator... access, Index N=sizeof...(access)>
consteval bool any_range() {
  if constexpr (std::same_as<strided_access, std::remove_cvref_t<first>>) return true;
  else if constexpr (N not_eq 0) return any_range<access...>();
  else return false;
};

template <access_operator T, access_operator ... Ts>
[[nodiscard]] bool check_index_sizes(auto&& arr, Index r, T first, Ts... rest) {
  constexpr bool test = integral<T>;
  if constexpr (test) {
    if (first >= arr.extent(r)) return false;
  }

  if constexpr (sizeof...(Ts) == 0) {
    return true;
  } else {
    return check_index_sizes(arr, r+1, rest...);
  }
}

template <Index N>
[[nodiscard]] consteval std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  out.fill(joker);  // I do not know if this is necessary, but it is just a pure compile-time cost
  return out;
}

template <Index M, Index N>
[[nodiscard]] constexpr auto tup(Index i) {
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

template <Index M, Index N, class mdspan_type>
[[nodiscard]] auto sub(mdspan_type &v, Index i) {
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
            Index M = get_access_size<access...>()>                            \
  [[nodiscard]] auto operator()(access... ind)                                 \
      ->std::conditional_t<                                                    \
          M == N, T &,                                                         \
          std::conditional_t<                                                  \
              not any_range<access...>() and                                   \
                  detail::is_always_exhaustive_v<                              \
                      decltype(std::experimental::submdspan(view, ind...))>,   \
              simple_view<T, N - M, false>, strided_view<T, N - M, false>>>    \
    requires(sizeof...(access) == N and not constant)                          \
  {                                                                            \
    ARTS_ASSERT(check_index_sizes(view, 0, ind...), "Out-of-bounds")           \
    if constexpr (M == N)                                                      \
      return view(ind...);                                                     \
    else if constexpr (any_range<access...>())                                 \
      return strided_view<T, N - M, false>{                                    \
          std::experimental::submdspan(view, ind...)}                          \
          .range_adaptor(std::forward<access>(ind)...);                        \
    else                                                                       \
      return std::experimental::submdspan(view, ind...);                       \
  }                                                                            \
  template <access_operator... access,                                         \
            Index M = get_access_size<access...>()>                   \
  [[nodiscard]] auto operator()(access... ind) const->std::conditional_t<      \
      M == N, const T &,                                                       \
      std::conditional_t<                                                      \
          not any_range<access...>() and                                       \
              detail::is_always_exhaustive_v<                                  \
                  decltype(std::experimental::submdspan(view, ind...))>,       \
          simple_view<T, N - M, true>, strided_view<T, N - M, true>>>          \
    requires(sizeof...(access) == N)                                           \
  {                                                                            \
    ARTS_ASSERT(check_index_sizes(view, 0, ind...), "Out-of-bounds")           \
    if constexpr (M == N)                                                      \
      return view(ind...);                                                     \
    else if constexpr (any_range<access...>())                                 \
      return strided_view<T, N - M, true>{                                     \
          std::experimental::submdspan(view, ind...)}                          \
          .range_adaptor(std::forward<access>(ind)...);                        \
    else                                                                       \
      return std::experimental::submdspan(view, ind...);                       \
  }

#define ITERATION_OPS                                                          \
  [[nodiscard]] auto operator[](Index i)                              \
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
  [[nodiscard]] auto operator[](Index i) const->std::conditional_t<   \
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
  }                                                                            \
  [[nodiscard]] auto operator[](strided_access range) {                        \
    strided_view<T, N, constant> out{view};                                    \
    out.template strided_dimension<0>(range);                                  \
    return out;                                                                \
  }                                                                            \
  [[nodiscard]] auto operator[](strided_access range) const {                  \
    strided_view<T, N, true> out{view};                                        \
    out.template strided_dimension<0>(range);                                  \
    return out;                                                                \
  }

#define ARTS_SIZES                                                             \
  [[nodiscard]] auto size() const noexcept { return view.size(); }             \
  [[nodiscard]] auto stride(Index i) const noexcept {                 \
    return view.stride(i);                                                     \
  }                                                                            \
  [[nodiscard]] auto extent(Index i) const noexcept {                 \
    return view.extent(i);                                                     \
  }                                                                            \
  [[nodiscard]] auto nelem() const noexcept                                    \
    requires(N == 1)                                                           \
  {                                                                            \
    return extent(0);                                                          \
  }                                                                            \
  [[nodiscard]] auto ncols() const noexcept                                    \
    requires(N >= 1)                                                           \
  {                                                                            \
    return extent(N - 1);                                                      \
  }                                                                            \
  [[nodiscard]] auto nrows() const noexcept                                    \
    requires(N >= 2)                                                           \
  {                                                                            \
    return extent(N - 2);                                                      \
  }                                                                            \
  [[nodiscard]] auto npages() const noexcept                                   \
    requires(N >= 3)                                                           \
  {                                                                            \
    return extent(N - 3);                                                      \
  }                                                                            \
  [[nodiscard]] auto nbooks() const noexcept                                   \
    requires(N >= 4)                                                           \
  {                                                                            \
    return extent(N - 4);                                                      \
  }                                                                            \
  [[nodiscard]] auto nshelves() const noexcept                                 \
    requires(N >= 5)                                                           \
  {                                                                            \
    return extent(N - 5);                                                      \
  }                                                                            \
  [[nodiscard]] auto nvitrines() const noexcept                                \
    requires(N >= 6)                                                           \
  {                                                                            \
    return extent(N - 6);                                                      \
  }                                                                            \
  [[nodiscard]] auto nlibraries() const noexcept                               \
    requires(N >= 7)                                                           \
  {                                                                            \
    return extent(N - 7);                                                      \
  }                                                                            \
  [[nodiscard]] std::array<Index, N> shape() const noexcept {         \
    std::array<Index, N> out;                                         \
    for (Index i = 0; i < N; i++)                                     \
      out[i] = extent(i);                                                      \
    return out;                                                                \
  }                                                                            \
  [[nodiscard]] bool empty() const { return view.empty(); }

#define ITERATORS(IN, OUT)                                                     \
  using iterator =                                                             \
      mditer<0, false, IN, OUT<T, std::max<Index>(N - 1, 1), false>>; \
  using const_iterator =                                                       \
      mditer<0, true, const IN,                                                \
             const OUT<T, std::max<Index>(N - 1, 1), true>>;          \
  using value_type = T;                                                        \
  [[nodiscard]] static constexpr auto rank() { return N; }                     \
  [[nodiscard]] static constexpr auto is_const() { return constant; }          \
  [[nodiscard]] auto begin()                                                   \
    requires(not constant)                                                     \
  {                                                                            \
    return iterator{*this};                                                    \
  }                                                                            \
  [[nodiscard]] auto end()                                                     \
    requires(not constant)                                                     \
  {                                                                            \
    return iterator{*this} + extent(0);                                        \
  }                                                                            \
  [[nodiscard]] auto begin() const { return const_iterator{*this}; }           \
  [[nodiscard]] auto end() const { return const_iterator{*this} + extent(0); } \
  [[nodiscard]] auto cbegin() const { return const_iterator{*this}; }          \
  [[nodiscard]] auto cend() const {                                            \
    return const_iterator{*this} + extent(0);                                  \
  }

template <Index M, bool constant, class mdspan_type,
          class submdspan_type>
class mditer {
  static constexpr auto N = mdspan_type::rank();
  using T = std::conditional_t<constant, const typename mdspan_type::value_type,
                               typename mdspan_type::value_type>;

  static_assert(N > 0, "Not for 0-rank mdspan");
  static_assert(M >= 0 and M < N, "Out-of-rank");

  Index pos{0};
  std::conditional_t<constant, const mdspan_type *, mdspan_type *> orig;

public:
  using difference_type = Index;
  using value_type = std::conditional_t<N == 1, T, submdspan_type>;

  mditer() = default;
  mditer(mditer&&) noexcept = default;
  mditer(const mditer&) = default;
  mditer& operator=(mditer&&) noexcept = default;
  mditer& operator=(const mditer&) = default;

  mditer(mdspan_type &x) requires(not constant) : orig(&x) {}
  mditer(const mdspan_type &x) requires(constant) : orig(&x) {}

  constexpr mditer& operator++() noexcept {pos++; return *this;}
  constexpr mditer& operator--() noexcept {pos--; return *this;}
  [[nodiscard]] constexpr mditer operator++(int) const noexcept {mditer out(*this); ++out.pos; return out;}
  [[nodiscard]] constexpr mditer operator--(int) const noexcept {mditer out(*this); --out.pos; return out;}

  constexpr mditer& operator+=(Index i) noexcept {pos+=i; return *this;}
  constexpr mditer& operator-=(Index i) noexcept {pos-=i; return *this;}
  [[nodiscard]] constexpr mditer operator+(Index i) const noexcept {mditer out(*this); out.pos+=i; return out;}
  [[nodiscard]] constexpr mditer operator-(Index i) const noexcept {mditer out(*this); out.pos-=i; return out;}
  [[nodiscard]] constexpr friend mditer operator+(Index i, const mditer& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend mditer operator-(Index i, const mditer& m) noexcept {return m - i;}

  [[nodiscard]] constexpr difference_type operator-(const mditer& other) const noexcept {return pos-other.pos;}
  [[nodiscard]] constexpr difference_type operator+(const mditer& other) const noexcept {return pos+other.pos;}

  [[nodiscard]] constexpr auto operator<=>(const mditer&) const noexcept = default;

  [[nodiscard]] auto operator*() const
      -> std::conditional_t<N == 1, std::add_lvalue_reference_t<T>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->operator[](pos);
    } else {
      return sub<M, N>(*orig, pos);
    }
  }

  [[nodiscard]] auto operator[](Index i) const
      -> std::conditional_t<N == 1, std::add_lvalue_reference_t<T>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->operator[](pos + i);
    } else {
      return sub<M, N>(*orig, pos + i);
    }
  }
};

template <Index N> std::array<Index, N> ones() {
  std::array<Index, N> x;
  x.fill(1);
  return x;
}

template <math_type T, Index N, bool constant> class simple_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::exhaustive_mdspan<T, N> view;

  template <math_type U, Index M> friend class simple_data;
  template <math_type U, Index M, bool c> friend class simple_view;
  template <math_type U, Index M, bool c> friend class strided_view;

public:
  constexpr simple_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr simple_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr simple_view(T* data, std::array<Index, N> sz) : view(data, sz) {}
  template <integral ... inds, Index M = sizeof...(inds)>
  simple_view(T* data, inds ... ind) requires(M == N) : view(data, ind...) {}
  simple_view(T& x) : view(&x, ones<N>()) {}

  constexpr simple_view() = default;
  constexpr simple_view(const simple_view&) = default;
  constexpr simple_view& operator=(const simple_view&) = default;
  constexpr simple_view(simple_view&&) noexcept = default;
  constexpr simple_view& operator=(simple_view&&) noexcept = default;

  template <matpack_convertible U>
  simple_view& operator=(const U& v) {return *this = strided_view<T, N, true>(v);}

  constexpr operator simple_view<T, N, true>() const requires(not constant) {return view;}
  constexpr operator strided_view<T, N, false>() {return view;}
  constexpr operator strided_view<T, N, true>() const {return view;}
  template <Index M>
  explicit constexpr operator simple_view<T, M, false>() requires(M > N) {
    std::array<Index, M> sz;
    for (Index i=0; i<M; i++) sz[i] = i < N ? extent(i) : 1;
    return simple_view<T, M, constant> {view.data_handle(), sz};
  }
  template <Index M>
  explicit constexpr operator strided_view<T, M, false>() requires(M > N) {
    return simple_view<T, M, constant>{*this};
  }
  template <Index M>
  explicit constexpr operator simple_view<T, M, true>() const requires(M > N) {
    std::array<Index, M> sz;
    for (Index i=0; i<M; i++) sz[i] = i < N ? extent(i) : 1;
    return simple_view<T, M, constant> {view.data_handle(), sz};
  }
  template <Index M>
  explicit constexpr operator strided_view<T, M, true>() const requires(M > N) {
    return simple_view<T, M, constant>{*this};
  }

  constexpr simple_view &operator=(const any_same_md<N> auto &sv)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == sv.shape(), "Mismatch shape")
    if (view.data_handle() not_eq sv.view.data_handle()) {
      std::copy(sv.begin(), sv.end(), begin());}
    return *this;
  }

  constexpr simple_view& operator=(const T& x) requires (not constant) {
    std::fill(data_handle(), data_handle()+size(), x);
    return *this;
  }

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES
  ITERATORS(simple_view, simple_view)

  //! ARTS Helper functions
  [[nodiscard]] constexpr T sum() const {
    return std::reduce(data_handle(), data_handle()+size(), T{});
  }
  
  [[nodiscard]] constexpr T operator*(const matpack_vector auto&x) const
  {
    ARTS_ASSERT(size() == x.size())
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }

  [[nodiscard]] constexpr auto transpose() const { return strided_view<T, N, true>{*this}.transpose(); }
  [[nodiscard]] constexpr auto transpose() { return strided_view<T, N, constant>{*this}.transpose(); }
  [[nodiscard]] constexpr auto real() const {return strided_view<T, N, true>{*this}.real();}
  [[nodiscard]] constexpr auto real() {return strided_view<T, N, constant>{*this}.real();}
  [[nodiscard]] constexpr auto imag() const {return strided_view<T, N, true>{*this}.imag();}
  [[nodiscard]] constexpr auto imag() {return strided_view<T, N, constant>{*this}.imag();}

  constexpr simple_view& operator/=(T x) requires(not constant) {
    std::transform(data_handle(), data_handle()+size(), data_handle(), [x](auto& y){return x / y;});
    return *this;
  }

  constexpr simple_view& operator*=(T x) requires(not constant) {
    std::transform(data_handle(), data_handle()+size(), data_handle(), [x](auto& y){return x * y;});
    return *this;
  }

  constexpr simple_view& operator+=(T x) requires(not constant) {
    std::transform(data_handle(), data_handle()+size(), data_handle(), [x](auto& y){return x + y;});
    return *this;
  }

  constexpr simple_view& operator-=(T x) requires(not constant) {
    std::transform(data_handle(), data_handle()+size(), data_handle(), [x](auto& y){return y - x;});
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator+=(const simple_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    std::transform(data_handle(), data_handle() + size(), x.data_handle(),
                   data_handle(), [](auto &a, auto &b) { return a + b; });
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator-=(const simple_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    std::transform(data_handle(), data_handle() + size(), x.data_handle(),
                   data_handle(), [](auto &a, auto &b) { return a - b; });
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator*=(const simple_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    std::transform(data_handle(), data_handle() + size(), x.data_handle(),
                   data_handle(), [](auto &a, auto &b) { return a * b; });
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator/=(const simple_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    std::transform(data_handle(), data_handle() + size(), x.data_handle(),
                   data_handle(), [](auto &a, auto &b) { return a / b; });
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator+=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    strided_view<T, N, false>{view} += x;
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator-=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    strided_view<T, N, false>{view} -= x;
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator*=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    strided_view<T, N, false>{view} *= x;
    return *this;
  }

  template <bool other_constant>
  constexpr simple_view &operator/=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    strided_view<T, N, false>{view} /= x;
    return *this;
  }

  template <class F>
  void transform_elementwise(F&& func) requires(not constant) {
    std::transform(data_handle(), data_handle()+size(), data_handle(), func);
  }
  [[nodiscard]] T max() const {return *std::max_element(data_handle(), data_handle()+size());}
  [[nodiscard]] T min() const {return *std::min_element(data_handle(), data_handle()+size());}

  [[nodiscard]] static constexpr bool is_exhaustive() noexcept {return true;}
  [[nodiscard]] T * data_handle() const noexcept {return view.data_handle();}
  [[nodiscard]] auto is_decreasing() const noexcept {return data_handle()+size() == std::adjacent_find(data_handle(), data_handle()+size(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a <= b;});}
  [[nodiscard]] auto is_increasing() const noexcept {return data_handle()+size() == std::adjacent_find(data_handle(), data_handle()+size(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a >= b;});}
  [[nodiscard]] auto is_sorted() const noexcept {return data_handle()+size() == std::adjacent_find(data_handle(), data_handle()+size(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a > b;});;}

  friend std::ostream &operator<<(std::ostream &os, const simple_view &sv) {
    return os << strided_view<T, N, true>{sv};
  }
};

template<math_type T, Index N, bool constant>
class strided_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::strided_mdspan<T, N> view;

  template <math_type U, Index M> friend class simple_data;
  template <math_type U, Index M, bool c> friend class simple_view;
  template <math_type U, Index M, bool c> friend class strided_view;

  template <Index dim>
  constexpr void strided_dimension(strided_access range)
    requires(dim >= 0 and dim < N)
  {
    range = strided_access(extent(dim), range);
    auto *d = view.data_handle() + stride(dim) * range.offset;
    std::array<Index, N> extmap, strmap;
    for (Index i = 0; i < N; i++) {
      extmap[i] = extent(i);
      strmap[i] = stride(i);
    }
    std::get<dim>(extmap) =
        range.extent < 0 ? std::get<dim>(extmap) : range.extent;
    std::get<dim>(strmap) *= range.stride;
    view = detail::strided_mdspan<T, N>{d, {extmap, strmap}};
  }

  template <Index i = 0, access_operator first,
            access_operator... access>
  constexpr strided_view &range_adaptor(first &&maybe_range, access &&...rest) {
    if constexpr (std::same_as<std::remove_cvref_t<first>, strided_access>)
      strided_dimension<i>(std::forward<first>(maybe_range));
    if constexpr (sizeof...(access) not_eq 0)
      range_adaptor<i + not integral<first>>(std::forward<access>(rest)...);
    return *this;
  }

public:
  constexpr strided_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {} 
  constexpr strided_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr strided_view(T& x) : view(detail::exhaustive_mdspan<T, N>{&x, ones<N>()}) {}

  template <matpack_convertible U>
  explicit strided_view(const U& v) :
    view(const_cast<T*>(mddata(v)), {mdshape(v), mdstride(v)}) {
  }

  template <matpack_convertible U>
  strided_view& operator=(const U& v) {return *this = strided_view(v);}

  constexpr strided_view() = default;
  constexpr strided_view(const strided_view&) = default;
  constexpr strided_view& operator=(const strided_view&) = default;
  constexpr strided_view(strided_view&&) noexcept = default;
  constexpr strided_view& operator=(strided_view&&) noexcept = default;

  explicit constexpr operator simple_view<T, N, false>() const requires(not constant) {return view;}
  explicit constexpr operator simple_view<T, N, true>() const requires(constant) {return view;}
  constexpr operator strided_view<T, N, true>() const requires(not constant) {return view;}
  template <Index M>
  explicit constexpr operator strided_view<T, M, true>() const requires(M > N) {
    std::array<Index, M> sz, st;
    for (Index i=0; i<M; i++) sz[i] = i < N ? extent(i) : 1;
    for (Index i=0; i<M; i++) st[i] = i < N ? stride(i) : 1;
    return detail::strided_mdspan<T, M>{view.data_handle(), {sz, st}};
  }
  template <Index M>
  explicit constexpr operator strided_view<T, M, false>() requires(M > N) {
    std::array<Index, M> sz, st;
    for (Index i=0; i<M; i++) sz[i] = i < N ? extent(i) : 1;
    for (Index i=0; i<M; i++) st[i] = i < N ? stride(i) : 1;
    return detail::strided_mdspan<T, M>{view.data_handle(), {sz, st}};
  }

  strided_view &operator=(const any_same_md<N> auto &sv)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == sv.shape(), "Mismatch shape")
    if (view.data_handle() not_eq sv.view.data_handle())
      std::copy(sv.begin(), sv.end(), begin());
    return *this;
  }

  constexpr strided_view& operator=(const T& x) requires(not constant) {
    if constexpr (N == 1)
      std::fill(begin(), end(), x);
    else
      for (auto y: *this) y = x;
    return *this;
  }

  [[nodiscard]] constexpr strided_view<T, N, true> transpose() const
    requires(N == 2)
  {
    using mapping_type = typename detail::strided_mdspan<T, N>::mapping_type;
    mapping_type transposed_map{
        std::array<Index, 2>{extent(1), extent(0)},
        std::array<Index, 2>{stride(1), stride(0)}};
    return detail::strided_mdspan<T, N>{view.data_handle(), transposed_map};
  }

  [[nodiscard]] constexpr strided_view transpose()
    requires(N == 2)
  {
    using mapping_type = typename detail::strided_mdspan<T, N>::mapping_type;
    mapping_type transposed_map{
        std::array<Index, 2>{extent(1), extent(0)},
        std::array<Index, 2>{stride(1), stride(0)}};
    return detail::strided_mdspan<T, N>{view.data_handle(), transposed_map};
  }

  [[nodiscard]] constexpr auto real() const
    requires(complex_type<T>)
  {
    using mapping_type =
        typename detail::strided_mdspan<complex_subtype<T>, N>::mapping_type;
    std::array<Index, N> extmap, strmap;
    for (Index i = 0; i < N; i++) {
      extmap[i] = extent(i);
      strmap[i] = 2 * stride(i);
    }
    return strided_view<complex_subtype<T>, N, true>{
        detail::strided_mdspan<complex_subtype<T>, N>{
            &reinterpret_cast<complex_subtype<T> *>(view.data_handle())[0],
            mapping_type{extmap, strmap}}};
  }

  [[nodiscard]] constexpr auto real()
    requires(complex_type<T>)
  {
    using mapping_type =
        typename detail::strided_mdspan<complex_subtype<T>, N>::mapping_type;
    std::array<Index, N> extmap, strmap;
    for (Index i = 0; i < N; i++) {
      extmap[i] = extent(i);
      strmap[i] = 2 * stride(i);
    }
    return strided_view<complex_subtype<T>, N, constant>{
        detail::strided_mdspan<complex_subtype<T>, N>{
            &reinterpret_cast<complex_subtype<T> *>(view.data_handle())[0],
            mapping_type{extmap, strmap}}};
  }

  [[nodiscard]] constexpr auto imag() const
    requires(complex_type<T>)
  {
    using mapping_type =
        typename detail::strided_mdspan<complex_subtype<T>, N>::mapping_type;
    std::array<Index, N> extmap, strmap;
    for (Index i = 0; i < N; i++) {
      extmap[i] = extent(i);
      strmap[i] = 2 * stride(i);
    }
    return strided_view<complex_subtype<T>, N, true>{
        detail::strided_mdspan<complex_subtype<T>, N>{
            &reinterpret_cast<complex_subtype<T> *>(view.data_handle())[1],
            mapping_type{extmap, strmap}}};
  }

  [[nodiscard]] constexpr auto imag()
    requires(complex_type<T>)
  {
    using mapping_type =
        typename detail::strided_mdspan<complex_subtype<T>, N>::mapping_type;
    std::array<Index, N> extmap, strmap;
    for (Index i = 0; i < N; i++) {
      extmap[i] = extent(i);
      strmap[i] = 2 * stride(i);
    }
    return strided_view<complex_subtype<T>, N, constant>{
        detail::strided_mdspan<complex_subtype<T>, N>{
            &reinterpret_cast<complex_subtype<T> *>(view.data_handle())[1],
            mapping_type{extmap, strmap}}};
  }

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES
  ITERATORS(strided_view, strided_view)

  //! ARTS Helper functions
  [[nodiscard]] constexpr T sum() const {
    if constexpr (N == 1) {
      return std::reduce(begin(), end(), T{});
    } else {
      return std::reduce(begin(), end(), T{}, [](auto& v) {return v.sum();});
    }
  }
  [[nodiscard]] constexpr T operator*(const matpack_vector auto &x) const requires(N == 1)
  {
    ARTS_ASSERT(shape() == x.shape())
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }

  constexpr strided_view& operator/=(T x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(begin(), end(), begin(), [x](auto& y){return x / y;});
    else
      for (auto y: *this) y /= x;
    return *this;
  }

  constexpr strided_view& operator*=(T x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(begin(), end(), begin(), [x](auto& y){return x * y;});
    else
      for (auto y: *this) y *= x;
    return *this;
  }

  constexpr strided_view &operator+=(T x)
    requires(not constant)
  {
    if constexpr (N == 1)
      std::transform(begin(), end(), begin(), [x](auto &y) { return x + y; });
    else
      for (auto y : *this)
        y += x;
    return *this;
  }

  constexpr strided_view &operator-=(T x)
    requires(not constant)
  {
    if constexpr (N == 1)
      std::transform(begin(), end(), begin(), [x](auto &y) { return y - x; });
    else
      for (auto y : *this)
        y -= x;
    return *this;
  }

  template <bool other_constant>
  constexpr strided_view &
  operator+=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    if constexpr (N == 1)
      std::transform(begin(), end(), x.begin(), begin(),
                     [](auto &a, auto &b) { return a + b; });
    else
      for (Index i = 0; i < extent(0); i++)
        operator[](i) += x[i];
    return *this;
  }

  template <bool other_constant>
  constexpr strided_view &
  operator-=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    if constexpr (N == 1)
      std::transform(begin(), end(), x.begin(), begin(),
                     [](auto &a, auto &b) { return a - b; });
    else
      for (Index i = 0; i < extent(0); i++)
        operator[](i) -= x[i];
    return *this;
  }

  template <bool other_constant>
  constexpr strided_view &
  operator*=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    ARTS_ASSERT(shape() == x.shape())
    if constexpr (N == 1)
      std::transform(begin(), end(), x.begin(), begin(),
                     [](auto &a, auto &b) { return a * b; });
    else
      for (Index i = 0; i < extent(0); i++)
        operator[](i) *= x[i];
    return *this;
  }

  template <bool other_constant>
  constexpr strided_view &
  operator/=(const strided_view<T, N, other_constant> &x)
    requires(not constant)
  {
    if constexpr (N == 1)
      std::transform(begin(), end(), x.begin(), begin(),
                     [](auto &a, auto &b) { return a / b; });
    else
      for (Index i = 0; i < extent(0); i++)
        operator[](i) /= x[i];
    return *this;
  }

  template <class F>
  void transform_elementwise(F&& func) requires(not constant) {
    if constexpr (N == 1)
      std::transform(begin(), end(), begin(), func);
    else
      for (auto sv: *this) sv.transform_elementwise(func);
  }
  [[nodiscard]] T max() const {
    if constexpr (N == 1)
      return *std::max_element(begin(), end());
    else {
      T out = std::numeric_limits<T>::lowest();
      for (auto &sv : *this)
        out = std::max(out, sv.max());
      return out;
    }
  }
  [[nodiscard]] T min() const {
    if constexpr (N == 1)
      return *std::min_element(begin(), end());
    else {
      T out = std::numeric_limits<T>::max();
      for (auto &sv : *this)
        out = std::min(out, sv.max());
      return out;
    }
  }

  [[nodiscard]] constexpr bool is_exhaustive() const noexcept {return view.is_exhaustive();}
  [[nodiscard]] auto data_handle() const noexcept {return simple_view<T, N, constant>{view}.data_handle();}
  [[nodiscard]] auto unsafe_data_handle() const noexcept {return view.data_handle();}
  [[nodiscard]] auto is_decreasing() const noexcept requires (N==1){return end() == std::adjacent_find(begin(), end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a <= b;});}
  [[nodiscard]] auto is_increasing() const noexcept requires (N==1) {return end() == std::adjacent_find(begin(), end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a >= b;});}
  [[nodiscard]] auto is_sorted() const noexcept requires (N==1) {return end() == std::adjacent_find(begin(), end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a > b;});}

  //! A common output operator
  friend std::ostream &operator<<(std::ostream &os, const strided_view &sv) {
    constexpr char space = N == 1 ? ' ' : '\n';

    bool first = true;
    for (auto&& v: sv) {
      if (not first) os << space; else first = false;
      os << std::forward<decltype(v)>(v);
    }
    return os;
  }
};

template <math_type T, Index N> class simple_data {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  static constexpr bool constant = false;

  std::vector<T> data;
  simple_view<T, N, constant> view;

  template <typename... arguments>
  static constexpr T defdata(arguments&&... args) noexcept
    requires(sizeof...(arguments) == N + 1)
  {
    return static_cast<T>(std::get<sizeof...(arguments) - 1>(
        std::make_tuple(std::forward<arguments>(args)...)));
  }

  template <typename... arguments>
  static constexpr void
  sz_data_impl(std::array<Index, N> &out,
                         integral auto &&x,
                         arguments&&... args [[maybe_unused]])
    requires(sizeof...(arguments) <= N)
  {
    std::get<N - sizeof...(arguments)>(out) =
        static_cast<Index>(std::forward<decltype(x)>(x));
    if constexpr (sizeof...(arguments) > 1)
      sz_data_impl(out, std::forward<arguments>(args)...);
  }

  template <typename... arguments>
  static constexpr std::array<Index, N>
  sz_data(arguments&&... args) noexcept
    requires(sizeof...(arguments) == N or sizeof...(arguments) == N + 1)
  {
    if constexpr (sizeof...(arguments) == N) {
      return {
          static_cast<Index>(std::forward<arguments>(args))...};
    } else {
      std::array<Index, N> out;
      sz_data_impl(out, std::forward<arguments>(args)...);
      return out;
    }
  }

  template <math_type U, Index M> friend class simple_data;
  template <math_type U, Index M, bool c> friend class simple_view;
  template <math_type U, Index M, bool c> friend class strided_view;

public:
  // Standard constructor
  constexpr simple_data(const std::array<Index, N> &sz = {},
                        const T &x = {})
      : data(std::reduce(sz.begin(), sz.end(), Index{1},
                         std::multiplies<>()),
             x),
        view(data.data(), sz) {}

  // Create the type completely default constructed
  template <integral... inds>
  constexpr simple_data(inds&&... ind)
    requires(sizeof...(inds) == N)
      : simple_data(sz_data(std::forward<inds>(ind)...)) {
  }

  // Create the type with known constant value
  template <typename... arguments>
  constexpr simple_data(arguments... args)
    requires(sizeof...(args) == N + 1)
      : simple_data(sz_data(std::forward<arguments>(args)...),
                    defdata(std::forward<arguments>(args)...)) {
  }

  // Create a linearly-spaced vector
  simple_data(T x0, Index sz, T dx) requires(N == 1) : data(sz), view(data.data(), sz) {
    auto gen = [x=x0, dx]() mutable {auto out=x; x+=dx; return out;};
    std::generate(data.begin(), data.end(), gen);
  }

  template <matpack_convertible U>
  explicit simple_data(const U& v) requires (N == rank<U>()) {
    *this = strided_view<T, N, true>(v);
  }

  constexpr simple_data(std::vector<T> &&x) noexcept
    requires(N == 1)
      : data(std::move(x)),
        view(data.data(), data.size()) {}

  simple_data(std::initializer_list<T> x) noexcept
    requires(N == 1) : simple_data(std::vector<T>{x}) {}

  constexpr simple_data(const simple_data& x) : data(x.data), view(data.data(), x.shape()) {}
  constexpr simple_data& operator=(const simple_data& x) { data=x.data; view=simple_view<T, N, false>{data.data(), x.shape()}; return *this; }
  constexpr simple_data(simple_data&&) noexcept = default;
  constexpr simple_data& operator=(simple_data&&) noexcept = default;

  template <matpack_convertible U>
  simple_data& operator=(const U& v) {return *this = strided_view<T, N, true>(v);}

  constexpr operator simple_view<T, N, false>&() {return view;}
  constexpr operator simple_view<T, N, true>() const {return view;}
  constexpr operator strided_view<T, N, false>() {return view;}
  constexpr operator strided_view<T, N, true>() const {return view;}
  template <Index M>
  explicit constexpr operator simple_view<T, M, true>() const requires(M > N) {
    return simple_view<T, M, true>{view};
  }
  template <Index M>
  explicit constexpr operator strided_view<T, M, true>() const requires(M > N) {
    return strided_view<T, M, true>{simple_view<T, M, true>{view}};
  }
  template <Index M>
  explicit constexpr operator simple_view<T, M, false>() requires(M > N) {
    return simple_view<T, M, false>{view};
  }
  template <Index M>
  explicit constexpr operator strided_view<T, M, false>() requires(M > N) {
    return strided_view<T, M, false>{simple_view<T, M, false>{view}};
  }

  constexpr simple_data(const any_same_view<N> auto &sv)
      : data(sv.size()), view(data.data(), sv.shape()) {
    std::copy(sv.begin(), sv.end(), begin());
  }

  constexpr simple_data& operator=(const any_same_view<N> auto& sv) {
    if (view.view.data_handle() not_eq sv.view.data_handle()) {
      if (auto s = sv.shape(); shape() not_eq s) resize(s);
      std::copy(sv.begin(), sv.end(), begin());
    }
    return *this;
  }

  // Resize operation
  constexpr void resize(std::array<Index, N> sz) {
    data.resize(std::reduce(sz.begin(), sz.end(), Index{1},
                            std::multiplies<>()));
    view = simple_view<T, N, false>{data.data(), sz};
  }

  // Resize operation
  template <integral... inds> constexpr void resize(inds &&...ind) {
    resize(std::array{Index(std::forward<inds>(ind))...});
  }

  // swap
  constexpr void swap(simple_data& other) noexcept {
    data.swap(other.data);
    std::swap(view, other.view);
  }

  // Reshape an rvalue
  template <integral... inds, Index M = sizeof...(inds)>
      [[nodiscard]] constexpr simple_data<T, M> reshape(inds... ind) &&
      requires(N not_eq M) {
        ARTS_ASSERT(static_cast<std::size_t>((ind * ...)) == view.size(),
                    "Mismatch in reshape")
        simple_data<T, M> out;
        out.data.swap(data);
        out.view = simple_view<T, M, false>{out.data.data(),
                                            std::forward<inds>(ind)...};
        view = simple_view<T, N, false>{data.data(),
                                        std::array<Index, N>{}};
        return out;
      }

  template <Index... inds, Index M = sizeof...(inds)>
  auto reduce_rank() &&
    requires(M < N)
  {
    std::array<Index, M> sz{inds...};
    for (auto &x : sz)
      x = extent(x);
    ARTS_ASSERT(
        size() == std::accumulate(sz.begin(), sz.end(), std::size_t{1}, std::multiplies<>()),
        size(), " is not ",
        std::accumulate(sz.begin(), sz.end(), std::size_t{1}, std::multiplies<>()))
    
    simple_data<T, M> out;
    out.data.swap(data);
    out.view = simple_view<T, M, false>{out.data.data(), sz};
    view =
        simple_view<T, N, false>{data.data(), std::array<Index, N>{}};
    return out;
  }

  // Flatten to a vector-like
  [[nodiscard]] constexpr simple_data<T, 1> flatten() && requires(N > 1) {
    return std::move(*this).reshape(data.size());
  }

  constexpr simple_data& operator=(const T& x) {
    std::fill(data.begin(), data.end(), x);
    return *this;
  }

  constexpr simple_data& operator+=(T x) {
    std::transform(data.begin(), data.end(), data.begin(), [x](auto& y){return x + y;});
    return *this;
  }

  constexpr simple_data& operator-=(T x) {
    std::transform(data.begin(), data.end(), data.begin(), [x](auto& y){return y - x;});
    return *this;
  }

  constexpr simple_data& operator/=(T x) {
    std::transform(data.begin(), data.end(), data.begin(), [x](auto& y){return x / y;});
    return *this;
  }

  constexpr simple_data& operator*=(T x) {
    std::transform(data.begin(), data.end(), data.begin(), [x](auto& y){return x * y;});
    return *this;
  }

  constexpr simple_data& operator+=(const simple_data& x) {ARTS_ASSERT(shape() == x.shape()) std::transform(data.begin(), data.end(), x.data.begin(), data.begin(), [](auto& a, auto& b){return a + b;}); return *this;}
  constexpr simple_data& operator-=(const simple_data& x) {ARTS_ASSERT(shape() == x.shape()) std::transform(data.begin(), data.end(), x.data.begin(), data.begin(), [](auto& a, auto& b){return a - b;}); return *this;}
  constexpr simple_data& operator*=(const simple_data& x) {ARTS_ASSERT(shape() == x.shape()) std::transform(data.begin(), data.end(), x.data.begin(), data.begin(), [](auto& a, auto& b){return a * b;}); return *this;}
  constexpr simple_data& operator/=(const simple_data& x) {ARTS_ASSERT(shape() == x.shape()) std::transform(data.begin(), data.end(), x.data.begin(), data.begin(), [](auto& a, auto& b){return a / b;}); return *this;}
  
  template <bool other_constant> constexpr simple_data& operator+=(const simple_view<T, N, other_constant>& x) {view += x; return *this;}
  template <bool other_constant> constexpr simple_data& operator-=(const simple_view<T, N, other_constant>& x) {view -= x; return *this;}
  template <bool other_constant> constexpr simple_data& operator*=(const simple_view<T, N, other_constant>& x) {view *= x; return *this;}
  template <bool other_constant> constexpr simple_data& operator/=(const simple_view<T, N, other_constant>& x) {view /= x; return *this;}

  template <bool other_constant> constexpr simple_data& operator-=(const strided_view<T, N, other_constant>& x) {view -= x; return *this;}
  template <bool other_constant> constexpr simple_data& operator+=(const strided_view<T, N, other_constant>& x) {view += x; return *this;}
  template <bool other_constant> constexpr simple_data& operator*=(const strided_view<T, N, other_constant>& x) {view *= x; return *this;}
  template <bool other_constant> constexpr simple_data& operator/=(const strided_view<T, N, other_constant>& x) {view /= x; return *this;}

  //! access operator
  ARTS_SIZES
  template<access_operator ... access> [[nodiscard]] constexpr auto operator()(access... ind) -> decltype(view(access{}...)) {return view(std::forward<access>(ind)...);}
  template<access_operator ... access> [[nodiscard]] constexpr auto operator()(access... ind) const -> decltype(view(access{}...)) {return view(std::forward<access>(ind)...);}
  [[nodiscard]] constexpr auto operator[](Index i) -> decltype(view[0]) {return view[i];}
  [[nodiscard]] constexpr auto operator[](Index i) const -> decltype(view[0]) {return view[i];}
  [[nodiscard]] constexpr auto operator[](strided_access i) {return view[i];}
  [[nodiscard]] constexpr auto operator[](strided_access i) const {return view[i];}
  using iterator = std::conditional_t<N==1, decltype(data.begin()), decltype(view.begin())>;
  using const_iterator = std::conditional_t<N==1, decltype(data.cbegin()), decltype(view.cbegin())>;
  using value_type = T;
  [[nodiscard]] static constexpr auto rank() { return N; }
  [[nodiscard]] static constexpr auto is_const() { return constant; }
  [[nodiscard]] iterator end() {if constexpr (N==1) return data.end(); else return view.end();}
  [[nodiscard]] const_iterator end() const {if constexpr (N==1) return data.end(); else return view.end();}
  [[nodiscard]] const_iterator cend() const {if constexpr (N==1) return data.cend(); else return view.cend();}
  [[nodiscard]] iterator begin() {if constexpr (N==1) return data.begin(); else return view.begin();}
  [[nodiscard]] const_iterator begin() const {if constexpr (N==1) return data.begin(); else return view.begin();}
  [[nodiscard]] const_iterator cbegin() const {if constexpr (N==1) return data.cbegin(); else return view.cbegin();}
  [[nodiscard]] auto rend() requires(N == 1) {return data.rend();}
  [[nodiscard]] auto rend() const requires(N == 1) {return data.rend();}
  [[nodiscard]] auto crend() const requires(N == 1) {return data.crend();}
  [[nodiscard]] auto rbegin() requires(N == 1) {return data.rbegin();}
  [[nodiscard]] auto rbegin() const requires(N == 1) {return data.rbegin();}
  [[nodiscard]] auto crbegin() const requires(N == 1) {return data.crbegin();}

  //! ARTS Helper functions
  [[nodiscard]] constexpr T operator*(const any_same_md<N> auto &x) const requires(N==1) {
    ARTS_ASSERT(size() == x.size())
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }
  [[nodiscard]] constexpr T sum() const {return std::reduce(data.begin(), data.end(), T{});}
  [[nodiscard]] constexpr auto transpose() const { return view.transpose(); }
  [[nodiscard]] constexpr auto transpose() { return view.transpose(); }
  [[nodiscard]] constexpr auto real() {return view.real();}
  [[nodiscard]] constexpr auto real() const {return view.real();}
  [[nodiscard]] constexpr auto imag() {return view.imag();}
  [[nodiscard]] constexpr auto imag() const {return view.imag();}

  template <class F>
  void transform_elementwise(F&& func) {
    std::transform(data.begin(), data.end(), data.begin(), func);
  }
  [[nodiscard]] T max() const {ARTS_ASSERT(size()) return *std::max_element(data.begin(), data.end());}
  [[nodiscard]] T min() const {ARTS_ASSERT(size()) return *std::min_element(data.begin(), data.end());}

  [[nodiscard]] static constexpr bool is_exhaustive() noexcept {return true;}
  [[nodiscard]] auto data_handle() const noexcept {return view.data_handle();}
  [[nodiscard]] auto is_decreasing() const noexcept {return data.end() == std::adjacent_find(data.begin(), data.end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a <= b;});}
  [[nodiscard]] auto is_increasing() const noexcept {return data.end() == std::adjacent_find(data.begin(), data.end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a >= b;});}
  [[nodiscard]] auto is_sorted() const noexcept {return data.end() == std::adjacent_find(data.begin(), data.end(), [](auto& a, auto& b){return std::isnan(a) or std::isnan(b) or a > b;});;}

  friend std::ostream &operator<<(std::ostream &os, const simple_data &sd) {
    return os << simple_view<T, N, true>{sd.view};
  }
};

#undef ACCESS_OPS
#undef ITERATION_OPS
#undef ARTS_SIZES

namespace static_tests {
template <typename T> consteval bool test_random_access_iterator() {
  static_assert(std::random_access_iterator<decltype(T{}.begin())>);
  static_assert(std::random_access_iterator<decltype(T{}.end())>);
  static_assert(std::random_access_iterator<decltype(T{}.cbegin())>);
  static_assert(std::random_access_iterator<decltype(T{}.cend())>);
  return true;
}

static_assert(test_random_access_iterator<std::vector<Numeric>>());
static_assert(test_random_access_iterator<simple_data<Numeric, 1>>());
static_assert(test_random_access_iterator<const simple_data<Numeric, 1>>());
static_assert(test_random_access_iterator<simple_data<Numeric, 2>>());
static_assert(test_random_access_iterator<const simple_data<Numeric, 2>>());
static_assert(test_random_access_iterator<simple_view<Numeric, 1, true>>());
static_assert(test_random_access_iterator<simple_view<Numeric, 1, false>>());
static_assert(test_random_access_iterator<const simple_view<Numeric, 1, true>>());
static_assert(test_random_access_iterator<const simple_view<Numeric, 1, false>>());
static_assert(test_random_access_iterator<simple_view<Numeric, 2, true>>());
static_assert(test_random_access_iterator<simple_view<Numeric, 2, false>>());
static_assert(test_random_access_iterator<const simple_view<Numeric, 2, true>>());
static_assert(test_random_access_iterator<const simple_view<Numeric, 2, false>>());
static_assert(test_random_access_iterator<strided_view<Numeric, 1, true>>());
static_assert(test_random_access_iterator<strided_view<Numeric, 1, false>>());
static_assert(test_random_access_iterator<const strided_view<Numeric, 1, true>>());
static_assert(test_random_access_iterator<const strided_view<Numeric, 1, false>>());
static_assert(test_random_access_iterator<strided_view<Numeric, 2, true>>());
static_assert(test_random_access_iterator<strided_view<Numeric, 2, false>>());
static_assert(test_random_access_iterator<const strided_view<Numeric, 2, true>>());
static_assert(test_random_access_iterator<const strided_view<Numeric, 2, false>>());
} // namespace static_tests
}  // namespace matpack::md

using Range = matpack::md::strided_access;
inline constexpr matpack::md::Joker joker = matpack::md::joker;

// Vectors
using Vector = matpack::md::simple_data<Numeric, 1>;
using FastVectorView = matpack::md::simple_view<Numeric, 1, false>;
using FastConstVectorView = matpack::md::simple_view<Numeric, 1, true>;
using VectorView = matpack::md::strided_view<Numeric, 1, false>;
using ConstVectorView = matpack::md::strided_view<Numeric, 1, true>;

// Matrices
using Matrix = matpack::md::simple_data<Numeric, 2>;
using FastMatrixView = matpack::md::simple_view<Numeric, 2, false>;
using FastConstMatrixView = matpack::md::simple_view<Numeric, 2, true>;
using MatrixView = matpack::md::strided_view<Numeric, 2, false>;
using ConstMatrixView = matpack::md::strided_view<Numeric, 2, true>;

// Tensor3(s)
using Tensor3 = matpack::md::simple_data<Numeric, 3>;
using FastTensor3View = matpack::md::simple_view<Numeric, 3, false>;
using FastConstTensor3View = matpack::md::simple_view<Numeric, 3, true>;
using Tensor3View = matpack::md::strided_view<Numeric, 3, false>;
using ConstTensor3View = matpack::md::strided_view<Numeric, 3, true>;

// Tensor4(s)
using Tensor4 = matpack::md::simple_data<Numeric, 4>;
using FastTensor4View = matpack::md::simple_view<Numeric, 4, false>;
using FastConstTensor4View = matpack::md::simple_view<Numeric, 4, true>;
using Tensor4View = matpack::md::strided_view<Numeric, 4, false>;
using ConstTensor4View = matpack::md::strided_view<Numeric, 4, true>;

// Tensor5(s)
using Tensor5 = matpack::md::simple_data<Numeric, 5>;
using FastTensor5View = matpack::md::simple_view<Numeric, 5, false>;
using FastConstTensor5View = matpack::md::simple_view<Numeric, 5, true>;
using Tensor5View = matpack::md::strided_view<Numeric, 5, false>;
using ConstTensor5View = matpack::md::strided_view<Numeric, 5, true>;

// Tensor6(s)
using Tensor6 = matpack::md::simple_data<Numeric, 6>;
using FastTensor6View = matpack::md::simple_view<Numeric, 6, false>;
using FastConstTensor6View = matpack::md::simple_view<Numeric, 6, true>;
using Tensor6View = matpack::md::strided_view<Numeric, 6, false>;
using ConstTensor6View = matpack::md::strided_view<Numeric, 6, true>;

// Tensor7(s)
using Tensor7 = matpack::md::simple_data<Numeric, 7>;
using FastTensor7View = matpack::md::simple_view<Numeric, 7, false>;
using FastConstTensor7View = matpack::md::simple_view<Numeric, 7, true>;
using Tensor7View = matpack::md::strided_view<Numeric, 7, false>;
using ConstTensor7View = matpack::md::strided_view<Numeric, 7, true>;

// Complex vectors
using ComplexVector = matpack::md::simple_data<Complex, 1>;
using FastComplexVectorView = matpack::md::simple_view<Complex, 1, false>;
using FastConstComplexVectorView = matpack::md::simple_view<Complex, 1, true>;
using ComplexVectorView = matpack::md::strided_view<Complex, 1, false>;
using ConstComplexVectorView = matpack::md::strided_view<Complex, 1, true>;

// Complex matrices
using ComplexMatrix = matpack::md::simple_data<Complex, 2>;
using FastComplexMatrixView = matpack::md::simple_view<Complex, 2, false>;
using FastConstComplexMatrixView = matpack::md::simple_view<Complex, 2, true>;
using ComplexMatrixView = matpack::md::strided_view<Complex, 2, false>;
using ConstComplexMatrixView = matpack::md::strided_view<Complex, 2, true>;

/** Compute y = alpha * A * x + beta * y
 *
 * Uses LAPACK for speedy computations
 *
 * @param[inout] y A vector view that has the size of N
 * @param[in] A A matrix view that has the size of M x N
 * @param[in] x A vector view that has the size of M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult_fast(matpack::md::simple_view<Numeric, 1, false> y,
               const matpack::md::simple_view<Numeric, 2, true> &A,
               const matpack::md::simple_view<Numeric, 1, true> &x,
               Numeric alpha = 1.0, Numeric beta = 0.0);

/** Compute y = A * x
 *
 * @param[out] y A vector view that has the size of N
 * @param[in] A A matrix view that has the size of M x N
 * @param[in] x A vector view that has the size of M
 */
void mult(matpack::md::strided_view<Numeric, 1, false> y,
          const matpack::md::strided_view<Numeric, 2, true> &A,
          const matpack::md::strided_view<Numeric, 1, true> &x);

/** Compute C = alpha * A * B + beta * C
 *
 * Uses LAPACK for speedy computations
 *
 * @param[inout] C A matrix view that has the size of N x M
 * @param[in] A A matrix view that has the size of N x K
 * @param[in] B A matrix view that has the size of K x M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult_fast(matpack::md::simple_view<Numeric, 2, false> C,
               const matpack::md::simple_view<Numeric, 2, true> &A,
               const matpack::md::simple_view<Numeric, 2, true> &B,
               Numeric alpha = 1.0, Numeric beta = 0.0);

/** Compute C = alpha * A * B + beta * C
 *
 * @param[out] C A matrix view that has the size of N x M
 * @param[in] A A matrix view that has the size of N x K
 * @param[in] B A matrix view that has the size of K x M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult(matpack::md::strided_view<Numeric, 2, false> C,
          const matpack::md::strided_view<Numeric, 2, true> &A,
          const matpack::md::strided_view<Numeric, 2, true> &B,
          Numeric alpha = 1.0, Numeric beta = 0.0);

/** Compute y = alpha * A * x + beta * y
 *
 * Uses LAPACK for speedy computations
 *
 * @param[inout] y A vector view that has the size of N
 * @param[in] A A matrix view that has the size of M x N
 * @param[in] x A vector view that has the size of M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult_fast(matpack::md::simple_view<Complex, 1, false> y,
               const matpack::md::simple_view<Complex, 2, true> &A,
               const matpack::md::simple_view<Complex, 1, true> &x,
               Complex alpha = 1.0, Complex beta = 0.0);

/** Compute y = A * x
 *
 * @param[out] y A vector view that has the size of N
 * @param[in] A A matrix view that has the size of M x N
 * @param[in] x A vector view that has the size of M
 */
void mult(matpack::md::strided_view<Complex, 1, false> y,
          const matpack::md::strided_view<Complex, 2, true> &A,
          const matpack::md::strided_view<Complex, 1, true> &x);

/** Compute C = alpha * A * B + beta * C
 *
 * Uses LAPACK for speedy computations
 *
 * @param[inout] C A matrix view that has the size of N x M
 * @param[in] A A matrix view that has the size of N x K
 * @param[in] B A matrix view that has the size of K x M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult_fast(matpack::md::simple_view<Complex, 2, false> C,
               const matpack::md::simple_view<Complex, 2, true> &A,
               const matpack::md::simple_view<Complex, 2, true> &B,
               Complex alpha = 1.0, Complex beta = 0.0);

/** Compute C = alpha * A * B + beta * C
 *
 * @param[out] C A matrix view that has the size of N x M
 * @param[in] A A matrix view that has the size of N x K
 * @param[in] B A matrix view that has the size of K x M
 * @param[in] alpha A constant (default 1.0)
 * @param[in] beta A constant (default 0.0)
 */
void mult(matpack::md::strided_view<Complex, 2, false> C,
          const matpack::md::strided_view<Complex, 2, true> &A,
          const matpack::md::strided_view<Complex, 2, true> &B,
          Complex alpha = 1.0, Complex beta = 0.0);

namespace detail {
//! Help deduce if we can transpose a matrix
template<typename T>
concept transposable_matrix = requires(T a) {
  { a.transpose() } -> matpack::md::matpack_matrix;
};
}  // namespace detail

/** Returns a view of the transpose of some object
 *
 * @param[in] A A matrix view that has the size of N x K
 * @return A const-correct matrix view of size K x N
 */
constexpr auto transpose(detail::transposable_matrix auto &&A) noexcept
{
  return std::forward<decltype(A)>(A).transpose();
}

namespace detail {
//! Implement a transform functionality
template <typename T, Index N, typename Function>
void transform_impl(matpack::md::any_mutable<T, N> auto& out, Function&& f, const matpack::md::any_md<T, N> auto& in) {
  out = in;
  out.transform_elementwise(std::forward<Function>(f));
}

//! A concept to say wheter or not we can call transform_impl
template <typename T>
concept transform_implable = requires(T a) {
  { a.rank() } -> std::integral;
  { a.is_const() } -> std::same_as<bool>;
};
} // namespace detail


/** Transforms the out variable by applying f on all in
 * 
 * @tparam OUTPUT An outputable multi-dimensional array
 * @tparam INPUT An inputable multi-dimensional array
 * @tparam Function Some function 
 */
template <typename OUTPUT, typename INPUT, typename Function>
void transform(OUTPUT &out, Function &&f, const INPUT &in)
  requires(
      detail::transform_implable<OUTPUT> and
      detail::transform_implable<INPUT> and
      not OUTPUT::is_const() and
      INPUT::rank() == OUTPUT::rank() and
      std::same_as<typename INPUT::value_type, typename OUTPUT::value_type> and
      requires(OUTPUT a, Function g) {a.transform_elementwise(g);})
{
  detail::transform_impl<typename INPUT::value_type, INPUT::rank()>(
      out, std::forward<Function>(f), in);
}

/** Creates a linearly spaced vector of size n between start and stop
 * 
 * @param out Vector of some type
 * @param start First value
 * @param stop Last value
 * @param n Number of steps
 */
template <typename T>
void nlinspace(matpack::md::simple_data<T, 1>& out, auto start, auto stop, Index n) {
  ARTS_ASSERT(1 < n)
  out.resize(n);
  T step = T(stop - start) / T(n - 1);
  for (Index i = 0; i < n - 1; i++) out[i] = start + (T)i * step;
  out[n - 1] = stop;
}

/** Creates a linearly spaced vector of size n between start and stop
 * 
 * @param out Vector of some type
 * @param start First value
 * @param stop Last value
 * @param n Number of steps
 */
template <typename T>
void linspace(matpack::md::simple_data<T, 1>& out, auto start, auto stop, auto step_) {
  Index n = (Index)floor((stop - start) / step_) + 1;
  if (n < 1) n = 1;
  out.resize(n);
  T step = T(step_);
  for (Index i = 0; i < n - 1; i++) out[i] = start + (T)i * step;
  out[n - 1] = stop;
}

namespace detail {
//! The type has max()
template <typename T>
concept maxable = requires(T a) {
  a.max();
};

//! The type has min()
template <typename T>
concept minable = requires(T a) {
  a.min();
};
}  // namespace detail

/** Returns the maximum
 * 
 * @param in Some type that has a max() member function
 * @return auto The return of said max() member function call
 */
auto max(const detail::maxable auto& in) {
  return in.max();
}

/** Returns the minimum
 * 
 * @param in Some type that has a min() member function
 * @return auto The return of said max() member function call
 */
auto min(const detail::minable auto& in) {
  return in.min();
}

/** Checks if the list is sorted ascendingly (i.e., always increases or remain the same)
 * 
 * @param a Some list
 * @return true if sorted correctly
 */
template <typename T>
auto is_sorted(const T& a) requires(requires(T x) {x.is_sorted();}) {
  return a.is_sorted();
}

/** Checks if the list is ascending always
 * 
 * @param a Some list
 * @return true if sorted correctly
 */
template <typename T>
auto is_increasing(const T& a) requires(requires(T x) {x.is_increasing();}) {
  return a.is_increasing();
}

/** Checks if the list is sorted decreasing always
 * 
 * @param a Some list
 * @return true if sorted correctly
 */
template <typename T>
auto is_decreasing(const T& a) requires(requires(T x) {x.is_decreasing();}) {
  return a.is_decreasing();
}