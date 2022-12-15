#include "debug.h"
#include "matpack.h"

#include <algorithm>
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

template <detail::size_t Rank>
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
template <typename T, detail::size_t N> class simple_data;
template <typename T, detail::size_t N, bool constant> class simple_view;
template <typename T, detail::size_t N, bool constant> class strided_view;
template <typename T, detail::size_t N, typename U>
concept mutable_view =
    std::same_as<std::remove_cvref_t<U>, simple_view<T, N, false>> or
    std::same_as<std::remove_cvref_t<U>, strided_view<T, N, false>>;
template <typename T, detail::size_t N, typename U>
concept const_view =
    std::same_as<std::remove_cvref_t<U>, simple_view<T, N, true>> or
    std::same_as<std::remove_cvref_t<U>, strided_view<T, N, true>>;
template <typename T, detail::size_t N, typename U>
concept any_view =
    mutable_view<T, N, U> or const_view<T, N, U>;
template <typename T, detail::size_t N, typename U>
concept any_md = any_view<T, N, U> or std::same_as<std::remove_cvref_t<U>, simple_data<T, N>>;

template <access_operator Access, typename... arguments>
[[nodiscard]] consteval detail::size_t get_access_size() noexcept {
  constexpr bool index = std::integral<Access>;
  if constexpr (sizeof...(arguments)) {
    return index + get_access_size<arguments...>();
  } else {
    return index;
  }
}

template <access_operator T, access_operator ... Ts>
[[nodiscard]] bool check_index_sizes(auto&& arr, detail::size_t r, T first, Ts... rest) {
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

template <detail::size_t N>
[[nodiscard]] consteval std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  // Is this necessary? out.fill(joker);
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
  [[nodiscard]] auto operator()(access... ind)                                 \
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
  [[nodiscard]] auto operator()(access... ind) const->std::conditional_t<      \
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
  [[nodiscard]] auto operator[](detail::size_t i)                              \
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
  [[nodiscard]] auto operator[](detail::size_t i) const->std::conditional_t<   \
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
  [[nodiscard]] auto size() const noexcept { return view.size(); }             \
  [[nodiscard]] auto extent(detail::size_t i) const noexcept {                 \
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
    return extent(0);                                                          \
  }                                                                            \
  [[nodiscard]] auto nrows() const noexcept                                    \
    requires(N >= 2)                                                           \
  {                                                                            \
    return extent(1);                                                          \
  }                                                                            \
  [[nodiscard]] auto npages() const noexcept                                   \
    requires(N >= 3)                                                           \
  {                                                                            \
    return extent(2);                                                          \
  }                                                                            \
  [[nodiscard]] auto nbooks() const noexcept                                   \
    requires(N >= 4)                                                           \
  {                                                                            \
    return extent(3);                                                          \
  }                                                                            \
  [[nodiscard]] auto nshelves() const noexcept                                 \
    requires(N >= 5)                                                           \
  {                                                                            \
    return extent(4);                                                          \
  }                                                                            \
  [[nodiscard]] auto nvitrines() const noexcept                                \
    requires(N >= 6)                                                           \
  {                                                                            \
    return extent(5);                                                          \
  }                                                                            \
  [[nodiscard]] auto nlibraries() const noexcept                               \
    requires(N >= 7)                                                           \
  {                                                                            \
    return extent(6);                                                          \
  }                                                                            \
  [[nodiscard]] std::array<detail::size_t, N> shape() const noexcept {         \
    std::array<detail::size_t, N> out;                                         \
    for (detail::size_t i = 0; i < N; i++)                                     \
      out[i] = extent(i);                                                      \
    return out;                                                                \
  }                                                                            \
  [[nodiscard]] bool empty() const { return view.empty(); }

#define ITERATORS(IN, OUT)                                                     \
  using iterator =                                                             \
      mditer<0, false, IN, OUT<T, std::max<detail::size_t>(N - 1, 1), false>>; \
  using const_iterator =                                                       \
      mditer<0, true, const IN,                                                \
             const OUT<T, std::max<detail::size_t>(N - 1, 1), true>>;          \
  using value_type = T;                                                        \
  [[nodiscard]] static constexpr auto rank() { return N; }                     \
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

template <detail::size_t M, bool constant, class mdspan_type,
          class submdspan_type>
class mditer {
  static constexpr auto N = mdspan_type::rank();
  using T = std::conditional_t<constant, const typename mdspan_type::value_type,
                               typename mdspan_type::value_type>;

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

  constexpr mditer& operator++() noexcept {pos++; return *this;}
  constexpr mditer& operator--() noexcept {pos--; return *this;}
  [[nodiscard]] constexpr mditer operator++(int) const noexcept {mditer out(*this); ++out.pos; return out;}
  [[nodiscard]] constexpr mditer operator--(int) const noexcept {mditer out(*this); --out.pos; return out;}

  constexpr mditer& operator+=(detail::size_t i) noexcept {pos+=i; return *this;}
  constexpr mditer& operator-=(detail::size_t i) noexcept {pos-=i; return *this;}
  [[nodiscard]] constexpr mditer operator+(detail::size_t i) const noexcept {mditer out(*this); out.pos+=i; return out;}
  [[nodiscard]] constexpr mditer operator-(detail::size_t i) const noexcept {mditer out(*this); out.pos-=i; return out;}
  [[nodiscard]] constexpr friend mditer operator+(detail::size_t i, const mditer& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend mditer operator-(detail::size_t i, const mditer& m) noexcept {return m - i;}

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

  [[nodiscard]] auto operator[](detail::size_t i) const
      -> std::conditional_t<N == 1, std::add_lvalue_reference_t<T>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->operator[](pos + i);
    } else {
      return sub<M, N>(*orig, pos + i);
    }
  }
};

template <typename T, detail::size_t N, bool constant> class simple_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::exhaustive_mdspan<T, N> view;

  template <typename U, detail::size_t M> friend class simple_data;
  template <typename U, detail::size_t M, bool c> friend class simple_view;
  template <typename U, detail::size_t M, bool c> friend class strided_view;

public:
  constexpr simple_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}
  constexpr simple_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {}

  constexpr simple_view() = default;
  constexpr simple_view(const simple_view&) = default;
  constexpr simple_view& operator=(const simple_view&) = default;
  constexpr simple_view(simple_view&&) noexcept = default;
  constexpr simple_view& operator=(simple_view&&) noexcept = default;

  constexpr operator simple_view<T, N, true>() const requires(not constant) {return view;}
  constexpr operator strided_view<T, N, constant>() {return view;}
  constexpr operator strided_view<T, N, true>() const {return view;}

  template <typename U>
  constexpr simple_view &operator=(const U &sv)
    requires(not constant and any_md<T, N, U>)
  {
    ARTS_ASSERT(shape() == sv.shape(), "Mismatch shape")
    if (view.data_handle() not_eq sv.view.data_handle())
      std::copy(sv.begin(), sv.end(), begin());
    return *this;
  }

  constexpr simple_view& operator=(const T& x) {
    std::fill(begin(), end(), x);
    return *this;
  }

  //! access operator
  ACCESS_OPS
  ITERATION_OPS
  ARTS_SIZES
  ITERATORS(simple_view, simple_view)

  //! ARTS Helper functions

  //! ARTS Helper functions
  [[nodiscard]] constexpr T sum() const {
    if constexpr (N == 1) {
      return std::reduce(begin(), end(), T{});
    } else {
      return std::reduce(begin(), end(), T{}, [](auto& v) {return v.sum();});
    }
  }
  template <typename U>
  [[nodiscard]] constexpr T operator*(const U&x) const
    requires(N == 1 and any_md<T, 1, U>)
  {
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }

  [[nodiscard]] static constexpr bool is_exhaustive() noexcept {return true;}

  friend std::ostream &operator<<(std::ostream &os, const simple_view &sv) {
    return os << strided_view<T, N, true>{sv};
  }
};

template<typename T, detail::size_t N, bool constant>
class strided_view {
  static_assert(N >= 1, "Must be vector-like or of greater rank");

  detail::strided_mdspan<T, N> view;

  template <typename U, detail::size_t M> friend class simple_data;
  template <typename U, detail::size_t M, bool c> friend class simple_view;
  template <typename U, detail::size_t M, bool c> friend class strided_view;

public:
  constexpr strided_view(detail::strided_mdspan<T, N> v) : view(std::move(v)) {} 
  constexpr strided_view(detail::exhaustive_mdspan<T, N> v) : view(std::move(v)) {}

  constexpr strided_view() = default;
  constexpr strided_view(const strided_view&) = default;
  constexpr strided_view& operator=(const strided_view&) = default;
  constexpr strided_view(strided_view&&) noexcept = default;
  constexpr strided_view& operator=(strided_view&&) noexcept = default;

  constexpr operator strided_view<T, N, true>() const requires(not constant) {return view;}

  template <typename U>
  strided_view &operator=(const U &sv)
    requires(not constant and any_md<T, N, U>)
  {
    ARTS_ASSERT(shape() == sv.shape(), "Mismatch shape")
    if (view.data_handle() not_eq sv.view.data_handle())
      std::copy(sv.begin(), sv.end(), begin());
    return *this;
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
  template <typename U>
  [[nodiscard]] constexpr T operator*(const U&x) const
    requires(N == 1 and any_md<T, 1, U>)
  {
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }

  [[nodiscard]] constexpr bool is_exhaustive() const noexcept {return view.is_exhaustive();}

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

template <typename T, detail::size_t N> class simple_data {
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

  template <typename U, detail::size_t M> friend class simple_data;
  template <typename U, detail::size_t M, bool c> friend class simple_view;
  template <typename U, detail::size_t M, bool c> friend class strided_view;

public:
  // Standard constructor
  constexpr simple_data(const std::array<detail::size_t, N> &sz = {},
                        const T &x = {})
      : data(std::reduce(sz.begin(), sz.end(), detail::size_t{1},
                         std::multiplies<>()),
             x),
        view(detail::exhaustive_mdspan<T, N>{data.data(), sz}) {}

  // Create the type completely default constructed
  template <std::integral... inds>
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

  constexpr simple_data(std::vector<T> &&x) noexcept
    requires(N == 1)
      : data(std::move(x)),
        view(detail::exhaustive_mdspan<T, N>{data.data(), data.size()}) {}

  constexpr simple_data(const simple_data& x) : data(x.data), view(data.data(), x.shape()) {}
  constexpr simple_data& operator=(const simple_data& x) { data=x.data; view = detail::exhaustive_mdspan<T, N>{data.data(), x.shape()}; return *this; }
  constexpr simple_data(simple_data&&) noexcept = default;
  constexpr simple_data& operator=(simple_data&&) noexcept = default;

  constexpr operator simple_view<T, N, false>&() {return view;}
  constexpr operator simple_view<T, N, true>() const {return view;}
  constexpr operator strided_view<T, N, false>() {return view;}
  constexpr operator strided_view<T, N, true>() const {return view;}

  // Resize operation
  constexpr void resize(std::array<detail::size_t, N> sz) {
    data.resize(std::reduce(sz.begin(), sz.end(), detail::size_t{1},
                            std::multiplies<>()));
    view = detail::exhaustive_mdspan<T, N>{data.data(), sz};
  }

  // Resize operation
  template <std::integral... inds> constexpr void resize(inds &&...ind) {
    resize(std::array{detail::size_t(std::forward<inds>(ind))...});
  }

  // swap
  constexpr void swap(simple_data& other) noexcept {
    data.swap(other.data);
    std::swap(view, other.view);
  }

  // Reshape an rvalue
  template <std::integral... inds, detail::size_t M = sizeof...(inds)>
  [[nodiscard]] constexpr simple_data<T, M> reshape(inds... ind) && requires(N not_eq M) {
    ARTS_ASSERT(static_cast<std::size_t>((ind * ...)) == view.size(),
                "Mismatch in reshape")
    simple_data<T, M> out;
    out.data.swap(data);
    out.view = detail::exhaustive_mdspan<T, M>{out.data.data(), std::forward<inds>(ind)...};
    view = detail::exhaustive_mdspan<T, N>{data.data(), std::array<detail::size_t, N>{}};
    return out;
  }

  // Flatten to a vector-like
  [[nodiscard]] constexpr simple_data<T, 1> flatten() && requires(N > 1) {
    return std::move(*this).reshape(data.size());
  }

  template <typename U>
  constexpr simple_data& operator=(const U& sv) requires(any_view<T, N, U>) {
    if (view.view.data_handle() not_eq sv.view.data_handle()) {
      if (auto s = sv.shape(); shape() not_eq s) resize(s);
      std::copy(sv.begin(), sv.end(), begin());
    }
    return *this;
  }

  constexpr simple_data& operator=(const T& x) {
    std::fill(data.begin(), data.end(), x);
    return *this;
  }

  //! access operator
  ARTS_SIZES
  template<access_operator ... access> [[nodiscard]] constexpr auto operator()(access&&... ind) -> decltype(view(access{}...)) {return view(std::forward<access>(ind)...);}
  template<access_operator ... access> [[nodiscard]] constexpr auto operator()(access&&... ind) const -> decltype(view(access{}...)) {return view(std::forward<access>(ind)...);}
  [[nodiscard]] constexpr auto operator[](detail::size_t i) -> decltype(view[0]) {return view[i];}
  [[nodiscard]] constexpr auto operator[](detail::size_t i) const -> decltype(view[0]) {return view[i];}
  using iterator = std::conditional_t<N==1, decltype(data.begin()), decltype(view.begin())>;
  using const_iterator = std::conditional_t<N==1, decltype(data.cbegin()), decltype(view.cbegin())>;
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
  [[nodiscard]] constexpr T sum() const {return std::reduce(data.begin(), data.end(), T{});}
  template <typename U>
  [[nodiscard]] constexpr T operator*(const U &x) const
    requires(N == 1 and any_md<T, 1, U>)
  {
    return std::transform_reduce(begin(), end(), x.begin(), T{});
  }

  [[nodiscard]] static constexpr bool is_exhaustive() noexcept {return true;}

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

using TMPVector = matpack::md::simple_data<Numeric, 1>;
using TMPMatrix = matpack::md::simple_data<Numeric, 2>;
using TMPTensor3 = matpack::md::simple_data<Numeric, 3>;

using Vec = matpack::md::simple_data<Numeric, 1>;
using Mat = matpack::md::simple_data<Numeric, 2>;
using Ten3 = matpack::md::simple_data<Numeric, 3>;

using VecVF = matpack::md::simple_view<Numeric, 1, false>;
using MatVF = matpack::md::simple_view<Numeric, 2, false>;
using Ten3VF = matpack::md::simple_view<Numeric, 3, false>;
using VecCVF = matpack::md::simple_view<Numeric, 1, true>;
using MatCVF = matpack::md::simple_view<Numeric, 2, true>;
using Ten3CVF = matpack::md::simple_view<Numeric, 3, true>;

using VecVS = matpack::md::strided_view<Numeric, 1, false>;
using MatVS = matpack::md::strided_view<Numeric, 2, false>;
using Ten3VS = matpack::md::strided_view<Numeric, 3, false>;
using VecCVS = matpack::md::strided_view<Numeric, 1, true>;
using MatCVS = matpack::md::strided_view<Numeric, 2, true>;
using Ten3CVS = matpack::md::strided_view<Numeric, 3, true>;

