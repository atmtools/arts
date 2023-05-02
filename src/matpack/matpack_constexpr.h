#pragma once

#include "matpack_concepts.h"
#include "matpack_iter.h"

#include <algorithm>
#include <concepts>
#include <type_traits>
#include <utility>

namespace matpack {
//! The constant matpack view
template <typename T, bool constant, Index... alldim>
struct matpack_constant_view {
  static constexpr Index N = sizeof...(alldim);
  static constexpr Index sz = (alldim * ...);
  static_assert(sz > 0);

  template <typename U, Index R>
  static constexpr bool magic_test_of_matpack_constant_view() {
    return N == R and std::same_as<U, T>;
  }

  using view_type = stdx::mdspan<T, stdx::extents<Index, alldim...>>;
  using data_type = matpack_constant_data<T, alldim...>;
  using const_matpack_constant_view = matpack_constant_view<T, true, alldim...>;
  using mut_matpack_constant_view = matpack_constant_view<T, false, alldim...>;

  view_type view;

  template <bool c, Index outer, Index... rest> struct inner_type {
    using type = matpack_constant_view<T, c, rest...>;
  };

  [[nodiscard]] constexpr matpack_constant_view(std::array<T, sz> &x) : view(x.data()) {}
  [[nodiscard]] constexpr matpack_constant_view(const std::array<T, sz> &x) requires(constant) : view(const_cast<T*>(x.data())) {}
  [[nodiscard]] constexpr matpack_constant_view(const view_type &x) : view(x) {}
  [[nodiscard]] constexpr matpack_constant_view(const mut_matpack_constant_view& x) requires(constant) : view(x.view) {}

  template <bool c> using inner_view = typename inner_type<c, alldim...>::type;

  template <integral... index>
  [[nodiscard]] constexpr T &operator()(index... i)
    requires(not constant and sizeof...(index) == N)
  {
    ARTS_ASSERT((std::cmp_less(static_cast<Index>(i), alldim) && ...),
                shape_help<N>{std::array{static_cast<Index>(i)...}}, " vs ",
                shape_help<N>{shape()})
    return view(i...);
  }

  template <integral... index>
  [[nodiscard]] constexpr const T &operator()(index... i) const
    requires(sizeof...(index) == N)
  {
    ARTS_ASSERT((std::cmp_less(static_cast<Index>(i), alldim) && ...),
                shape_help<N>{std::array{static_cast<Index>(i)...}}, " vs ",
                shape_help<N>{shape()})
    return view(i...);
  }

  template <typename... Jokers>
  [[nodiscard]] constexpr inner_view<false> operator()(Index i, Jokers... v)
    requires(not constant and sizeof...(Jokers) == N - 1 and N > 1 and
             (std::same_as<std::remove_cvref_t<Jokers>, Joker> and ...))
  {
    ARTS_ASSERT(i < view.extent(0), " vs ", extent(0))
    return stdx::submdspan(view, i, v...);
  }

  template <typename... Jokers>
  [[nodiscard]] constexpr inner_view<true> operator()(Index i, Jokers... v) const
    requires(sizeof...(Jokers) == N - 1 and N > 1 and
             (std::same_as<std::remove_cvref_t<Jokers>, Joker> and ...))
  {
    ARTS_ASSERT(i < extent(0), i, " vs ", extent(0))
    return stdx::submdspan(view, i, v...);
  }

  [[nodiscard]] constexpr auto operator[](Index i)
      -> std::conditional_t<N == 1, T &, inner_view<false>>
    requires(not constant)
  {
    ARTS_ASSERT(i < extent(0), " vs ", extent(0))
    if constexpr (N == 1)
      return view[i];
    else
      return sub<0, N>(view, i);
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<N == 1, const T &, inner_view<true>> {
    ARTS_ASSERT(i < extent(0), " vs ", extent(0))
    if constexpr (N == 1)
      return view[i];
    else
      return sub<0, N>(view, i);
  }

  using value_type = T;
  [[nodiscard]] static constexpr auto shape() { return std::array{alldim...}; }
  [[nodiscard]] constexpr auto extent(Index i) const { return view.extent(i); }
  [[nodiscard]] constexpr auto stride(Index i) const { return view.stride(i); }
  [[nodiscard]] static constexpr auto size() { return sz; }
  [[nodiscard]] static constexpr auto rank() { return N; }
  [[nodiscard]] static constexpr auto is_const() { return constant; }

  using iterator = matpack_mditer<0, false, matpack_constant_view, inner_view<false>>;
  using const_iterator = matpack_mditer<0, true, const matpack_constant_view, const inner_view<true>>;

  [[nodiscard]] constexpr auto begin() requires(not constant) { if constexpr (N == 1) return view.data_handle(); else return iterator{*this}; }
  [[nodiscard]] constexpr auto end() requires(not constant) { if constexpr (N == 1) return view.data_handle() + sz; else return iterator{*this} + extent(0); }
  [[nodiscard]] constexpr auto begin() const { if constexpr (N == 1) return view.data_handle(); else return const_iterator{*this}; }
  [[nodiscard]] constexpr auto end() const { if constexpr (N == 1) return view.data_handle() + sz; else return const_iterator{*this} + extent(0); }

  [[nodiscard]] constexpr auto elem_begin() requires(not constant) { return view.data_handle(); }
  [[nodiscard]] constexpr auto elem_end() requires(not constant) { return view.data_handle() + sz; }
  [[nodiscard]] constexpr auto elem_begin() const { return view.data_handle(); }
  [[nodiscard]] constexpr auto elem_end() const { return view.data_handle() + sz; }

  [[nodiscard]] constexpr bool operator==(const data_type &other) const {
    return view == other.view();
  }
  [[nodiscard]] constexpr bool operator!=(const data_type &other) const {
    return view != other.view();
  }
  [[nodiscard]] constexpr bool operator==(const mut_matpack_constant_view &other) const {
    return std::equal(elem_begin(), elem_end(), other.elem_begin());
  }
  [[nodiscard]] constexpr bool operator!=(const mut_matpack_constant_view &other) const {
    return not std::equal(elem_begin(), elem_end(), other.elem_begin());
  }
  [[nodiscard]] constexpr bool operator==(const const_matpack_constant_view &other) const {
    return std::equal(elem_begin(), elem_end(), other.elem_begin());
  }
  [[nodiscard]] constexpr bool operator!=(const const_matpack_constant_view &other) const {
    return not std::equal(elem_begin(), elem_end(), other.elem_begin());
  }

  friend std::ostream &operator<<(std::ostream &os, const matpack_constant_view &mv) {
    constexpr char extra = N == 1 ? ' ' : '\n';
    bool first = true;
    for (auto &&v : mv) {
      if (not first)
        os << extra;
      else
        first = false;
      os << std::forward<decltype(v)>(v);
    }
    return os;
  }
};

//! The constant matpack data
template <typename T, Index... alldim> struct matpack_constant_data {
  static constexpr Index N = sizeof...(alldim);
  static constexpr Index sz = (alldim * ...);
  static_assert(sz > 0);

  template <typename U, Index R>
  static constexpr bool magic_test_of_matpack_constant_data() {
    return N == R and std::same_as<U, T>;
  }

  using mut_view_type = matpack_constant_view<T, false, alldim...>;
  using const_view_type = matpack_constant_view<T, true, alldim...>;
  template <bool c> using inner_view = typename mut_view_type::template inner_view<c>;

  std::array<T, sz> data{};

  [[nodiscard]] explicit constexpr operator matpack_data<T, N>() const {
    if constexpr (N == 1) return matpack_data<T, 1>{data};
    else return matpack_data<T, 1>{data}.reshape(alldim...);
  }

  [[nodiscard]] constexpr auto view() { return mut_view_type{data}; }
  [[nodiscard]] constexpr auto view() const { return const_view_type{data}; }

  template <integral... index>
  [[nodiscard]] constexpr T &operator()(index... i)
    requires(sizeof...(index) == N)
  {
    return view()(i...);
  }

  template <integral... index>
  [[nodiscard]] constexpr T operator()(index... i) const
    requires(sizeof...(index) == N)
  {
    return view()(i...);
  }

  template <typename... Jokers>
  [[nodiscard]] constexpr inner_view<false> operator()(Index i, Jokers... v)
    requires(sizeof...(Jokers) == N - 1 and N > 1 and
             (std::same_as<std::remove_cvref_t<Jokers>, Joker> and ...))
  {
    return view()(i, v...);
  }

  template <typename... Jokers>
  [[nodiscard]] constexpr inner_view<true> operator()(Index i,
                                                      Jokers... v) const
    requires(sizeof...(Jokers) == N - 1 and N > 1 and
             (std::same_as<std::remove_cvref_t<Jokers>, Joker> and ...))
  {
    return view()(i, v...);
  }

  [[nodiscard]] constexpr auto operator[](Index i)
      -> std::conditional_t<N == 1, T &, inner_view<false>> {
    return view()[i];
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<N == 1, const T &, inner_view<true>> {
    return view()[i];
  }

  using value_type = T;
  [[nodiscard]] static constexpr auto shape() { return std::array{alldim...}; }
  [[nodiscard]] constexpr auto extent(Index i) const { return view().extent(i); }
  [[nodiscard]] constexpr auto stride(Index i) const { return view().stride(i); }
  [[nodiscard]] static constexpr auto size() { return sz; }
  [[nodiscard]] static constexpr auto rank() { return N; }

  [[nodiscard]] constexpr auto begin() { if constexpr (N == 1) return data.begin(); else return view().begin(); }
  [[nodiscard]] constexpr auto end() { if constexpr (N == 1) return data.end(); else return view().end(); }
  [[nodiscard]] constexpr auto begin() const { if constexpr (N == 1) return data.begin(); else return view().begin(); }
  [[nodiscard]] constexpr auto end() const { if constexpr (N == 1) return data.end(); else return view().end(); }

  [[nodiscard]] constexpr auto elem_begin() { return data.begin(); }
  [[nodiscard]] constexpr auto elem_end() { return data.end(); }
  [[nodiscard]] constexpr auto elem_begin() const { return data.begin(); }
  [[nodiscard]] constexpr auto elem_end() const { return data.end(); }

  [[nodiscard]] constexpr bool operator==(const matpack_constant_data &other) const { return data == other.data; }
  [[nodiscard]] constexpr bool operator!=(const matpack_constant_data &other) const { return data != other.data; }
  [[nodiscard]] constexpr bool operator==(const mut_view_type &other) const { return view() == other; }
  [[nodiscard]] constexpr bool operator!=(const mut_view_type &other) const { return view() != other; }
  [[nodiscard]] constexpr bool operator==(const const_view_type &other) const { return view() == other; }
  [[nodiscard]] constexpr bool operator!=(const const_view_type &other) const { return view() != other; }

  matpack_constant_data& operator+=(const matpack_constant_data& o) {
    std::transform(elem_begin(), elem_end(), o.elem_begin(), elem_begin(), [](auto a, auto b){return a + b;});
    return *this;
  }
  matpack_constant_data& operator+=(T x) {
    std::transform(elem_begin(), elem_end(), elem_begin(), [b=x](auto a){return a + b;});
    return *this;
  }
  matpack_constant_data& operator-=(const matpack_constant_data& o) {
    std::transform(elem_begin(), elem_end(), o.elem_begin(), elem_begin(), [](auto a, auto b){return a - b;});
    return *this;
  }
  matpack_constant_data& operator-=(T x) {
    std::transform(elem_begin(), elem_end(), elem_begin(), [b=x](auto a){return a - b;});
    return *this;
  }
  matpack_constant_data& operator/=(const matpack_constant_data& o) {
    std::transform(elem_begin(), elem_end(), o.elem_begin(), elem_begin(), [](auto a, auto b){return a / b;});
    return *this;
  }
  matpack_constant_data& operator/=(T x) {
    std::transform(elem_begin(), elem_end(), elem_begin(), [b=x](auto a){return a / b;});
    return *this;
  }
  matpack_constant_data& operator*=(const matpack_constant_data& o) {
    std::transform(elem_begin(), elem_end(), o.elem_begin(), elem_begin(), [](auto a, auto b){return a * b;});
    return *this;
  }
  matpack_constant_data& operator*=(T x) {
    std::transform(elem_begin(), elem_end(), elem_begin(), [b=x](auto a){return a * b;});
    return *this;
  }

  template <Index i> constexpr T& get() & {return std::get<i>(data);}

  template <Index i> constexpr const T& get() const & {return std::get<i>(data);}

  template <Index i> constexpr T&& get() && {return std::get<i>(data);}

  friend std::ostream &operator<<(std::ostream &os, const matpack_constant_data &mv) {
    return os << mv.view();
  }
};
} // namespace matpack

using Vector2 = matpack::matpack_constant_data<Numeric, 2>;
using Vector3 = matpack::matpack_constant_data<Numeric, 3>;

//! Make the constant data structured, so "[a,b,c] = Vector3{1,2,3};" works
namespace std {
  template <matpack::any_matpack_constant_data T> struct tuple_size<T> : std::integral_constant<size_t, T::size()> { };
  template <std::size_t I, matpack::any_matpack_constant_data T> struct tuple_element<I, T> { using type = typename T::value_type; };
}  // namespace std
