#pragma once

#include <concepts>
#include <functional>
#include <type_traits>

#include "experimental/__p0009_bits/extents.hpp"
#include "matpack_mdspan_common.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_elemwise_mditer.h"

namespace matpack {
template <typename T, Size... dims>
struct [[nodiscard]] cdata_t {
  constexpr static bool matpack_magic_cdata = true;

  constexpr static Size N     = sizeof...(dims);
  constexpr static Size ndata = (dims * ...);

  using extents_type     = view_t<T, N>::extents_type;
  using layout_type      = view_t<T, N>::layout_type;
  using accessor_type    = view_t<T, N>::accessor_type;
  using mapping_type     = view_t<T, N>::mapping_type;
  using element_type     = view_t<T, N>::element_type;
  using value_type       = view_t<T, N>::value_type;
  using index_type       = view_t<T, N>::index_type;
  using size_type        = view_t<T, N>::size_type;
  using rank_type        = view_t<T, N>::rank_type;
  using data_handle_type = view_t<T, N>::data_handle_type;
  using reference        = view_t<T, N>::reference;

  static_assert(not std::is_const_v<T>);
  static_assert(ndata > 0);

  std::array<T, ndata> data;

  static constexpr std::array<Index, N> shape_ = {dims...};
  static constexpr const std::array<Index, N>& shape() { return shape_; }

  constexpr stdx::mdspan<const T, stdx::extents<Size, dims...>> cview() const {
    return stdx::mdspan<const T, stdx::extents<Size, dims...>>{
        const_cast<T*>(data.data())};
  }

  constexpr stdx::mdspan<T, stdx::extents<Size, dims...>> cview() {
    return stdx::mdspan<T, stdx::extents<Size, dims...>>{
        const_cast<T*>(data.data())};
  }

  constexpr view_t<const T, N> view() const {
    return mdview_t<T, N>{const_cast<T*>(data.data()), shape_};
  }

  constexpr view_t<T, N> view() {
    return mdview_t<T, N>{const_cast<T*>(data.data()), shape_};
  }

  constexpr auto base_md() const { return cview(); }
  constexpr auto base_md() { return cview(); }

  constexpr operator view_t<T, N>() { return view(); }
  constexpr operator view_t<const T, N>() const { return view(); }
  constexpr operator strided_view_t<T, N>() { return view(); }
  constexpr operator strided_view_t<const T, N>() const { return view(); }
  explicit constexpr operator data_t<T, N>() const {
    return data_t<T, N>{view()};
  }

  /*

   Common operators using both data and view

   */

  constexpr void swap(cdata_t& other) noexcept { data.swap(other.data); }

  /*

  Operations purely on the view

   */

  static constexpr auto rank() { return view_t<T, N>::rank(); }
  static constexpr auto rank_dynamic() { return view_t<T, N>::rank_dynamic(); }
  static constexpr auto static_extent() {
    return view_t<T, N>::static_extent();
  }

  static constexpr bool is_const      = false;
  static constexpr bool is_exhaustive = view_t<T, N>::is_always_exhaustive();
  static constexpr bool is_strided    = view_t<T, N>::is_always_strided();
  static constexpr bool is_unique     = view_t<T, N>::is_always_unique();

  //! Sizes
  [[nodiscard]] constexpr auto ncols() const { return extent(N - 1); }
  [[nodiscard]] constexpr auto nrows() const { return extent(N - 2); }
  [[nodiscard]] constexpr auto npages() const { return extent(N - 3); }
  [[nodiscard]] constexpr auto nbooks() const { return extent(N - 4); }
  [[nodiscard]] constexpr auto nshelves() const { return extent(N - 5); }
  [[nodiscard]] constexpr auto nvitrines() const { return extent(N - 6); }
  [[nodiscard]] constexpr auto nlibraries() const { return extent(N - 7); }

  template <typename Self, Index M>
  constexpr auto view_as(this Self&& self, const std::array<Index, M>& exts) {
    return std::forward<Self>(self).view().view_as(exts);
  }

  template <typename Self, integral... inds>
  constexpr auto view_as(this Self&& self, inds... exts) {
    return std::forward<Self>(self).view().view_as(std::forward<inds>(exts)...);
  }

  template <access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](Acc&&... i)
    requires(sizeof...(Acc) <= N)
  {
    if constexpr (sizeof...(Acc) == N and (integral<Acc> and ...)) {
      return cview()[std::forward<Acc>(i)...];
    } else {
      return view()[std::forward<Acc>(i)...];
    }
  }

  template <access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](Acc&&... i) const
    requires(sizeof...(Acc) <= N)
  {
    if constexpr (sizeof...(Acc) == N and (integral<Acc> and ...)) {
      return cview()[std::forward<Acc>(i)...];
    } else {
      return view()[std::forward<Acc>(i)...];
    }
  }

  template <typename U>
  constexpr cdata_t& operator+=(U&& x) {
    view() += std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr cdata_t& operator-=(U&& x) {
    view() -= std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr cdata_t& operator*=(U&& x) {
    view() *= std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr cdata_t& operator/=(U&& x) {
    view() /= std::forward<U>(x);
    return *this;
  }

  template <typename Self>
  constexpr decltype(auto) elem_at(this Self&& self, Index i) {
    return std::forward<Self>(self).data[i];
  }

  template <typename Self>
  constexpr auto elem_begin(this Self&& self) {
    return std::forward<Self>(self).data.begin();
  }

  template <typename Self>
  constexpr auto elem_end(this Self&& self) {
    return std::forward<Self>(self).data.end();
  }

  template <typename Self>
  constexpr auto begin(this Self&& self) {
    return std::forward<Self>(self).view().begin();
  }

  template <typename Self>
  constexpr auto end(this Self&& self) {
    return std::forward<Self>(self).view().end();
  }

  template <typename Self>
  constexpr auto imag(this Self&& self) {
    return std::forward<Self>(self).view().imag();
  }

  template <typename Self>
  constexpr auto real(this Self&& self) {
    return std::forward<Self>(self).view().real();
  }

  template <typename Self>
  constexpr auto accessor(this Self&& self) {
    return std::forward<Self>(self).view().accessor();
  }

  template <typename Self>
  constexpr auto data_handle(this Self&& self) {
    return std::forward<Self>(self).data.data();
  }

  template <typename Self>
  constexpr auto empty(this Self&& self) {
    return std::forward<Self>(self).view().empty();
  }

  template <typename Self>
  constexpr auto extent(this Self&& self, Index i) {
    return std::forward<Self>(self).view().extent(i);
  }

  template <typename Self>
  constexpr auto extents(this Self&& self) {
    return std::forward<Self>(self).view().extents();
  }

  template <typename Self>
  constexpr auto mapping(this Self&& self) {
    return std::forward<Self>(self).view().mapping();
  }

  template <typename Self>
  constexpr auto size(this Self&& self) {
    return std::forward<Self>(self).view().size();
  }

  template <typename Self>
  constexpr auto stride(this Self&& self, integral auto i) {
    return std::forward<Self>(self).view().stride(i);
  }

  //! Tuple interface:
  template <Index i>
  constexpr T& get() & {
    return std::get<i>(data);
  }

  template <Index i>
  constexpr const T& get() const& {
    return std::get<i>(data);
  }

  template <Index i>
  constexpr T&& get() && {
    return std::get<i>(data);
  }
};

template <any_cdata T>
constexpr T operator+(const T& x, const T& y) {
  T z  = x;
  z   += y;
  return z;
}

template <any_cdata T>
constexpr T operator-(const T& x, const T& y) {
  T z  = x;
  z   -= y;
  return z;
}

template <any_cdata T>
constexpr T normalized(T x) {
  x /= hypot(x);
  return x;
}
}  // namespace matpack

template <typename T, Size... dims>
struct std::formatter<matpack::cdata_t<T, dims...>> {
  std::formatter<matpack::view_t<const T, sizeof...(dims)>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::cdata_t<T, dims...>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
  }
};

//! Make the constant data structured, so "[a,b,c] = Vector3{1,2,3};" works
namespace std {
template <matpack::any_cdata T>
struct tuple_size<T> : std::integral_constant<Size, T::ndata> {};

template <std::size_t I, matpack::any_cdata T>
struct tuple_element<I, T> {
  using type = matpack::value_type<T>;
};
}  // namespace std

using Vector2 = matpack::cdata_t<Numeric, 2>;
using Vector3 = matpack::cdata_t<Numeric, 3>;
using Vector4 = matpack::cdata_t<Numeric, 4>;
