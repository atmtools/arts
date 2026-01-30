#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <span>
#include <type_traits>
#include <vector>

#include "matpack_mdspan_common_types.h"
#include "matpack_mdspan_elemwise_mditer.h"
#include "matpack_mdspan_strided_view_t.h"

namespace matpack {
//! A wrapper for mdspan that's always C-like in layout
template <class T, Size N>
struct view_t final : public mdview_t<T, N> {
  constexpr static bool matpack_magic_view = true;

  using base = mdview_t<T, N>;

  constexpr mdview_t<T, N> base_md() const { return *this; }

  using extents_type     = base::extents_type;
  using layout_type      = base::layout_type;
  using accessor_type    = base::accessor_type;
  using mapping_type     = base::mapping_type;
  using element_type     = base::element_type;
  using value_type       = base::value_type;
  using index_type       = base::index_type;
  using size_type        = base::size_type;
  using rank_type        = base::rank_type;
  using data_handle_type = base::data_handle_type;
  using reference        = base::reference;

  using base::is_always_exhaustive;
  using base::is_always_strided;
  using base::is_always_unique;
  using base::rank;
  using base::rank_dynamic;
  using base::static_extent;

  using base::accessor;
  using base::data_handle;
  using base::empty;
  using base::extent;
  using base::extents;
  using base::mapping;
  using base::size;
  using base::stride;

  constexpr void base_set(base v) noexcept { static_cast<base&>(*this) = v; }

  static constexpr bool is_const      = std::is_const_v<T>;
  static constexpr bool is_exhaustive = is_always_exhaustive();
  static constexpr bool is_strided    = is_always_strided();
  static constexpr bool is_unique     = is_always_unique();

  constexpr view_t() : base() {};
  constexpr view_t(const view_t& x) : base(x) {};
  constexpr view_t(view_t&& x) noexcept : base(std::move(x)) {};

  //! Create from ANY mdspan - it is your responsibility to ensure that this is runtime valid
  template <class OtherElementType,
            class OtherExtents,
            class OtherLayoutPolicy,
            class OtherAccessor>
  explicit constexpr view_t(const stdx::mdspan<OtherElementType,
                                               OtherExtents,
                                               OtherLayoutPolicy,
                                               OtherAccessor>& other)
      : base(other) {}

  constexpr view_t(const data_t<value_type, N>& x)
    requires(is_const)
      : base(x.base_md()) {}

  constexpr view_t(const view_t<value_type, N>& x)
    requires(is_const)
      : base(x.base_md()) {}

  explicit constexpr view_t(T& v) : base(&v, std::array<Index, 1>{1}) {}
  explicit constexpr view_t(const T& v)
    requires(is_const)
      : base(&v, std::array<Index, 1>{1}) {}

  explicit constexpr view_t(std::vector<value_type>& v)
    requires(N == 1)
      : base(const_cast<value_type*>(v.data()), v.size()) {}

  explicit constexpr view_t(const std::vector<value_type>& v)
    requires(N == 1 and is_const)
      : base(const_cast<value_type*>(v.data()), v.size()) {}

  [[nodiscard]] constexpr auto shape() const {
    std::array<Index, N> out;
    for (Size i = 0; i < N; i++) out[i] = this->extent(i);
    return out;
  }

  template <typename Self, access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](this Self&& self,
                                                    Acc&&... i)
    requires(sizeof...(Acc) <= N)
  {
    using U = decltype(std::forward<Self>(self));

    if constexpr (sizeof...(Acc) == N and (integral<Acc> and ...)) {
      if constexpr (const_forwarded<U>)
        return std::as_const(
            std::forward<Self>(self).base::operator[](std::forward<Acc>(i)...));
      else
        return std::forward<Self>(self).base::operator[](
            std::forward<Acc>(i)...);
    } else {
      if constexpr (const_forwarded<U>)
        return left_mdsel_t<std::add_const_t<T>, N, is_exhaustive, Acc...>(
            left_sub<N>(std::forward<Self>(self), std::forward<Acc>(i)...));
      else
        return left_mdsel_t<T, N, is_exhaustive, Acc...>(
            left_sub<N>(std::forward<Self>(self), std::forward<Acc>(i)...));
    }
  }

  template <typename Self>
  [[nodiscard]] constexpr auto begin(this Self&& self) {
    if constexpr (N == 1) {
      return std::forward<Self>(self).elem_begin();
    } else {
      return left_mditer(self);
    }
  }

  template <typename Self>
  [[nodiscard]] constexpr auto end(this Self&& self) {
    if constexpr (N == 1) {
      return std::forward<Self>(self).elem_end();
    } else {
      const auto n = self.extent(0);
      return std::forward<Self>(self).begin() + n;
    }
  }

  template <typename Self, Size M>
  constexpr auto view_as(this Self&& self, const std::array<Index, M>& exts) {
    assert(self.size() == mdsize(exts));
    return view_t<
        std::conditional_t<const_forwarded<Self>, std::add_const_t<T>, T>,
        M>{mdview_t<T, M>{self.data_handle(), exts}};
  }

  template <typename Self, integral... NewExtents>
  constexpr auto view_as(this Self&& self, NewExtents... exts) {
    return std::forward<Self>(self).view_as(
        std::array{static_cast<Index>(exts)...});
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_begin(this Self&& self) {
    using U = decltype(std::forward<Self>(self));

    if constexpr (const_forwarded<U>)
      return const_cast<const T*>(std::forward<Self>(self).data_handle());
    else
      return std::forward<Self>(self).data_handle();
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_end(this Self&& self) {
    const auto n = self.size();
    return std::forward<Self>(self).elem_begin() + n;
  }

  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) elem_at(this Self&& self, Size i) {
    assert(self.size() > i);
    return *(std::forward<Self>(self).elem_begin() + i);
  }

  constexpr view_t& operator=(const view_t& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if (this != &r) stdr::copy(elemwise_range(r), elem_begin());

    return *this;
  }

  constexpr view_t& operator=(const value_type& x)
    requires(not is_const)
  {
    std::fill(elem_begin(), elem_end(), x);
    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr view_t& operator=(const R& r)
    requires(not(is_const or std::same_as<view_t, R> or
                 std::same_as<value_type, R>))
  {
    assert(shape() == r.shape());

    if constexpr (std::same_as<T, matpack::value_type<R>>) {
      stdr::copy(elemwise_range(r), elem_begin());
    } else {
      stdr::transform(elemwise_range(r), elem_begin(), [](const auto& a) -> T {
        return a;
      });
    }

    return *this;
  }

  template <mdvalue_type_compatible<T> U>
  constexpr view_t& operator=(const U& x)
    requires(not(is_const or std::same_as<view_t, U> or
                 std::same_as<value_type, U> or ranked_convertible_md<U, T, N>))
  {
    const auto sh = shape();

    assert(sh == mdshape(x));

    for (Size i = 0; i < size(); i++)
      elem_at(i) = static_cast<T>(mdvalue(x, sh, i));

    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr view_t& operator+=(const R& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if constexpr (std::same_as<T, matpack::value_type<R>>) {
      stdr::transform(elemwise_range(*this),
                      elemwise_range(r),
                      elem_begin(),
                      [](auto&& a, auto&& b) -> T { return a + b; });
    } else {
      stdr::transform(
          elemwise_range(*this),
          elemwise_range(r),
          elem_begin(),
          [](auto&& a, auto&& b) -> T { return a + static_cast<T>(b); });
    }

    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr view_t& operator-=(const R& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if constexpr (std::same_as<T, matpack::value_type<R>>) {
      stdr::transform(elemwise_range(*this),
                      elemwise_range(r),
                      elem_begin(),
                      [](auto&& a, auto&& b) -> T { return a - b; });
    } else {
      stdr::transform(
          elemwise_range(*this),
          elemwise_range(r),
          elem_begin(),
          [](auto&& a, auto&& b) -> T { return a - static_cast<T>(b); });
    }

    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr view_t& operator*=(const R& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if constexpr (std::same_as<T, matpack::value_type<R>>) {
      stdr::transform(elemwise_range(*this),
                      elemwise_range(r),
                      elem_begin(),
                      [](auto&& a, auto&& b) -> T { return a * b; });
    } else {
      stdr::transform(
          elemwise_range(*this),
          elemwise_range(r),
          elem_begin(),
          [](auto&& a, auto&& b) -> T { return a * static_cast<T>(b); });
    }

    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr view_t& operator/=(const R& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());
    if constexpr (std::same_as<T, matpack::value_type<R>>) {
      stdr::transform(elemwise_range(*this),
                      elemwise_range(r),
                      elem_begin(),
                      [](auto&& a, auto&& b) -> T { return a / b; });
    } else {
      stdr::transform(
          elemwise_range(*this),
          elemwise_range(r),
          elem_begin(),
          [](auto&& a, auto&& b) -> T { return a / static_cast<T>(b); });
    }
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator+=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x += v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator-=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x -= v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator*=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x *= v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator/=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x /= v; });
    return *this;
  }

  template <typename Self>
  constexpr auto real(this Self&& self) {
    return strided_view_t<T, N>(self).real();
  }

  template <typename Self>
  constexpr auto imag(this Self&& self) {
    return strided_view_t<T, N>(self).imag();
  }

  //! std::vector-like
  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) front(this Self&& self) {
    assert(not self.empty());
    return std::forward<Self>(self).elem_at(0);
  }

  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) back(this Self&& self) {
    assert(not self.empty());
    return std::forward<Self>(self).elem_at(self.size() - 1);
  }

  //! Sizes
  [[nodiscard]] constexpr auto ncols() const { return extent(N - 1); }
  [[nodiscard]] constexpr auto nrows() const { return extent(N - 2); }
  [[nodiscard]] constexpr auto npages() const { return extent(N - 3); }
  [[nodiscard]] constexpr auto nbooks() const { return extent(N - 4); }
  [[nodiscard]] constexpr auto nshelves() const { return extent(N - 5); }
  [[nodiscard]] constexpr auto nvitrines() const { return extent(N - 6); }
  [[nodiscard]] constexpr auto nlibraries() const { return extent(N - 7); }

  constexpr void swap(view_t& other) noexcept
    requires(not is_const)
  {
    std::swap(static_cast<base&>(*this), static_cast<base&>(other));
  }

  //! Used by std::indirectly_writable concept, remove comment if linter works
  const view_t& operator=(const view_t& r) const
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if (this != &r) {
      stdr::copy(elemwise_range(const_cast<view_t&>(*this)),
                 elemwise_range(const_cast<view_t&>(r)));
    }

    return *this;
  }

  friend constexpr void swap(const view_t& a, const view_t& b) noexcept
    requires(not is_const)
  {
    assert(a.shape() == b.shape());

    if (a.data_handle() != b.data_handle()) {
      stdr::swap_ranges(elemwise_range(const_cast<view_t&>(a)),
                        elemwise_range(const_cast<view_t&>(b)));
    }
  }

  constexpr operator std::span<T, std::dynamic_extent>()
    requires(N == 1)
  {
    return {begin(), end()};
  }
  constexpr operator std::span<const T, std::dynamic_extent>() const
    requires(N == 1)
  {
    return {begin(), end()};
  }
};
}  // namespace matpack

template <class T, Size N>
inline constexpr bool
    std::ranges::enable_borrowed_range<matpack::view_t<T, N>> = true;

std::string to_string(const matpack::view_t<const Numeric, 2>&,
                      format_tags& tags,
                      const std::span<const Size>);
std::string to_string(const matpack::view_t<const Complex, 2>&,
                      format_tags& tags,
                      const std::span<const Size>);

std::string to_string(const matpack::view_t<Numeric, 2>&,
                      format_tags& tags,
                      const std::span<const Size>);
std::string to_string(const matpack::view_t<Complex, 2>&,
                      format_tags& tags,
                      const std::span<const Size>);

template <typename T, Size N>
struct std::formatter<matpack::view_t<T, N>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    tags.depth = N;
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::view_t<T, N>& md,
                              FmtContext& ctx) const {
    if (md.empty()) return ctx.out();

    auto mat_view =
        md.view_as(md.size() / md.shape().back(), md.shape().back());

    if constexpr (requires { to_string(mat_view, tags); }) {
      auto nl = md.shape();
      for (Size i = N - 2; i < N; i--) nl[i] *= nl[i + 1];
      tags.format(ctx, to_string(mat_view, tags, nl));
    } else {
      std::formatter<matpack::strided_view_t<T, N>> fmt;
      fmt.tags = tags;
      fmt.format(md, ctx);
    }

    return ctx.out();
  }
};

using VectorView  = matpack::view_t<Numeric, 1>;
using MatrixView  = matpack::view_t<Numeric, 2>;
using Tensor3View = matpack::view_t<Numeric, 3>;
using Tensor4View = matpack::view_t<Numeric, 4>;
using Tensor5View = matpack::view_t<Numeric, 5>;
using Tensor6View = matpack::view_t<Numeric, 6>;
using Tensor7View = matpack::view_t<Numeric, 7>;

using ComplexVectorView  = matpack::view_t<Complex, 1>;
using ComplexMatrixView  = matpack::view_t<Complex, 2>;
using ComplexTensor3View = matpack::view_t<Complex, 3>;
using ComplexTensor4View = matpack::view_t<Complex, 4>;
using ComplexTensor5View = matpack::view_t<Complex, 5>;
using ComplexTensor6View = matpack::view_t<Complex, 6>;
using ComplexTensor7View = matpack::view_t<Complex, 7>;

using ConstVectorView  = matpack::view_t<const Numeric, 1>;
using ConstMatrixView  = matpack::view_t<const Numeric, 2>;
using ConstTensor3View = matpack::view_t<const Numeric, 3>;
using ConstTensor4View = matpack::view_t<const Numeric, 4>;
using ConstTensor5View = matpack::view_t<const Numeric, 5>;
using ConstTensor6View = matpack::view_t<const Numeric, 6>;
using ConstTensor7View = matpack::view_t<const Numeric, 7>;

using ConstComplexVectorView  = matpack::view_t<const Complex, 1>;
using ConstComplexMatrixView  = matpack::view_t<const Complex, 2>;
using ConstComplexTensor3View = matpack::view_t<const Complex, 3>;
using ConstComplexTensor4View = matpack::view_t<const Complex, 4>;
using ConstComplexTensor5View = matpack::view_t<const Complex, 5>;
using ConstComplexTensor6View = matpack::view_t<const Complex, 6>;
using ConstComplexTensor7View = matpack::view_t<const Complex, 7>;
