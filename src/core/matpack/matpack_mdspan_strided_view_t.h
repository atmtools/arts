#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <type_traits>

#include "matpack_mdspan_common.h"
#include "matpack_mdspan_mditer.h"

namespace matpack {
template <class T, Size N>
struct strided_view_t final : public mdstrided_t<T, N> {
  constexpr static bool matpack_magic_strided_view = true;

  using base                          = mdstrided_t<T, N>;

  constexpr mdstrided_t<T, N> base_md() const { return *this; }

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

  // Basic copy/move/default
  constexpr strided_view_t() : base() {};
  constexpr strided_view_t(const strided_view_t& x) : base(x) {};
  constexpr strided_view_t(strided_view_t&& x) noexcept : base(std::move(x)) {};

  // From non-const mdspan
  constexpr strided_view_t(mdview_t<value_type, N> x) : base(x) {};
  constexpr strided_view_t(mdstrided_t<value_type, N> x)
      : base(std::move(x)) {};

  // From const mdspan
  constexpr strided_view_t(mdview_t<const value_type, N> x)
    requires(is_const)
      : base(x) {};
  constexpr strided_view_t(mdstrided_t<const value_type, N> x)
    requires(is_const)
      : base(std::move(x)) {};

  // From non-const data holders
  constexpr strided_view_t(data_t<value_type, N>& x) : base(x.view) {};
  constexpr strided_view_t(view_t<value_type, N> x)
      : strided_view_t(x.base_md()) {}

  // Fron const othdata holdersers
  constexpr strided_view_t(const data_t<value_type, N>& x)
    requires(is_const)
      : base(x.view) {}
  constexpr strided_view_t(const view_t<const value_type, N>& x)
    requires(is_const)
      : base(static_cast<mdview_t<T, N>>(x)) {}

  explicit constexpr strided_view_t(T& v) : base(&v, std::array<Index, 1>{1}) {}

  template <typename Self>
  constexpr view_t<T, N> unsafe_view(this Self&& self) {
    assert(self.base::is_exhaustive());
    return mdview_t<T, N>{self.base_md()};
  }

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
    if constexpr (sizeof...(Acc) == N and (integral<Acc> and ...))
      return std::forward<Self>(self).base::operator[](std::forward<Acc>(i)...);
    else
      return left_mdsel_t<T, N, is_const, is_exhaustive, Acc...>(
          left_sub<N>(std::forward<Self>(self), std::forward<Acc>(i)...));
  }

  template <typename Self>
  [[nodiscard]] constexpr auto begin(this Self&& self) {
    return left_mditer(std::forward<Self>(self));
  }

  template <typename Self>
  [[nodiscard]] constexpr auto end(this Self&& self) {
    const auto n = self.extent(0);
    return left_mditer(std::forward<Self>(self)) + n;
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_begin(this Self&& self) {
    return elemwise_mditer(std::forward<Self>(self));
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_end(this Self&& self) {
    const auto n = self.size();
    return elemwise_mditer(std::forward<Self>(self)) + n;
  }

  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) elem_at(this Self&& self, Index i) {
    assert(static_cast<Index>(self.size()) > i);

    if constexpr (N > 1) {
      return std::forward<Self>(self).base::operator[](
          mdpos<N>(self.shape(), i));
    } else {
      return std::forward<Self>(self).base::operator[](i);
    }
  }

  constexpr strided_view_t& operator=(const strided_view_t& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());

    if (this != &r) stdr::copy(elemwise_range(r), elem_begin());

    return *this;
  }

  constexpr strided_view_t& operator=(const value_type& x)
    requires(not is_const)
  {
    std::fill(elem_begin(), elem_end(), x);
    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr strided_view_t& operator=(const R& r)
    requires(not(is_const or std::same_as<strided_view_t, R> or
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
  constexpr strided_view_t& operator=(const U& x)
    requires(not(is_const or std::same_as<strided_view_t, U> or
                 std::same_as<value_type, U> or ranked_convertible_md<U, T, N>))
  {
    const auto sh = shape();

    assert(sh == mdshape(x));

    for (Size i = 0; i < size(); i++)
      elem_at(i) = static_cast<T>(mdvalue(x, sh, i));

    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr strided_view_t& operator+=(const R& r)
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
  constexpr strided_view_t& operator-=(const R& r)
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
  constexpr strided_view_t& operator*=(const R& r)
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
  constexpr strided_view_t& operator/=(const R& r)
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
  constexpr strided_view_t& operator+=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x += v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator-=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x -= v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator*=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x *= v; });
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator/=(U&& r)
    requires(not is_const)
  {
    stdr::for_each(elemwise_range(*this),
                   [v = std::forward_like<T>(r)](T& x) { return x /= v; });
    return *this;
  }

  template <typename Self>
  constexpr auto real(this Self&& self)
    requires(complex_type<T>)
  {
    using U      = complex_subtype_t<T>;
    using cU     = ComplexLayout<U>;
    using ccU    = std::conditional_t<std::is_const_v<Self>, const cU, cU>;
    using elem_t = std::conditional_t<std::is_const_v<Self>, const U, U>;

    std::array<Index, N> strides{};
    for (Size i = 0; i < N; i++) strides[i] = 2 * self.stride(i);

    return strided_view_t<elem_t, N>{
        mdstrided_t<U, N>(&(reinterpret_cast<ccU*>(self.data_handle())->real),
                          {self.shape(), strides})};
  }

  template <typename Self>
  constexpr auto imag(this Self&& self)
    requires(complex_type<T>)
  {
    using U      = complex_subtype_t<T>;
    using cU     = ComplexLayout<U>;
    using ccU    = std::conditional_t<std::is_const_v<Self>, const cU, cU>;
    using elem_t = std::conditional_t<std::is_const_v<Self>, const U, U>;

    std::array<Index, N> strides{};
    for (Size i = 0; i < N; i++) strides[i] = 2 * self.stride(i);

    return strided_view_t<elem_t, N>{
        mdstrided_t<U, N>(&(reinterpret_cast<ccU*>(self.data_handle())->imag),
                          {self.shape(), strides})};
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
};

static_assert(
    std::random_access_iterator<left_mditer<strided_view_t<Numeric, 1>>>);
static_assert(
    std::random_access_iterator<left_mditer<strided_view_t<Numeric, 2>>>);
static_assert(
    std::random_access_iterator<left_mditer<strided_view_t<Numeric, 9>>>);
static_assert(
    std::random_access_iterator<elemwise_mditer<strided_view_t<Numeric, 1>>>);
static_assert(
    std::random_access_iterator<elemwise_mditer<strided_view_t<Numeric, 2>>>);
static_assert(
    std::random_access_iterator<elemwise_mditer<strided_view_t<Numeric, 9>>>);

static_assert(strided_view_t<Numeric, 10>::is_strided);
static_assert(not strided_view_t<Numeric, 10>::is_const);
static_assert(strided_view_t<const Numeric, 10>::is_const);
static_assert(not strided_view_t<Numeric, 10>::is_exhaustive);
static_assert(strided_view_t<Numeric, 10>::is_unique);
}  // namespace matpack

template <typename T, Size N>
struct std::formatter<matpack::strided_view_t<T, N>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::strided_view_t<T, N>& v,
                              FmtContext& ctx) const {
    using std::ranges::views::take, std::ranges::views::drop;

    const Index n              = v.shape()[0];
    const std::string_view sep = inner_fmt().tags.sep(N != 1);

    tags.add_if_bracket(ctx, '[');
    if constexpr (N > 1) tags.add_if_bracket(ctx, '\n');

    if (n > 0) {
      tags.format(ctx, v[0]);

      if (tags.short_str and n > 8) {
        for (auto&& a : v | take(3) | drop(1)) tags.format(ctx, sep, a);
        tags.format(ctx, sep, "..."sv);
        for (auto&& a : v | drop(n - 3)) tags.format(ctx, sep, a);
      } else {
        for (auto&& a : v | drop(1)) tags.format(ctx, sep, a);
      }
    }

    if constexpr (N > 1) tags.add_if_bracket(ctx, '\n');
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

using StridedVectorView  = matpack::strided_view_t<Numeric, 1>;
using StridedMatrixView  = matpack::strided_view_t<Numeric, 2>;
using StridedTensor3View = matpack::strided_view_t<Numeric, 3>;
using StridedTensor4View = matpack::strided_view_t<Numeric, 4>;
using StridedTensor5View = matpack::strided_view_t<Numeric, 5>;
using StridedTensor6View = matpack::strided_view_t<Numeric, 6>;
using StridedTensor7View = matpack::strided_view_t<Numeric, 7>;

using StridedComplexVectorView  = matpack::strided_view_t<Complex, 1>;
using StridedComplexMatrixView  = matpack::strided_view_t<Complex, 2>;
using StridedComplexTensor3View = matpack::strided_view_t<Complex, 3>;
using StridedComplexTensor4View = matpack::strided_view_t<Complex, 4>;
using StridedComplexTensor5View = matpack::strided_view_t<Complex, 5>;
using StridedComplexTensor6View = matpack::strided_view_t<Complex, 6>;
using StridedComplexTensor7View = matpack::strided_view_t<Complex, 7>;

using StridedConstVectorView  = matpack::strided_view_t<const Numeric, 1>;
using StridedConstMatrixView  = matpack::strided_view_t<const Numeric, 2>;
using StridedConstTensor3View = matpack::strided_view_t<const Numeric, 3>;
using StridedConstTensor4View = matpack::strided_view_t<const Numeric, 4>;
using StridedConstTensor5View = matpack::strided_view_t<const Numeric, 5>;
using StridedConstTensor6View = matpack::strided_view_t<const Numeric, 6>;
using StridedConstTensor7View = matpack::strided_view_t<const Numeric, 7>;

using StridedConstComplexVectorView = matpack::strided_view_t<const Complex, 1>;
using StridedConstComplexMatrixView = matpack::strided_view_t<const Complex, 2>;
using StridedConstComplexTensor3View =
    matpack::strided_view_t<const Complex, 3>;
using StridedConstComplexTensor4View =
    matpack::strided_view_t<const Complex, 4>;
using StridedConstComplexTensor5View =
    matpack::strided_view_t<const Complex, 5>;
using StridedConstComplexTensor6View =
    matpack::strided_view_t<const Complex, 6>;
using StridedConstComplexTensor7View =
    matpack::strided_view_t<const Complex, 7>;
