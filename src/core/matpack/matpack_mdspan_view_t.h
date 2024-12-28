#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iterator>
#include <ranges>
#include <type_traits>

#include "matpack_mdspan_elemwise_mditer.h"
#include "matpack_mdspan_strided_view_t.h"

namespace matpack {
//! A wrapper for mdspan that's always C-like in layout
template <class T, Size N>
struct view_t final : public mdview_t<T, N> {
  using base = mdview_t<T, N>;

  mdview_t<T, N> base_md() const { return *this; }

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

  void base_set(base v) noexcept { static_cast<base&>(*this) = v; }

  static constexpr bool is_const      = std::is_const_v<T>;
  static constexpr bool is_exhaustive = is_always_exhaustive();
  static constexpr bool is_strided    = is_always_strided();
  static constexpr bool is_unique     = is_always_unique();

  constexpr view_t() : base() {};
  constexpr view_t(const view_t& x) : base(x) {};
  constexpr view_t(view_t&& x) noexcept : base(std::move(x)) {};

  constexpr view_t(mdview_t<value_type, N> x) : base(x) {};
  constexpr view_t(mdview_t<const value_type, N> x)
    requires(is_const)
      : base(x) {};

  constexpr view_t(const data_t<value_type, N>& x)
    requires(is_const)
      : base(x.view) {}

  explicit constexpr view_t(T& v) : base(&v, std::array<Index, 1>{1}) {}

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
  [[nodiscard]] constexpr decltype(auto) operator[](this Self&& self, Acc&&... i) {
    if constexpr (sizeof...(Acc) == N and
                  (std::integral<std::remove_cvref_t<Acc>> and ...))
      return std::forward<Self>(self).base::operator[](std::forward<Acc>(i)...);
    else
      return left_mdsel_t<T, N, is_const, is_exhaustive, Acc...>{
          left_sub<N>(std::forward<Self>(self), std::forward<Acc>(i)...)};
  }

  template <typename Self>
  [[nodiscard]] constexpr auto begin(this Self&& self) {
    if constexpr (N == 1) {
      return std::forward<Self>(self).elem_begin();
    } else {
      return left_mditer(std::forward<Self>(self));
    }
  }

  template <typename Self>
  [[nodiscard]] constexpr auto end(this Self&& self) {
    if constexpr (N == 1) {
      return std::forward<Self>(self).elem_end();
    } else {
      const auto n = self.extent(0);
      return left_mditer(std::forward<Self>(self)) + n;
    }
  }

  template <typename Self, Size M>
  constexpr view_t<
      std::conditional_t<std::is_const_v<Self>, std::add_const_t<T>, T>,
      M>
  view_as(this Self&& self, const std::array<Index, M>& exts) {
    assert(self.size() == mdsize(exts));
    return mdview_t<T, M>{self.data_handle(), exts};
  }

  template <typename Self, std::integral... NewExtents>
  constexpr auto view_as(this Self&& self, NewExtents... exts) {
    return std::forward<Self>(self).view_as(
        std::array{static_cast<Index>(exts)...});
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_begin(this Self&& self) {
    return std::forward<Self>(self).data_handle();
  }

  template <typename Self>
  [[nodiscard]] constexpr auto elem_end(this Self&& self) {
    const auto n = self.size();
    return std::forward<Self>(self).data_handle() + n;
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
    stdr::copy(elemwise_range(r), elem_begin());
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
    unary_transform(r, [](const auto& v) { return static_cast<T>(v); });
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

  /** Reduce self by binary reduce operation.
   * 
   * Example use: see matpack_mdspan_helpers_reduce.h
   * 
   * Possible usecase: summing all elements in a object.
   * 
   * By default, this will sum all elements in the object.
   * 
   * @param self The self object
   * @param init The initial value, defaults to 0 for most common types
   * @param op The reduction operation, defaults to std::plus<>{}
   * @return As specified by std::reduce
   */
  template <typename Self, class U = T, class BinaryOp = std::plus<>>
  constexpr auto reduce(this Self&& self, U&& init = {}, BinaryOp&& op = {}) {
    return std::reduce(self.elem_begin(),
                       self.elem_end(),
                       std::forward<U>(init),
                       std::forward<BinaryOp>(op));
  }

  /** Reduce self by unary transform followed by binary reduce.
   * 
   * Example use: see matpack_mdspan_helpers_reduce.h
   * 
   * Possible usecase: sum ignoring all NaN by using unary transform to set NaN to 0, then binary reduce by sum.
   * 
   * Note that by default this function acts exactly like the reduce() member function, but the
   * latter does not invoke a transformation.  So unless you change the default unary transform
   * operation, you should use reduce() instead.
   * 
   * @param self The self object
   * @param init The initial value, defaults to 0 for most common types
   * @param reduce The binary reduce operation, defaults to std::plus<>{}
   * @param transform The unary transform operation, defaults to std::identity{}
   * @return As specified by std::transform_reduce
   */
  template <typename Self,
            class U        = T,
            class BinaryOp = std::plus<>,
            class UnaryOp  = std::identity>
  constexpr auto unary_transform_reduce(this Self&& self,
                                        U&& init            = {},
                                        BinaryOp&& reduce   = {},
                                        UnaryOp&& transform = {}) {
    return std::transform_reduce(self.elem_begin(),
                                 self.elem_end(),
                                 std::forward<U>(init),
                                 std::forward<BinaryOp>(reduce),
                                 std::forward<UnaryOp>(transform));
  }

  /** Reduce self and another by binary transform followed by binary reduce.
   * 
   * Example use: see matpack_mdspan_helpers_reduce.h
   * 
   * Possible usecase: the sum of the mutliplication of two objects.  I.e., the dot-product
   * 
   * By default, this will sum the elementwise product of the two objects.
   * 
   * @param self The self object
   * @param other The other object
   * @param init The initial value, defaults to 0 for most common types
   * @param reduce The binary reduce operation, defaults to std::plus<>{}
   * @param transform The binary transform operation, defaults to std::multiplies{}
   * @return As specified by std::transform_reduce
   */
  template <typename Self,
            elemwise_mditerable Other,
            class U         = T,
            class BinaryOp1 = std::plus<>,
            class BinaryOp2 = std::multiplies<>>
  constexpr auto binary_transform_reduce(this Self&& self,
                                         Other&& other,
                                         U&& init              = {},
                                         BinaryOp1&& reduce    = {},
                                         BinaryOp2&& transform = {}) {
    auto r = elemwise_range(std::forward<Other>(other));
    assert(stdr::size(r) == self.size());
    return std::transform_reduce(stdr::begin(r),
                                 stdr::end(r),
                                 std::forward<Self>(self).elem_begin(),
                                 std::forward<U>(init),
                                 std::forward<BinaryOp1>(reduce),
                                 std::forward<BinaryOp2>(transform));
  }

  /** Transform self by unary operation on other range.
   * 
   * Example use: See operator overloads of this class.
   * 
   * Possible usecase: Add a constant with operator+= to self, the other range is then also self.
   * 
   * See stdr::transform for more information on projections.
   * 
   * @param self The self object
   * @param other The other range
   * @param unary_op The operation to apply to each element
   * @param proj A projection to apply to each element before applying the unary operation, defaults to no projection
   * @return As defined by stdr::transform
   */
  template <typename Self,
            elemwise_mditerable Other,
            std::copy_constructible F,
            class Proj = std::identity>
  constexpr auto unary_transform(this Self&& self,
                                 Other&& other,
                                 F&& unary_op,
                                 Proj&& proj = {}) {
    auto r = elemwise_range(std::forward<Other>(other));
    assert(stdr::size(r) == self.size());
    return stdr::transform(r,
                           std::forward<Self>(self).elem_begin(),
                           std::forward<F>(unary_op),
                           std::forward<Proj>(proj));
  }

  /** Transform self by binary operation on other ranges.
   * 
   * Example use: See operator overloads of this class.
   * 
   * Possible usecase: Add elementwise with operator+= to self, one of the other ranges is then also self.
   * 
   * See stdr::transform for more information on projections.
   * 
   * @param self The self object
   * @param other1 The other range
   * @param other2 The other range
   * @param binary_op The operation to apply to the pair of each element
   * @param proj1 A projection to apply to each element of other1 before the binary operation, defaults to no projection
   * @param proj2 A projection to apply to each element of other2 before the binary operation, defaults to no projection
   * @return As defined by stdr::transform
   */
  template <typename Self,
            elemwise_mditerable Other1,
            elemwise_mditerable Other2,
            std::copy_constructible F,
            class Proj1 = std::identity,
            class Proj2 = std::identity>
  constexpr auto binary_transform(this Self&& self,
                                  Other1&& other1,
                                  Other2&& other2,
                                  F&& binary_op,
                                  Proj1&& proj1 = {},
                                  Proj2&& proj2 = {}) {
    auto r1 = elemwise_range(std::forward<Other1>(other1));
    auto r2 = elemwise_range(std::forward<Other2>(other2));
    assert(stdr::size(r1) == self.size() and stdr::size(r2) == self.size());
    return stdr::transform(r1,
                           r2,
                           std::forward<Self>(self).elem_begin(),
                           std::forward<F>(binary_op),
                           std::forward<Proj1>(proj1),
                           std::forward<Proj2>(proj2));
  }

  template <elemwise_mditerable R>
  constexpr view_t& operator+=(const R& r1) {
    binary_transform(*this, r1, std::plus<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr view_t& operator-=(const R& r1) {
    binary_transform(*this, r1, std::minus<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr view_t& operator*=(const R& r1) {
    binary_transform(*this, r1, std::multiplies<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr view_t& operator/=(const R& r1) {
    binary_transform(*this, r1, std::divides<>{});
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator+=(U&& r) {
    const auto add = [v = std::forward_like<T>(r)](auto&& x) { return x + v; };
    unary_transform(*this, add);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator-=(U&& r) {
    const auto red = [v = std::forward_like<T>(r)](auto&& x) { return x - v; };
    unary_transform(*this, red);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator/=(U&& r) {
    const auto div = [v = std::forward_like<T>(r)](auto&& x) { return x / v; };
    unary_transform(*this, div);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr view_t& operator*=(U&& r) {
    const auto mul = [v = std::forward_like<T>(r)](auto&& x) { return x * v; };
    unary_transform(*this, mul);
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
};

static_assert(std::random_access_iterator<left_mditer<view_t<Numeric, 1>>>);
static_assert(std::random_access_iterator<left_mditer<view_t<Numeric, 2>>>);
static_assert(std::random_access_iterator<left_mditer<view_t<Numeric, 9>>>);
static_assert(std::random_access_iterator<elemwise_mditer<view_t<Numeric, 1>>>);
static_assert(std::random_access_iterator<elemwise_mditer<view_t<Numeric, 2>>>);
static_assert(std::random_access_iterator<elemwise_mditer<view_t<Numeric, 9>>>);

static_assert(view_t<Numeric, 10>::is_strided);
static_assert(not view_t<Numeric, 10>::is_const);
static_assert(view_t<const Numeric, 10>::is_const);
static_assert(view_t<Numeric, 10>::is_exhaustive);
static_assert(view_t<Numeric, 10>::is_unique);
}  // namespace matpack

template <typename T, Size N>
struct std::formatter<matpack::view_t<T, N>> {
  std::formatter<matpack::strided_view_t<const T, N>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::view_t<T, N>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
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
