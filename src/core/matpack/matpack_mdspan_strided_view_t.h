#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iterator>
#include <type_traits>

#include "matpack_mdspan_common.h"
#include "matpack_mdspan_mditer.h"

namespace matpack {
template <class T, Size N>
struct strided_view_t final : public mdstrided_t<T, N> {
  using base = mdstrided_t<T, N>;

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

  constexpr strided_view_t()                          = default;
  constexpr strided_view_t(const strided_view_t&)     = default;
  constexpr strided_view_t(strided_view_t&&) noexcept = default;

  constexpr strided_view_t(const mdstrided_t<T, N>& x) : base(x) {};
  constexpr strided_view_t(mdstrided_t<T, N>&& x) : base(std::move(x)) {};
  constexpr strided_view_t(const mdview_t<T, N>& x) : base(x) {};
  constexpr strided_view_t(mdview_t<T, N>&& x) : base(std::move(x)) {};

  explicit constexpr strided_view_t(base::data_handle_type p,
                                    const std::array<Index, N>& ext)
      : base(mdview_t<T, N>(p, ext)) {}
  explicit constexpr strided_view_t(base::data_handle_type p,
                                    const base::extents_type& ext)
      : base(p, ext) {}
  explicit constexpr strided_view_t(base::data_handle_type p,
                                    const base::mapping_type& m)
      : base(p, m) {}
  explicit constexpr strided_view_t(base::data_handle_type p,
                                    const base::mapping_type& m,
                                    const base::accessor_type& a)
      : base(p, m, a) {}

  constexpr strided_view_t(const strided_view_t<value_type, N>& x)
    requires(is_const)
      : base(x) {}

  constexpr strided_view_t(data_t<value_type, N>& x) : base(x.view) {}
  constexpr strided_view_t(const data_t<value_type, N>& x) : base(x.view) {}

  constexpr strided_view_t<value_type, N>(view_t<value_type, N>& x)
    requires(not is_const)
      : base(x) {}
  constexpr strided_view_t<const value_type, N>(view_t<value_type, N>& x)
    requires(is_const)
      : base(x) {}
  constexpr strided_view_t<const value_type, N>(const view_t<value_type, N>& x)
    requires(is_const)
      : base(x) {}
  constexpr strided_view_t<const value_type, N>(
      const view_t<const value_type, N>& x)
    requires(is_const)
      : base(x) {}

  explicit constexpr strided_view_t(T& v) {
    std::array<Index, N> exts{};
    exts.fill(1);
    base_set(base(&v, exts));
  }

  explicit constexpr strided_view_t(std::vector<value_type>& v)
    requires(N == 1 and not is_const)
      : base(const_cast<value_type*>(v.data()), v.size()) {}

  explicit constexpr strided_view_t(const std::vector<value_type>& v)
    requires(N == 1 and is_const)
      : base(const_cast<value_type*>(v.data()), v.size()) {}

  template <Size M>
  explicit constexpr strided_view_t(const strided_view_t<T, M>& x)
    requires(M < N)
  {
    std::array<Index, N> exts{}, strides{};
    exts.fill(1);
    strides.fill(1);
    for (Size i = 0; i < M; i++) {
      exts[i]    = x.extent(i);
      strides[i] = x.stride(i);
    }
    base_set(base(x.data_handle(), exts, strides));
  }

  template <typename Self>
  constexpr view_t<T, N> unsafe_view(this Self&& self) {
    assert(self.base::is_exhaustive());
    return view_t<T, N>{self.base::data_handle(), self.base::extents()};
  }

  [[nodiscard]] constexpr auto shape() const {
    std::array<Index, N> out;
    for (Size i = 0; i < N; i++) out[i] = this->extent(i);
    return out;
  }

  template <typename Self, access_operator... Acc>
  [[nodiscard]] constexpr auto operator[](this Self&& self, Acc&&... i)
      -> left_mdsel_t<T, N, is_const, is_exhaustive, Acc...> {
    if constexpr (sizeof...(Acc) == N and
                  (std::integral<std::remove_cvref_t<Acc>> and ...))
      return std::forward<Self>(self).base::operator[](std::forward<Acc>(i)...);
    else
      return left_mdsel_t<T, N, is_const, is_exhaustive, Acc...>{
          left_sub<N>(std::forward<Self>(self), std::forward<Acc>(i)...)};
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
      auto sh = self.shape();
      return std::forward<Self>(self).base::operator[](mdpos<N>(sh, i));
    } else {
      return std::forward<Self>(self).base::operator[](i);
    }
  }

  constexpr strided_view_t& operator=(const strided_view_t& r)
    requires(not is_const)
  {
    assert(shape() == r.shape());
    stdr::copy(elemwise_range(r), elem_begin());
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
    unary_transform(r, [](const auto& v) { return static_cast<T>(v); });
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
  constexpr decltype(auto) reduce(this Self&& self,
                                  U&& init      = {},
                                  BinaryOp&& op = std::plus<>{}) {
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
                                        BinaryOp&& reduce   = std::plus<>{},
                                        UnaryOp&& transform = std::identity{}) {
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
  constexpr decltype(auto) binary_transform_reduce(
      this Self&& self,
      Other&& other,
      U&& init              = {},
      BinaryOp1&& reduce    = std::plus<>{},
      BinaryOp2&& transform = std::multiplies<>{}) {
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
  constexpr decltype(auto) unary_transform(this Self&& self,
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
  constexpr decltype(auto) binary_transform(this Self&& self,
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
  constexpr strided_view_t& operator+=(const R& r1) {
    binary_transform(*this, r1, std::plus<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr strided_view_t& operator-=(const R& r1) {
    binary_transform(*this, r1, std::minus<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr strided_view_t& operator*=(const R& r1) {
    binary_transform(*this, r1, std::multiplies<>{});
    return *this;
  }

  template <elemwise_mditerable R>
  constexpr strided_view_t& operator/=(const R& r1) {
    binary_transform(*this, r1, std::divides<>{});
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator+=(U&& r) {
    const auto add = [v = std::forward_like<T>(r)](auto&& x) { return x + v; };
    unary_transform(*this, add);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator-=(U&& r) {
    const auto red = [v = std::forward_like<T>(r)](auto&& x) { return x - v; };
    unary_transform(*this, red);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator/=(U&& r) {
    const auto div = [v = std::forward_like<T>(r)](auto&& x) { return x / v; };
    unary_transform(*this, div);
    return *this;
  }

  template <std::convertible_to<T> U>
  constexpr strided_view_t& operator*=(U&& r) {
    const auto mul = [v = std::forward_like<T>(r)](auto&& x) { return x * v; };
    unary_transform(*this, mul);
    return *this;
  }

  template <typename Self>
  constexpr decltype(auto) real(this Self&& self)
    requires(complex_type<T>)
  {
    using U   = std::remove_cvref_t<decltype(T{}.real())>;
    using cU  = ComplexLayout<U>;
    using ccU = std::conditional_t<std::is_const_v<Self>, const cU, cU>;

    std::array<Index, N> strides{};
    for (Size i = 0; i < N; i++) strides[i] = 2 * self.stride(i);

    return strided_view_t<U, N>(
        &(reinterpret_cast<ccU*>(self.data_handle())->real),
        {self.shape(), strides});
  }

  template <typename Self>
  constexpr decltype(auto) imag(this Self&& self)
    requires(complex_type<T>)
  {
    using U   = std::remove_cvref_t<decltype(T{}.real())>;
    using cU  = ComplexLayout<U>;
    using ccU = std::conditional_t<std::is_const_v<Self>, const cU, cU>;

    std::array<Index, N> strides{};
    for (Size i = 0; i < N; i++) strides[i] = 2 * self.stride(i);

    return strided_view_t<U, N>(
        &(reinterpret_cast<ccU*>(self.data_handle())->imag),
        {self.shape(), strides});
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
