#pragma once

#include <concepts>
#include <functional>
#include <type_traits>

#include "matpack_mdspan_common.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_elemwise_mditer.h"

namespace matpack {
template <typename T, Size... dims>
struct cdata_t {
  constexpr static bool magic_cdata_t_v = true;
  constexpr static Size N               = sizeof...(dims);
  constexpr static Size ndata           = (dims * ...);

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
  static_assert(N > 0);

  std::array<T, ndata> data;

  static constexpr std::array<Index, N> shape_ = {dims...};
  static constexpr const std::array<Index, N>& shape() { return shape_; }

  view_t<const T, N> view() const {
    return mdview_t<T, N>{const_cast<T*>(data.data()), shape_};
  }

  view_t<T, N> view() {
    return mdview_t<T, N>{const_cast<T*>(data.data()), shape_};
  }

  auto base_md() const { return view().base_md(); }
  auto base_md() { return view().base_md(); }

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

  template <typename Self, std::integral... inds>
  constexpr auto view_as(this Self&& self, inds... exts) {
    return std::forward<Self>(self).view().view_as(std::forward<inds>(exts)...);
  }

  template <access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](Acc&&... i) {
    return view()[std::forward<Acc>(i)...];
  }

  template <access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](Acc&&... i) const {
    return view()[std::forward<Acc>(i)...];
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
  constexpr auto stride(this Self&& self, std::integral auto i) {
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
