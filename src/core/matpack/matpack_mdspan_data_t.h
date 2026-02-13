#pragma once

#include <configtypes.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include "matpack_mdspan_common_sizes.h"
#include "matpack_mdspan_elemwise_mditer.h"
#include "matpack_mdspan_view_t.h"

namespace matpack {
namespace {
template <typename T, Size N, Size I = 0, typename arg, typename... args>
constexpr std::pair<std::array<Index, N>, T> to_pair_i(
    std::array<Index, N>& ret, arg&& x, args&&... y) {
  if constexpr (sizeof...(args) > 0) {
    ret[I] = static_cast<Index>(x);
    return to_pair_i<T, N, I + 1>(ret, std::forward<args>(y)...);
  } else {
    return {ret, std::forward_like<T>(x)};
  }
}

template <typename T, Size N, typename... args>
constexpr std::pair<std::array<Index, N>, T> to_pair(args&&... y) {
  std::array<Index, N> ret{};
  return to_pair_i<T, N, 0>(ret, std::forward<args>(y)...);
}
}  // namespace

template <typename T, Size N>
class [[nodiscard]] data_t {
  static_assert(not std::is_const_v<T>);
  static_assert(N > 0);

  std::vector<T> data;
  view_t<T, N> view;

  void reset_view(const std::array<Index, N>& sz) noexcept {
    view.base_set({data.data(), sz});
  }

  template <typename U, Size M>
  friend class data_t;
  template <typename U, Size M>
  friend struct view_t;
  template <typename U, Size M>
  friend struct strided_view_t;

 public:
  using base = view_t<T, N>::base;

  constexpr auto base_md() const { return view.base_md(); }

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

  explicit constexpr data_t(std::pair<std::array<Index, N>, T> x)
      : data(mdsize(x.first), x.second), view(base{data.data(), x.first}) {}

  constexpr data_t(std::array<Index, N> sz, T x = {})
      : data_t(std::pair<std::array<Index, N>, T>{sz, x}) {}

  constexpr data_t() : data_t(std::array<Index, N>{}) {}
  constexpr data_t(const data_t& x)
      : data(x.data), view(base{data.data(), x.shape()}) {}
  constexpr data_t(data_t&& x) noexcept
      : data(std::move(x.data)), view(base{data.data(), x.shape()}) {
    x.clear();
  }

  template <integral... inds>
  explicit constexpr data_t(inds... ind)
    requires(sizeof...(inds) == N)
      : data_t(std::array{static_cast<Index>(ind)...}) {}

  template <class... args>
  constexpr data_t(args... arg)
    requires(sizeof...(args) == (N + 1))
      : data_t(to_pair<T, N>(std::forward<args>(arg)...)) {}

  explicit constexpr data_t(std::vector<T>&& r)
    requires(N == 1)
      : data(std::move(r)), view(data) {}

  data_t(std::initializer_list<T> r)
    requires(N == 1)
      : data_t(std::vector<T>(stdr::begin(r), stdr::end(r))) {}

  data_t(std::from_range_t, stdr::input_range auto&& r)
    requires(N == 1)
  {
    std::vector<T> v;
    v.reserve(stdr::size(r));
    for (auto&& x : r) v.emplace_back(static_cast<T>(x));
    *this = std::move(v);
  }

  template <mdvalue_type_compatible<T> U>
  explicit constexpr data_t(const U& x) : data_t(mdshape(x)) {
    view = x;
  }

  template <typename U>
  explicit constexpr data_t(const view_t<U, N>& x) : data_t(x.shape()) {
    view = x;
  }

  template <typename U>
  explicit constexpr data_t(const strided_view_t<U, N>& x) : data_t(x.shape()) {
    view = x;
  }

  template <ranked_convertible_md<T, N> Real, ranked_convertible_md<T, N> Imag>
  constexpr data_t(const Real& real, const Imag& imag)
    requires(complex_type<T> and
             std::same_as<typename Real::value_type, complex_subtype_t<T>> and
             std::same_as<typename Imag::value_type, complex_subtype_t<T>>)
      : data_t(real.shape()) {
    assert(real.shape() == imag.shape());
    auto r = real | by_elem;
    auto i = imag | by_elem;
    std::transform(stdr::begin(r),
                   stdr::end(r),
                   stdr::begin(i),
                   data.begin(),
                   [](auto a, auto b) { return T{a, b}; });
  }

  constexpr operator view_t<T, N>&() { return view; }
  constexpr operator const view_t<T, N>&() const { return view; }
  constexpr operator std::span<T>()
    requires(N == 1)
  {
    return {begin(), end()};
  }
  constexpr operator std::span<const T>() const
    requires(N == 1)
  {
    return {begin(), end()};
  }

  /*

   Common operators using both data and view

   */

  constexpr data_t& resize(std::array<Index, N> sz) {
    data.resize(mdsize<N>(sz));
    reset_view(sz);
    return *this;
  }

  template <integral... inds>
  constexpr data_t& resize(inds... ind)
    requires(sizeof...(inds) == N)
  {
    return resize({static_cast<Index>(ind)...});
  }

  constexpr data_t& operator=(data_t&& x) noexcept {
    data = std::move(x.data);
    reset_view(x.view.shape());
    x.clear();
    return *this;
  }

  constexpr data_t& operator=(const data_t& x) noexcept {
    data = x.data;
    reset_view(x.view.shape());
    return *this;
  }

  constexpr data_t& operator=(const value_type& x) {
    view = x;
    return *this;
  }

  template <ranked_convertible_md<T, N> R>
  constexpr data_t& operator=(const R& r)
    requires(not(std::same_as<data_t, R> or std::same_as<value_type, R>))
  {
    resize(r.shape());
    view = r;
    return *this;
  }

  template <mdvalue_type_compatible<T> U>
  constexpr data_t& operator=(const U& x)
    requires(not(std::same_as<data_t, U> or std::same_as<value_type, U> or
                 ranked_convertible_md<U, T, N>))
  {
    *this = data_t(x);
    return *this;
  }

  template <Size M>
  constexpr data_t<T, M> reshape(const std::array<Index, M>& exts) && {
    data_t<T, M> out{};
    out.data = std::move(data);
    out.reset_view(exts);
    reset_view({});
    return out;
  }

  template <integral... inds>
  constexpr auto reshape(inds... exts) && {
    return std::move(*this).template reshape<sizeof...(inds)>(
        std::array{static_cast<Index>(exts)...});
  }

  constexpr auto flatten() && {
    const auto n = size();
    return std::move(*this).reshape(n);
  }

  constexpr void swap(data_t& other) noexcept {
    data.swap(other.data);
    view.swap(other.view);
  }

  /*

  Operations purely to mimic std::vector

   */

  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) front(this Self&& self) {
    assert(not self.empty());
    return std::forward<Self>(self)[0];
  }

  template <typename Self>
  [[nodiscard]] constexpr decltype(auto) back(this Self&& self) {
    assert(not self.empty());
    const Size n = self.size();
    return std::forward<Self>(self)[n - 1];
  }

  template <std::convertible_to<T> U>
  constexpr void push_back(U&& x)
    requires(N == 1)
  {
    data.push_back(static_cast<T>(std::forward<U>(x)));
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  constexpr void pop_back()
    requires(N == 1)
  {
    assert(not data.empty());
    data.pop_back();
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  constexpr void erase(std::vector<T>::iterator pos)
    requires(N == 1)
  {
    assert(not data.empty());
    data.erase(pos);
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  constexpr void erase(std::vector<T>::const_iterator pos)
    requires(N == 1)
  {
    assert(not data.empty());
    data.erase(pos);
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  constexpr void erase(std::vector<T>::iterator first,
                       std::vector<T>::iterator last)
    requires(N == 1)
  {
    assert(not data.empty());
    data.erase(first, last);
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  constexpr void erase(std::vector<T>::const_iterator first,
                       std::vector<T>::const_iterator last)
    requires(N == 1)
  {
    assert(not data.empty());
    data.erase(first, last);
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
  }

  void reserve(Size n) { data.reserve(n); }

  template <class... Args>
  constexpr T& emplace_back(Args&&... args)
    requires(N == 1)
  {
    auto& out = data.emplace_back(std::forward<Args>(args)...);
    reset_view(std::array<Index, N>{static_cast<Index>(data.size())});
    return out;
  }

  void clear() {
    data.clear();
    reset_view({});
  }

  /*

  Operations purely to mimic std::mdspan

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

  template <typename Self, Size M>
  constexpr auto view_as(this Self&& self, const std::array<Index, M>& exts) {
    return std::forward<Self>(self).view.view_as(exts);
  }

  template <typename Self, integral... inds>
  constexpr auto view_as(this Self&& self, inds... exts) {
    return std::forward<Self>(self).view.view_as(std::forward<inds>(exts)...);
  }

  template <typename Self, access_operator... Acc>
  [[nodiscard]] constexpr decltype(auto) operator[](this Self&& self,
                                                    Acc&&... i)
    requires(sizeof...(Acc) <= N)
  {
    return std::forward<Self>(self).view[std::forward<Acc>(i)...];
  }

  template <typename U>
  constexpr data_t& operator+=(U&& x) {
    view += std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr data_t& operator-=(U&& x) {
    view -= std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr data_t& operator*=(U&& x) {
    view *= std::forward<U>(x);
    return *this;
  }

  template <typename U>
  constexpr data_t& operator/=(U&& x) {
    view /= std::forward<U>(x);
    return *this;
  }

  template <typename Self>
  constexpr decltype(auto) elem_at(this Self&& self, Index i) {
    return std::forward<Self>(self).data[i];
  }

  template <typename Self>
  constexpr auto shape(this Self&& self) {
    return std::forward<Self>(self).view.shape();
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
    return std::forward<Self>(self).view.begin();
  }

  template <typename Self>
  constexpr auto end(this Self&& self) {
    return std::forward<Self>(self).view.end();
  }

  template <typename Self>
  constexpr auto imag(this Self&& self) {
    return std::forward<Self>(self).view.imag();
  }

  template <typename Self>
  constexpr auto real(this Self&& self) {
    return std::forward<Self>(self).view.real();
  }

  template <typename Self>
  constexpr auto accessor(this Self&& self) {
    return std::forward<Self>(self).view.accessor();
  }

  template <typename Self>
  constexpr auto data_handle(this Self&& self) {
    return std::forward<Self>(self).data.data();
  }

  template <typename Self>
  constexpr auto empty(this Self&& self) {
    return std::forward<Self>(self).view.empty();
  }

  template <typename Self>
  constexpr auto extent(this Self&& self, Index i) {
    return std::forward<Self>(self).view.extent(i);
  }

  template <typename Self>
  constexpr auto extents(this Self&& self) {
    return std::forward<Self>(self).view.extents();
  }

  template <typename Self>
  constexpr auto mapping(this Self&& self) {
    return std::forward<Self>(self).view.mapping();
  }

  template <typename Self>
  constexpr auto size(this Self&& self) {
    return std::forward<Self>(self).view.size();
  }

  template <typename Self>
  constexpr auto stride(this Self&& self, integral auto i) {
    return std::forward<Self>(self).view.stride(i);
  }
};
}  // namespace matpack

template <typename T, Size N>
struct std::formatter<matpack::data_t<T, N>> {
  std::formatter<matpack::view_t<const T, N>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::data_t<T, N>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
  }
};

using Vector  = matpack::data_t<Numeric, 1>;
using Matrix  = matpack::data_t<Numeric, 2>;
using Tensor3 = matpack::data_t<Numeric, 3>;
using Tensor4 = matpack::data_t<Numeric, 4>;
using Tensor5 = matpack::data_t<Numeric, 5>;
using Tensor6 = matpack::data_t<Numeric, 6>;
using Tensor7 = matpack::data_t<Numeric, 7>;

using ComplexVector  = matpack::data_t<Complex, 1>;
using ComplexMatrix  = matpack::data_t<Complex, 2>;
using ComplexTensor3 = matpack::data_t<Complex, 3>;
using ComplexTensor4 = matpack::data_t<Complex, 4>;
using ComplexTensor5 = matpack::data_t<Complex, 5>;
using ComplexTensor6 = matpack::data_t<Complex, 6>;
using ComplexTensor7 = matpack::data_t<Complex, 7>;

using IndexVector = matpack::data_t<Index, 1>;
