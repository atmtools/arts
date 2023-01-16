#pragma once

#include "matpack_complex.h"
#include "matpack_concepts.h"
#include "matpack_view.h"

#include <exception>
#include <tuple>

namespace matpack {
namespace detail {
template <Index I, typename T, Index N, typename value, typename... rest>
void front_array_impl(std::array<T, N> &out, value &&v, rest &&...x) {
  if constexpr (I < N) {
    std::get<I>(out) = static_cast<T>(std::forward<value>(v));
    front_array_impl<I + 1, T, N>(out, std::forward<rest>(x)...);
  }
}
} // namespace detail

template <typename T, Index N, typename... input,
          Index M = sizeof...(input)>
constexpr std::array<T, N> front_array(input... in)
  requires(N <= M and N > 0)
{
  using namespace detail;
  std::array<T, N> out;
  front_array_impl<0, T, N>(out, in...);
  return out;
}

template<typename... rest, Index N=sizeof...(rest)>
auto get_last(rest&&... r) requires(N>0) {
  return std::get<N-1>(std::tuple{std::forward<rest>(r)...});
}

template <typename T, Index N>
class matpack_data {
  std::vector<T> data;

  using view_type = matpack_view<T, N, false, false>;
  view_type view;

  using me_view = typename view_type::me_view;
  using ms_view = typename view_type::ms_view;
  using ce_view = typename view_type::ce_view;
  using cs_view = typename view_type::cs_view;

  template <typename U, Index M> friend class matpack_data;
  template <typename U, Index M, bool c, bool s> friend class matpack_view;

public:
  constexpr matpack_data(const std::array<Index, N> &sz = {}, const T &x = {})
      : data(mdsize<N>(sz), x),
        view(data.data(), sz) {}

  template<integral... inds, Index M = sizeof...(inds)>
  constexpr matpack_data(inds... sz) requires(M == N) : matpack_data(std::array<Index, N>{static_cast<Index>(std::forward<inds>(sz))...}) {}

  template <typename... arguments>
  constexpr matpack_data(integral auto dim0, arguments ... args)
    requires(sizeof...(arguments) == N)
      : matpack_data(front_array<Index, N>(dim0, args...), get_last(args...)) {
  }

  matpack_data(const matpack_data& x) : data(x.data), view(data.data(), x.shape()) {}

  constexpr matpack_data(const me_view& x) : matpack_data(x.shape()) {view = x;}
  constexpr matpack_data(const ce_view& x) : matpack_data(x.shape()) {view = x;}
  constexpr matpack_data(const ms_view& x) : matpack_data(x.shape()) {view = x;}
  constexpr matpack_data(const cs_view& x) : matpack_data(x.shape()) {view = x;}

  template<Index M> explicit constexpr matpack_data(const matpack::matpack_data<T, M>& x) requires(M < N) : matpack_data(upview<N, M>(x.shape())) { view = view_type{x}; }

  [[nodiscard]] constexpr auto size() const noexcept { return view.size(); }
  [[nodiscard]] constexpr auto extent(Index i) const { return view.extent(i); }
  [[nodiscard]] constexpr auto stride(Index i) const { return view.stride(i); }
  [[nodiscard]] constexpr auto shape() const { return view.shape(); }
  [[nodiscard]] constexpr auto inner_shape() const { return view.inner_shape(); }
  [[nodiscard]] constexpr auto strides() const { return view.strides(); }
  [[nodiscard]] constexpr auto inner_strides() const { return view.inner_strides(); }
  [[nodiscard]] constexpr auto inner_map() const { return view.inner_map(); }
  [[nodiscard]] constexpr auto nelem() const requires(N == 1) { return view.nelem(); }
  [[nodiscard]] constexpr auto ncols() const requires(N >= 2) { return view.ncols(); }
  [[nodiscard]] constexpr auto nrows() const requires(N >= 2) { return view.nrows(); }
  [[nodiscard]] constexpr auto npages() const requires(N >= 3) { return view.npages(); }
  [[nodiscard]] constexpr auto nbooks() const requires(N >= 4) { return view.nbooks(); }
  [[nodiscard]] constexpr auto nshelves() const requires(N >= 5) { return view.nshelves(); }
  [[nodiscard]] constexpr auto nvitrines() const requires(N >= 6) { return view.nvitrines(); }
  [[nodiscard]] constexpr auto nlibraries() const requires(N >= 7) { return view.nlibraries(); }
  [[nodiscard]] constexpr auto empty() const noexcept { return data.empty(); }

  constexpr void resize(const std::array<Index, N>& sz) { if (shape() not_eq sz) *this = matpack_data(sz); }

  template<integral... inds, Index M = sizeof...(inds)> constexpr void resize(inds&&... sz) requires(M == N) {resize(std::array<Index, N>{static_cast<Index>(std::forward<inds>(sz))...});}
  
  template<Index M> constexpr matpack_data<T, M> reshape(const std::array<Index, M>& sz) && {
    using other_view_type = typename matpack_data<T, M>::view_type;

    if (size() != mdsize<M>(sz)) std::terminate();
    ARTS_ASSERT(size() == mdsize<M>(sz), size(), " vs ", mdsize<M>(sz))
    
    matpack_data<T, M> out;
    out.data = std::move(data);
    out.view.secret_set(other_view_type{out.data.data(), sz});
    
    view.secret_set(view_type{nullptr, constant_array<N, 0>()});
    return out;
  }

  template<integral... inds, Index M = sizeof...(inds)> constexpr auto reshape(inds&&... sz) && { return std::move(*this).template reshape<M>(std::array<Index, M>{static_cast<Index>(std::forward<inds>(sz))...}); }

  template<Index... inds> constexpr auto reduce_rank() && requires(sizeof...(inds) < N ) { return std::move(*this).reshape(extent(inds)...); }

  constexpr matpack_data<T, 1> flatten() && { return std::move(*this).reshape(size()); }

  constexpr matpack_data& operator=(const matpack_data& x) {
    if (this not_eq &x) {
      resize(x.shape());
      view = x.view;
    }
    return *this;
  }

  constexpr matpack_data(matpack_data&& x) noexcept : data(std::move(x.data)), view(data.data(), x.shape()) {}

  constexpr matpack_data& operator=(matpack_data&& x) noexcept {
    data = std::move(x.data);
    view.secret_set(view_type(data.data(), x.shape()));
    x.view.secret_set(view_type{nullptr, constant_array<N, 0>()});
    return *this;
  }

  template <access_operator... access> constexpr auto operator()(access&&... ind) -> decltype(view(std::forward<access>(ind)...)) {return view(std::forward<access>(ind)...); }
  template <access_operator... access> constexpr auto operator()(access&&... ind) const -> decltype(view(std::forward<access>(ind)...)) const {return view(std::forward<access>(ind)...); }
  template <access_operator access> constexpr auto operator[](access&& ind) -> decltype(view[std::forward<access>(ind)]) {return view[std::forward<access>(ind)]; }
  template <access_operator access> constexpr auto operator[](access&& ind) const -> decltype(view[std::forward<access>(ind)]) const {return view[std::forward<access>(ind)]; }

  using value_type = T;
  [[nodiscard]] static constexpr auto rank() { return N; }
  [[nodiscard]] static constexpr auto is_const() { return false; }
  [[nodiscard]] constexpr auto data_handle() const { return view.data_handle(); }
  [[nodiscard]] constexpr auto unsafe_data_handle() const { return view.unsafe_data_handle(); }

  [[nodiscard]] constexpr auto begin() {if constexpr (N == 1) return data.begin(); else return view.begin();}
  [[nodiscard]] constexpr auto end() {if constexpr (N == 1) return data.end(); else return view.end();}
  [[nodiscard]] constexpr auto begin() const {if constexpr (N == 1) return data.begin(); else return view.begin();}
  [[nodiscard]] constexpr auto end() const {if constexpr (N == 1) return data.end(); else return view.end();}
  [[nodiscard]] constexpr auto cbegin() const {return begin();}
  [[nodiscard]] constexpr auto cend() const {return end();}

  [[nodiscard]] constexpr T& elem_at(Index ind) { return data[ind]; }
  [[nodiscard]] constexpr const T& elem_at(Index ind) const { return data[ind]; }
  [[nodiscard]] constexpr T& elem_at(const std::array<Index, N>& pos) { return view.elem_at(pos); }
  [[nodiscard]] constexpr const T& elem_at(const std::array<Index, N>& pos) const { return view.elem_at(pos); }
   
  [[nodiscard]] constexpr auto elem_begin() { return data.begin(); }
  [[nodiscard]] constexpr auto elem_end() { return data.end(); }
  [[nodiscard]] constexpr auto elem_begin() const { return data.begin(); }
  [[nodiscard]] constexpr auto elem_end() const { return data.end(); }
  [[nodiscard]] constexpr auto elem_cbegin() const {return elem_begin();}
  [[nodiscard]] constexpr auto elem_cend() const {return elem_end();}

  friend std::ostream& operator<<(std::ostream& os, const matpack_data& md) {
    return os << md.view;
  }

  [[nodiscard]] constexpr auto real() requires(complex_type<T>) { return view.real(); }
  [[nodiscard]] constexpr auto imag() requires(complex_type<T>) { return view.imag(); }
  [[nodiscard]] constexpr auto real() const requires(complex_type<T>) { return view.real(); }
  [[nodiscard]] constexpr auto imag() const requires(complex_type<T>) { return view.imag(); }
  [[nodiscard]] constexpr auto diagonal() requires(N == 2) { return view.diagonal(); }
  [[nodiscard]] constexpr auto diagonal() const requires(N == 2) { return matpack_view<T, 1, true, false>{view}.diagonal(); }

  constexpr matpack_data& operator=(std::convertible_to<T> auto x) {view = x; return *this;}
  constexpr matpack_data& operator+=(std::convertible_to<T> auto x) {view += x; return *this;}
  constexpr matpack_data& operator-=(std::convertible_to<T> auto x) {view -= x; return *this;}
  constexpr matpack_data& operator*=(std::convertible_to<T> auto x) {view *= x; return *this;}
  constexpr matpack_data& operator/=(std::convertible_to<T> auto x) {view /= x; return *this;}
  
  template <bool c, bool s> constexpr matpack_data& operator+=(const matpack_view<T, N, c, s>& x) {view += x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator-=(const matpack_view<T, N, c, s>& x) {view -= x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator*=(const matpack_view<T, N, c, s>& x) {view *= x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator/=(const matpack_view<T, N, c, s>& x) {view /= x; return *this;}
  constexpr matpack_data& operator+=(const matpack_data& x) {view += x.view; return *this;}
  constexpr matpack_data& operator-=(const matpack_data& x) {view -= x.view; return *this;}
  constexpr matpack_data& operator*=(const matpack_data& x) {view *= x.view; return *this;}
  constexpr matpack_data& operator/=(const matpack_data& x) {view /= x.view; return *this;}

  template <bool c, bool s> constexpr bool operator==(const matpack_view<T, N, c, s>& x) const { return view == x; }
  template <bool c, bool s> constexpr bool operator!=(const matpack_view<T, N, c, s>& x) const { return view != x; }
  constexpr bool operator==(const matpack_data& x) const { return view == x.view; }
  constexpr bool operator!=(const matpack_data& x) const { return view != x.view; }

  constexpr void swap(matpack_data& other) noexcept {
    data.swap(other.data);
    view.swap(other.view);
  }

  template <matpack_convertible<T, N> U>
  constexpr matpack_data(const U& x) : matpack_data(mdshape(x)) { view = x; }

  template <matpack_convertible<T, N> U>
  constexpr matpack_data& operator=(const U& x) {
    *this = matpack_data(x);
    return *this;
  }

  constexpr matpack_data(std::vector<T>&& a) requires(N == 1) : data(std::move(a)), view(data.data(), std::array<Index, N>{static_cast<Index>(data.size())}) {}
  constexpr matpack_data(std::initializer_list<T> a) requires(N == 1) : matpack_data(std::vector<T>{a}) {}

  static constexpr bool is_always_exhaustive() noexcept { return true; }
};
} // namespace matpack

template <typename T, Index N>
std::string describe(const matpack::matpack_data<T, N>& m) {
  using namespace matpack;
  return var_string("matpack_data of rank ", N, " of shape ", shape_help<N>(m.shape()));
}

using Range = matpack::matpack_strided_access;

using Vector = matpack::matpack_data<Numeric, 1>;
using Matrix = matpack::matpack_data<Numeric, 2>;
using Tensor3 = matpack::matpack_data<Numeric, 3>;
using Tensor4 = matpack::matpack_data<Numeric, 4>;
using Tensor5 = matpack::matpack_data<Numeric, 5>;
using Tensor6 = matpack::matpack_data<Numeric, 6>;
using Tensor7 = matpack::matpack_data<Numeric, 7>;

using VectorView = matpack::matpack_view<Numeric, 1, false, true>;
using MatrixView = matpack::matpack_view<Numeric, 2, false, true>;
using Tensor3View = matpack::matpack_view<Numeric, 3, false, true>;
using Tensor4View = matpack::matpack_view<Numeric, 4, false, true>;
using Tensor5View = matpack::matpack_view<Numeric, 5, false, true>;
using Tensor6View = matpack::matpack_view<Numeric, 6, false, true>;
using Tensor7View = matpack::matpack_view<Numeric, 7, false, true>;

using ConstVectorView = matpack::matpack_view<Numeric, 1, true, true>;
using ConstMatrixView = matpack::matpack_view<Numeric, 2, true, true>;
using ConstTensor3View = matpack::matpack_view<Numeric, 3, true, true>;
using ConstTensor4View = matpack::matpack_view<Numeric, 4, true, true>;
using ConstTensor5View = matpack::matpack_view<Numeric, 5, true, true>;
using ConstTensor6View = matpack::matpack_view<Numeric, 6, true, true>;
using ConstTensor7View = matpack::matpack_view<Numeric, 7, true, true>;

using ExhaustiveVectorView = matpack::matpack_view<Numeric, 1, false, false>;
using ExhaustiveMatrixView = matpack::matpack_view<Numeric, 2, false, false>;
using ExhaustiveTensor3View = matpack::matpack_view<Numeric, 3, false, false>;
using ExhaustiveTensor4View = matpack::matpack_view<Numeric, 4, false, false>;
using ExhaustiveTensor5View = matpack::matpack_view<Numeric, 5, false, false>;
using ExhaustiveTensor6View = matpack::matpack_view<Numeric, 6, false, false>;
using ExhaustiveTensor7View = matpack::matpack_view<Numeric, 7, false, false>;

using ExhaustiveConstVectorView = matpack::matpack_view<Numeric, 1, true, false>;
using ExhaustiveConstMatrixView = matpack::matpack_view<Numeric, 2, true, false>;
using ExhaustiveConstTensor3View = matpack::matpack_view<Numeric, 3, true, false>;
using ExhaustiveConstTensor4View = matpack::matpack_view<Numeric, 4, true, false>;
using ExhaustiveConstTensor5View = matpack::matpack_view<Numeric, 5, true, false>;
using ExhaustiveConstTensor6View = matpack::matpack_view<Numeric, 6, true, false>;
using ExhaustiveConstTensor7View = matpack::matpack_view<Numeric, 7, true, false>;

using ComplexVector = matpack::matpack_data<Complex, 1>;
using ComplexMatrix = matpack::matpack_data<Complex, 2>;
using ComplexTensor3 = matpack::matpack_data<Complex, 3>;
using ComplexTensor4 = matpack::matpack_data<Complex, 4>;
using ComplexTensor5 = matpack::matpack_data<Complex, 5>;
using ComplexTensor6 = matpack::matpack_data<Complex, 6>;
using ComplexTensor7 = matpack::matpack_data<Complex, 7>;

using ComplexVectorView = matpack::matpack_view<Complex, 1, false, true>;
using ComplexMatrixView = matpack::matpack_view<Complex, 2, false, true>;
using ComplexTensor3View = matpack::matpack_view<Complex, 3, false, true>;
using ComplexTensor4View = matpack::matpack_view<Complex, 4, false, true>;
using ComplexTensor5View = matpack::matpack_view<Complex, 5, false, true>;
using ComplexTensor6View = matpack::matpack_view<Complex, 6, false, true>;
using ComplexTensor7View = matpack::matpack_view<Complex, 7, false, true>;

using ConstComplexVectorView = matpack::matpack_view<Complex, 1, true, true>;
using ConstComplexMatrixView = matpack::matpack_view<Complex, 2, true, true>;
using ConstComplexTensor3View = matpack::matpack_view<Complex, 3, true, true>;
using ConstComplexTensor4View = matpack::matpack_view<Complex, 4, true, true>;
using ConstComplexTensor5View = matpack::matpack_view<Complex, 5, true, true>;
using ConstComplexTensor6View = matpack::matpack_view<Complex, 6, true, true>;
using ConstComplexTensor7View = matpack::matpack_view<Complex, 7, true, true>;

using ExhaustiveComplexVectorView = matpack::matpack_view<Complex, 1, false, false>;
using ExhaustiveComplexMatrixView = matpack::matpack_view<Complex, 2, false, false>;
using ExhaustiveComplexTensor3View = matpack::matpack_view<Complex, 3, false, false>;
using ExhaustiveComplexTensor4View = matpack::matpack_view<Complex, 4, false, false>;
using ExhaustiveComplexTensor5View = matpack::matpack_view<Complex, 5, false, false>;
using ExhaustiveComplexTensor6View = matpack::matpack_view<Complex, 6, false, false>;
using ExhaustiveComplexTensor7View = matpack::matpack_view<Complex, 7, false, false>;

using ExhaustiveConstComplexVectorView = matpack::matpack_view<Complex, 1, true, false>;
using ExhaustiveConstComplexMatrixView = matpack::matpack_view<Complex, 2, true, false>;
using ExhaustiveConstComplexTensor3View = matpack::matpack_view<Complex, 3, true, false>;
using ExhaustiveConstComplexTensor4View = matpack::matpack_view<Complex, 4, true, false>;
using ExhaustiveConstComplexTensor5View = matpack::matpack_view<Complex, 5, true, false>;
using ExhaustiveConstComplexTensor6View = matpack::matpack_view<Complex, 6, true, false>;
using ExhaustiveConstComplexTensor7View = matpack::matpack_view<Complex, 7, true, false>;
