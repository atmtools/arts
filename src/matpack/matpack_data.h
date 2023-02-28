#pragma once

#include "matpack_complex.h"
#include "matpack_concepts.h"
#include "matpack_view.h"

#include <exception>
#include <tuple>

namespace matpack {
namespace detail {
//! Sets I in out to v --- implementation detail...
template <std::size_t I=0, typename T, std::size_t N, typename value, typename... rest>
void front_array_impl(std::array<T, N> &out, value &&v, rest &&...x) {
  if constexpr (I < N) {
    std::get<I>(out) = static_cast<T>(std::forward<value>(v));
    front_array_impl<I + 1>(out, std::forward<rest>(x)...);
  }
}
} // namespace detail

//! Extracts the first N elements out of a variadic input by statically casting
//! them to T
template <typename T, Index N, typename... input, Index M = sizeof...(input)>
constexpr std::array<T, N> front_array(input... in)
  requires(N <= M and N > 0)
{
  using namespace detail;
  std::array<T, N> out;
  front_array_impl(out, in...);
  return out;
}

//! Get the last element out of a variadic input
template<typename... rest, Index N=sizeof...(rest)>
auto get_last(rest&&... r) requires(N>0) {
  return std::get<N-1>(std::tuple{std::forward<rest>(r)...});
}

//! The basic data type
template <typename T, Index N>
class matpack_data {
  //! Allocates all the data
  std::vector<T> data;

  //! The basic type of which we view this data
  using view_type = matpack_view<T, N, false, false>;

  //! The view of the data
  view_type view;

  //! A mutable exhaustive matpack_view
  using me_view = typename view_type::me_view;

  //! A mutable strided matpack_view
  using ms_view = typename view_type::ms_view;

  //! A constant exhaustive matpack_view
  using ce_view = typename view_type::ce_view;

  //! A constant strided matpack_view
  using cs_view = typename view_type::cs_view;

  //! Any view type
  template <bool c, bool s>
  using any_view = matpack_view<T, N, c, s>;

  //! Allow all other matpack_data to view private information about this type
  template <typename U, Index M> friend class matpack_data;

  //! Allow all matpack views to view private information about this class
  template <typename U, Index M, bool c, bool s> friend class matpack_view;

public:
  constexpr operator view_type&() {return view;}
  constexpr operator const view_type&() const {return view;}

  /** @brief Core constructor
   * 
   * @param[in] sz The size object determining the extent of each dimension
   * @param[in] x The default value of all input
   */
  constexpr matpack_data(const std::array<Index, N> &sz = {}, const T &x = {})
      : data(mdsize<N>(sz), x),
        view(data.data(), sz) {}

  /** Variadic input of the size of this object
   *
   * The object is default initialized to T{}
   * 
   * @param[in] sz Sizes to determine the extent of each dimension
   */
  template<integral... inds, Index M = sizeof...(inds)>
  constexpr matpack_data(inds... sz) requires(M == N) : matpack_data(std::array<Index, N>{static_cast<Index>(std::forward<inds>(sz))...}) {}

  /** Variadic input of the size and default value of this object
   * 
   * @param[in] dim0 Left-most dimension extent
   * @param[in] args The extent of all other dimensions followed by the initial value
   */
  template <typename... arguments>
  constexpr matpack_data(integral auto dim0, arguments ... args)
    requires(sizeof...(arguments) == N)
      : matpack_data(front_array<Index, N>(dim0, args...), get_last(args...)) {
  }

  /** Construct a new matpack data object
   * 
   * Needs to be manual as the view should view the copied data
   * and not the old data's view
   *
   * @param[in] x An old matpack data object
   */
  matpack_data(const matpack_data& x) : data(x.data), view(data.data(), x.shape()) {}

  /** Construct a new matpack data object
   *
   * @param[in] x An old matpack view object
   */
  explicit constexpr matpack_data(const me_view& x) : matpack_data(x.shape()) {view = x;}

  /** Construct a new matpack data object
   *
   * @param[in] x An old matpack view object
   */
  explicit constexpr matpack_data(const ce_view& x) : matpack_data(x.shape()) {view = x;}

  /** Construct a new matpack data object
   *
   * @param[in] x An old matpack view object
   */
  explicit constexpr matpack_data(const ms_view& x) : matpack_data(x.shape()) {view = x;}

  /** Construct a new matpack data object
   *
   * @param[in] x An old matpack view object
   */
  explicit constexpr matpack_data(const cs_view& x) : matpack_data(x.shape()) {view = x;}

  /** Construct a new matpack data object
   *
   * The input is a smaller matpack data object, and the right-most
   * dimensions of the new object will be all ones whereas the left-most
   * extents will be the same as the old object's extents
   *
   * @param[in] x An old matpack view object that is smaller
   */
  template <Index M>
  explicit constexpr matpack_data(const matpack::matpack_data<T, M> &x)
    requires(M < N)
      : matpack_data(upview<N, M>(x.shape())) {
    view = view_type{x};
  }

  //! Return the total number of elements that can be accessed in this object
  [[nodiscard]] constexpr auto size() const noexcept { return view.size(); }

  //! Return the extent of dimenion i
  [[nodiscard]] constexpr auto extent(Index i) const { return view.extent(i); }

  //! Return the stride of dimenion i
  [[nodiscard]] constexpr auto stride(Index i) const { return view.stride(i); }

  //! Return the shape of the object
  [[nodiscard]] constexpr auto shape() const { return view.shape(); }

  //! Return the shape of the object skipping the left-most extent
  [[nodiscard]] constexpr auto inner_shape() const { return view.inner_shape(); }

  //! Return the strides of the object
  [[nodiscard]] constexpr auto strides() const { return view.strides(); }

  //! Return the strides of the object skipping the left-most stride
  [[nodiscard]] constexpr auto inner_strides() const { return view.inner_strides(); }

  //! Return a full mdspan mapping of this object
  [[nodiscard]] constexpr auto inner_map() const { return view.inner_map(); }

  //! Return the size of the 1-dimensional object
  [[nodiscard]] constexpr auto nelem() const requires(N == 1) { return view.nelem(); }

  //! Return the right-most extent of the multidimensional object
  [[nodiscard]] constexpr auto ncols() const requires(N >= 2) { return view.ncols(); }

  //! Return the right-most bar 1 extent of the multidimensional object
  [[nodiscard]] constexpr auto nrows() const requires(N >= 2) { return view.nrows(); }

  //! Return the right-most bar 2 extent of the multidimensional object
  [[nodiscard]] constexpr auto npages() const requires(N >= 3) { return view.npages(); }

  //! Return the right-most bar 3 extent of the multidimensional object
  [[nodiscard]] constexpr auto nbooks() const requires(N >= 4) { return view.nbooks(); }

  //! Return the right-most bar 4 extent of the multidimensional object
  [[nodiscard]] constexpr auto nshelves() const requires(N >= 5) { return view.nshelves(); }

  //! Return the right-most bar 5 extent of the multidimensional object
  [[nodiscard]] constexpr auto nvitrines() const requires(N >= 6) { return view.nvitrines(); }

  //! Return the right-most bar 6 extent of the multidimensional object
  [[nodiscard]] constexpr auto nlibraries() const requires(N >= 7) { return view.nlibraries(); }

  //! Return if the object contains no data
  [[nodiscard]] constexpr auto empty() const noexcept { return data.empty(); }

  //! Resize this object and default constructs the new one IF the new shape is different
  constexpr void resize(const std::array<Index, N>& sz) { if (shape() not_eq sz) *this = matpack_data(sz); }

  //! Resize this object and default constructs the new one IF the new shape is different
  template<integral... inds, Index M = sizeof...(inds)> constexpr void resize(inds&&... sz) requires(M == N) {resize(std::array<Index, N>{static_cast<Index>(std::forward<inds>(sz))...});}
  
  //! Reshape this object to another size of matpack data.  The new object must have the same size as the old one had
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

  //! Reshape this object to another size of matpack data.  The new object must have the same size as the old one had
  template<integral... inds, Index M = sizeof...(inds)> constexpr auto reshape(inds&&... sz) && { return std::move(*this).template reshape<M>(std::array<Index, M>{static_cast<Index>(std::forward<inds>(sz))...}); }

  //! Reduce the rank of the object.  The new object must have the same size as the old one had
  template<Index... inds> constexpr auto reduce_rank() && requires(sizeof...(inds) < N ) { return std::move(*this).reshape(extent(inds)...); }

  //! Return this object as a 1-dimensional object
  constexpr matpack_data<T, 1> flatten() && { return std::move(*this).reshape(size()); }

  constexpr matpack_data& operator=(const matpack_data& x) {
    if (this not_eq &x) {
      resize(x.shape());
      view = x.view;
    }
    return *this;
  }

  /** Construct a new matpack data object
   * 
   * Needs to be manual as the view should view the moved data
   * and not the old data's view
   *
   * @param[in] x An old matpack data object
   */
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

  //! The value type of this matpack data is public information
  using value_type = T;

  //! Return the rank of this object
  [[nodiscard]] static constexpr auto rank() { return N; }

  //! Return the constness of this object
  [[nodiscard]] static constexpr auto is_const() { return false; }

  //! Return a pointer to this object's allocated data
  [[nodiscard]] constexpr auto data_handle() const { return view.data_handle(); }

  //! Return a pointer to this object's allocated data
  [[nodiscard]] constexpr auto unsafe_data_handle() const { return view.unsafe_data_handle(); }

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto begin() {if constexpr (N == 1) return data.begin(); else return view.begin();}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto end() {if constexpr (N == 1) return data.end(); else return view.end();}

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto begin() const {if constexpr (N == 1) return data.begin(); else return view.begin();}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto end() const {if constexpr (N == 1) return data.end(); else return view.end();}

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto cbegin() const {return begin();}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto cend() const {return end();}

  //! Access the data elements regardless of the matpack data rank
  template <integral... Inds>
  [[nodiscard]] constexpr T& elem_at(Inds... inds) requires(sizeof...(Inds) == N) {
    if constexpr (N == 1) return data[std::get<0>(std::tuple{inds...})];
    else return view(inds...);
  }

  //! Access the data elements regardless of the matpack data rank  template <integral... Inds>
  template <integral... Inds>
  [[nodiscard]] constexpr T elem_at(Inds... inds) const requires(sizeof...(Inds) == N) {
    if constexpr (N == 1) return data[std::get<0>(std::tuple{inds...})];
    else return view(inds...);
  }
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_begin() { return data.begin(); }

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_end() { return data.end(); }
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_begin() const { return data.begin(); }

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_end() const { return data.end(); }
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_cbegin() const {return elem_begin();}

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_cend() const {return elem_end();}

  friend std::ostream& operator<<(std::ostream& os, const matpack_data& md) {
    return os << md.view;
  }

  //! Return a view of the real part of this object
  [[nodiscard]] constexpr auto real() requires(complex_type<T>) { return view.real(); }

  //! Return a view of the imaginary part of this object
  [[nodiscard]] constexpr auto imag() requires(complex_type<T>) { return view.imag(); }

  //! Return a view of the real part of this object
  [[nodiscard]] constexpr auto real() const requires(complex_type<T>) { return view.real(); }

  //! Return a view of the imaginary part of this object
  [[nodiscard]] constexpr auto imag() const requires(complex_type<T>) { return view.imag(); }

  //! Return a view of the diagonal of this object
  [[nodiscard]] constexpr auto diagonal() requires(N == 2) { return view.diagonal(); }

  //! Return a view of the diagonal of this object
  [[nodiscard]] constexpr auto diagonal() const requires(N == 2) { return matpack_view<T, 1, true, false>{view}.diagonal(); }

  constexpr matpack_data& operator=(std::convertible_to<T> auto x) {std::fill(data.begin(), data.end(), static_cast<T>(x)); return *this;}
  constexpr matpack_data& operator+=(std::convertible_to<T> auto x) {std::transform(data.begin(), data.end(), data.begin(), [y=static_cast<T>(x)](auto z){return z + y;}); return *this;}
  constexpr matpack_data& operator-=(std::convertible_to<T> auto x) {std::transform(data.begin(), data.end(), data.begin(), [y=static_cast<T>(x)](auto z){return z - y;}); return *this;}
  constexpr matpack_data& operator*=(std::convertible_to<T> auto x) {std::transform(data.begin(), data.end(), data.begin(), [y=static_cast<T>(x)](auto z){return z * y;}); return *this;}
  constexpr matpack_data& operator/=(std::convertible_to<T> auto x) {std::transform(data.begin(), data.end(), data.begin(), [y=static_cast<T>(x)](auto z){return z / y;}); return *this;}
  
  template <bool c, bool s> constexpr matpack_data& operator=(const any_view<c, s>& x) {
    if constexpr (not s) {
      if (data_handle() != x.unsafe_data_handle()) {
        resize(x.shape());
        view = x;
      }
    } else {
      *this = matpack_data{x};
    }
    return *this;
  }

  template <bool c, bool s> constexpr matpack_data& operator+=(const any_view<c, s>& x) {view += x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator-=(const any_view<c, s>& x) {view -= x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator*=(const any_view<c, s>& x) {view *= x; return *this;}
  template <bool c, bool s> constexpr matpack_data& operator/=(const any_view<c, s>& x) {view /= x; return *this;}
  constexpr matpack_data& operator+=(const matpack_data& x) {view += x.view; return *this;}
  constexpr matpack_data& operator-=(const matpack_data& x) {view -= x.view; return *this;}
  constexpr matpack_data& operator*=(const matpack_data& x) {view *= x.view; return *this;}
  constexpr matpack_data& operator/=(const matpack_data& x) {view /= x.view; return *this;}

  template <bool c, bool s> constexpr bool operator==(const any_view<c, s>& x) const { return view == x; }
  template <bool c, bool s> constexpr bool operator!=(const any_view<c, s>& x) const { return view != x; }
  constexpr bool operator==(const matpack_data& x) const { return view == x.view; }
  constexpr bool operator!=(const matpack_data& x) const { return view != x.view; }

  //! Swaps the data around
  constexpr void swap(matpack_data& other) noexcept {
    data.swap(other.data);
    view.swap(other.view);
  }

  //! Construct this object from another that is conceptually convertible to this object's type
  template <matpack_convertible<T, N> U>
  explicit constexpr matpack_data(const U& x) : matpack_data(mdshape(x)) { view = x; }

  template <matpack_convertible<T, N> U>
  constexpr matpack_data& operator=(const U& x) {
    *this = matpack_data(x);
    return *this;
  }

  //! Allow a specialization to construct this object from a standard vector
  constexpr matpack_data(std::vector<T>&& a) requires(N == 1) : data(std::move(a)), view(data.data(), std::array<Index, N>{static_cast<Index>(data.size())}) {}

  //! Allow a specialization to construct this object from a standard initializer list
  constexpr matpack_data(std::initializer_list<T> a) requires(N == 1) : matpack_data(std::vector<T>{a}) {}

  //! Return that this object is always exhaustive
  static constexpr bool is_always_exhaustive() noexcept { return true; }
};

/** Describe the matpack data type
 * 
 * @param m Any matpack_data type
 * @return std::string of the description
 */
template <typename T, Index N>
std::string describe(const matpack_data<T, N>& m) {
  using namespace matpack;
  return var_string("matpack_data of rank ", N, " of shape ", shape_help<N>(m.shape()));
}
} // namespace matpack

//! A vector of Numeric
using Vector = matpack::matpack_data<Numeric, 1>;

//! A matrix of Numeric
using Matrix = matpack::matpack_data<Numeric, 2>;

//! A tensor of Numeric of rank 3
using Tensor3 = matpack::matpack_data<Numeric, 3>;

//! A tensor of Numeric of rank 4
using Tensor4 = matpack::matpack_data<Numeric, 4>;

//! A tensor of Numeric of rank 5
using Tensor5 = matpack::matpack_data<Numeric, 5>;

//! A tensor of Numeric of rank 6
using Tensor6 = matpack::matpack_data<Numeric, 6>;

//! A tensor of Numeric of rank 7
using Tensor7 = matpack::matpack_data<Numeric, 7>;

//! A vector of Complex
using ComplexVector = matpack::matpack_data<Complex, 1>;

//! A matrix of Complex
using ComplexMatrix = matpack::matpack_data<Complex, 2>;
