#pragma once

#include "matpack_concepts.h"
#include "matpack_iter.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <debug.h>
#include <functional>
#include <numeric>
#include <ostream>
#include <string>
#include <type_traits>

namespace matpack {
//! A strided accessor to get the range of a type
struct matpack_strided_access {
  //! Offset of the object; must be positive or 0
  Index offset{0};

  //! Extent of the object; negative means all
  Index extent{-1};

  //! Stride of the object
  Index stride{1};

  //! Construct this object with an offset, extent, and stride
  constexpr matpack_strided_access(Index i0=0, Index len=-1, Index d=1) noexcept : offset(i0), extent(len), stride(d) {}

  //! Construct this object with an offset, full extent, and stride
  constexpr matpack_strided_access(Index i0, Joker, Index d=1) noexcept : offset(i0), stride(d) {}

  //! Construct this object with full extent and a stride
  constexpr matpack_strided_access(Joker, Index d=1) noexcept : stride(d) {}

  //! Construct this object with a max size from another strided access
  constexpr matpack_strided_access(Index max_size, const matpack_strided_access& r)  ARTS_NOEXCEPT
      : offset(r.offset),
        extent(r.extent),
        stride(r.stride) {
    // Start must be >= 0:
    ARTS_ASSERT(0 <= offset, "offset=", offset);
    // ... and < max_size:
    ARTS_ASSERT(offset < max_size or (max_size == 0 and offset == 0), "offset=", offset, "; max_size=", max_size);

    // Stride must be != 0:
    ARTS_ASSERT(0 != stride, "stride=", stride);

    // Convert negative extent (joker) to explicit extent
    if (extent < 0) {
      if (0 < stride)
        extent = 1 + (max_size - 1 - offset) / stride;
      else
        extent = 1 + (0 - offset) / stride;
    } else {
      // Check that extent is ok:
      ARTS_ASSERT(0 <= (offset + (extent - 1) * stride));
      ARTS_ASSERT((offset + (extent - 1) * stride) < max_size);
    }
  }

  //! Combine two strided access objects
  constexpr matpack_strided_access(const matpack_strided_access& p, const matpack_strided_access& n) ARTS_NOEXCEPT
      : offset(p.offset + n.offset * p.stride),
        extent(n.extent),
        stride(p.stride* n.stride) {
    // We have to juggle here a bit with previous, new, and resulting
    // quantities. I.e.;
    // p.mstride: Previous stride
    // n.mstride: New stride (as specified)
    // mstride:   Resulting stride (old*new)

    // Get the previous final element:
    Index const prev_fin = p.offset + (p.extent - 1) * p.stride;

    // Resulting start must be >= previous start:
    ARTS_ASSERT(p.offset <= offset);
    // and <= prev_fin, except for Joker:
    ARTS_ASSERT(offset <= prev_fin || extent == -1);

    // Resulting stride must be != 0:
    ARTS_ASSERT(0 != stride,
                "Problem: stride==",
                stride,
                " for offset==",
                offset,
                " and extent==",
                extent);

    // Convert negative extent (joker) to explicit extent
    if (extent < 0) {
      if (0 < stride)
        extent = 1 + (prev_fin - offset) / stride;
      else
        extent = 1 + (p.offset - offset) / stride;
    } else {
      ARTS_ASSERT(p.offset <= (offset + (extent - 1) * stride));
      ARTS_ASSERT((offset + (extent - 1) * stride) <= prev_fin);
    }
  }

  //! FIXME: Remove once mdspan has better type system 
  constexpr operator Joker() const noexcept {return joker;}

  //! FIXME: Remove once mdspan has better type system
  constexpr operator Joker() noexcept {return joker;}
  
  friend std::ostream &operator<<(std::ostream &os,
                                  const matpack_strided_access &msa) {
    return os << "Range(" << msa.offset << ", " << msa.extent << ", "
              << msa.stride << ")";
  }
};

//! Count the number of integral there are in the access operators
template <access_operator... access, Index N=sizeof...(access)>
[[nodiscard]] constexpr Index index_access_operator_count() noexcept {
  constexpr std::array<bool, N> indices{integral<access>...};
  return std::reduce(indices.begin(), indices.end(), Index{0});
}

//! Helper index
template <access_operator... access>
inline constexpr Index num_index = index_access_operator_count<access...>();

//! Check if the access operators contain any strided access
template <access_operator... access, Index N=sizeof...(access)>
constexpr bool any_range() {
  constexpr std::array<bool, N> ranges{is_matpack_strided_access<access>...};
  return std::any_of(ranges.begin(), ranges.end(), [](auto b){return b;});
};

//! Helper bool
template <access_operator... access>
inline constexpr bool has_any_range = any_range<access...>();

//! Check if all the indices of these access operators are to the left
template <access_operator... access, Index N=sizeof...(access)>
constexpr bool all_index_left() {
  constexpr std::array<bool, N> integrals{not integral<access>...};
  return std::is_sorted(integrals.begin(), integrals.end());
}

/** Check that the access operation is not out-of-bounds
 * 
 * @param arr The original extents
 * @param r The current extent position
 * @param first The way the current extent is accessed
 * @param rest All following access operators
 * @return true if the access operators are OK
 * @return false otherwise
 */
template <access_operator T, access_operator ... Ts>
[[nodiscard]] constexpr bool check_index_sizes(const auto& arr, Index r, T first, Ts... rest) {
  constexpr bool test = integral<T>;
  if constexpr (test) {
    if (static_cast<Index>(first) >= static_cast<Index>(arr.extent(r))) return false;
  }

  if constexpr (sizeof...(Ts) == 0) {
    return true;
  } else {
    return check_index_sizes(arr, r+1, rest...);
  }
}

//! Helper namespace to determine the type of an access operation on a mdspan view
namespace access_retval {

//! Multiple access operators possible
template <typename T, Index N, bool constant, bool strided,
          access_operator... access>
struct helper {
  //! Count of how much the dimensionality should be reduced
  static constexpr Index M = num_index<access...>;

  //! Are there any strided access in this type?
  static constexpr bool ranged = has_any_range<access...>;

  //! Are all the index-wise access operations to the left?
  static constexpr bool left_access = all_index_left<access...>();

  using mutable_type =
      std::conditional_t<M == N, T &,
                         matpack_view<T, N - M, constant, ranged or strided or not left_access>>;

  using constant_type =
      std::conditional_t<M == N, const T &,
                         matpack_view<T, N - M, true, ranged or strided or not left_access>>;
};

//! Single access operator helper
template <typename T, Index N, bool constant, bool strided, access_operator access>
struct subhelper {
  //! Is the accessor a strided access?
  static constexpr bool ranged = has_any_range<access>;

  using mutable_type = std::conditional_t<
      ranged, matpack_view<T, N, constant, true>,
      std::conditional_t<N == 1, T &,
                         matpack_view<T, N - 1, constant, strided>>>;

  using constant_type = std::conditional_t<
      ranged, matpack_view<T, N, true, true>,
      std::conditional_t<N == 1, const T &,
                         matpack_view<T, N - 1, true, strided>>>;
};
} // namespace access_retval

//! Helper for when the access type is mutable
template <typename T, Index N, bool constant, bool strided, access_operator... access>
using mutable_access_type = typename access_retval::helper<T, N, constant, strided, access...>::mutable_type;

//! Helper for when the access type is constant
template <typename T, Index N, bool constant, bool strided, access_operator... access>
using constant_access_type = typename access_retval::helper<T, N, constant, strided, access...>::constant_type;

//! Helper for when the access type is left mutable
template <typename T, Index N, bool constant, bool strided, access_operator access>
using mutable_left_access_type = typename access_retval::subhelper<T, N, constant, strided, access>::mutable_type;

//! Helper for when the access type is constant
template <typename T, Index N, bool constant, bool strided, access_operator access>
using constant_left_access_type = typename access_retval::subhelper<T, N, constant, strided, access>::constant_type;

//! Helper struct to pretty-print the shape of the object
template <std::size_t N> struct shape_help {
std::array<Index, N> shape;
constexpr shape_help(const std::array<Index, N>& x) : shape(x) {}
friend std::ostream& operator<<(std::ostream& os, const shape_help& sh) {
  bool first=true;
  os << '(';
  for (auto& x: sh.shape) {
    if (not first) os << ", ";
    first = false;
    os << x;
  }
  return os << ')';
}
};

//! Return a constant Index array of size N with value v
template <Index N, Index v> 
constexpr std::array<Index, N> constant_array() {
  std::array<Index, N> x;
  x.fill(v);
  return x;
}

//! Helper to view a smaller size as a larger size
template <Index N, Index M>
std::array<Index, N> upview(const std::array<Index, M> &low) {
  std::array<Index, N> sz = constant_array<N, 1>();
  for (Index i = 0; i < M; i++)
  sz[i] = low[i];
  return sz;
}

//! The basic view type
template <typename T, Index N, bool constant, bool strided>
class matpack_view {
  static_assert(N > 0);

  //! The strided view type from mdspan
  using strided_view = strided_mdspan<T, N>;

  //! The exhaustive view type from mdspan
  using exhaustive_view = exhaustive_mdspan<T, N>;

  //! The view type of this type
  using view_type = std::conditional_t<strided, strided_view, exhaustive_view>;

  //! A mutable exhaustive version of this template class
  using me_view = matpack_view<T, N, false, false>;

  //! A const exhaustive version of this template class
  using ce_view = matpack_view<T, N, true, false>;

  //! A mutable strided version of this template class
  using ms_view = matpack_view<T, N, false, true>;

  //! A const strided version of this template class
  using cs_view = matpack_view<T, N, true, true>;

  //! Helper to return correct mutable type
  template <access_operator... access>
  using mutable_access = mutable_access_type<T, N, constant, strided, access...>;

  //! Helper to return correct constant type
  template <access_operator... access>
  using constant_access = constant_access_type<T, N, constant, strided, access...>;

  //! Helper to return correct mutable type
  template <access_operator access>
  using mutable_left_access = mutable_access_type<T, N, constant, strided, access>;

  //! Helper to return correct constant type
  template <access_operator access>
  using constant_left_access = constant_access_type<T, N, constant, strided, access>;

  //! The pure data view from mdspan
  view_type view;

  //! Allow any matpack data type to access the private parts of this object
  template <typename U, Index M> friend class matpack_data;

  //! Allow any other matpack view type to access the private parts of this object
  template <typename U, Index M, bool c, bool s> friend class matpack_view;

  //! Helper function to change the stride of this object's view
  template <Index dim>
  constexpr void strided_dimension(matpack_strided_access range) requires(dim >= 0 and dim < N and strided) {
    range = matpack_strided_access(extent(dim), range);
    auto *d = view.data_handle() + stride(dim) * range.offset;
    std::array<Index, N> extmap=shape(), strmap=strides();
    std::get<dim>(extmap) =
        range.extent < 0 ? std::get<dim>(extmap) : range.extent;
    std::get<dim>(strmap) *= range.stride;
    view = strided_view{d, {extmap, strmap}};
  }

  //! Helper function to change the stride of this object's view
  template <Index i = 0, access_operator first, access_operator... access>
  constexpr matpack_view &range_adaptor(first &&maybe_range, access &&...rest) requires(strided) {
    if constexpr (is_matpack_strided_access<first>)
      strided_dimension<i>(std::forward<first>(maybe_range));
    if constexpr (sizeof...(access) not_eq 0)
      range_adaptor<i + not integral<first>>(std::forward<access>(rest)...);
    return *this;
  }

  //! Sets this view to another view type, ignoring shape-checks
  constexpr void secret_set(view_type x) { view = std::move(x); }

  //! Allow getting the view type out of here in private contexts
  operator view_type() {return view;}

public:
  constexpr matpack_view() = default;

  //! Construct this from an exhaustive view
  constexpr matpack_view(exhaustive_view x) : view(std::move(x)) {}

  //! Construct this from an exhaustive view
  constexpr matpack_view(strided_view x) : view(std::move(x)) {}

  //! Construct this a pointer and a size
  constexpr matpack_view(T* ptr, std::array<Index, N> sz) : matpack_view(exhaustive_view{ptr, sz}) {}

  //! Construct this from another compatible view
  constexpr matpack_view(const me_view& x) : view(x.view) {}

  //! Construct this from another compatible view
  constexpr matpack_view(const ce_view& x) requires(constant) : view(x.view) {}

  //! Construct this from another compatible view
  constexpr matpack_view(const ms_view& x) requires(strided) : view(x.view) {}

  //! Construct this from another compatible view
  constexpr matpack_view(const cs_view& x) requires(constant and strided) : view(x.view) {}

  //! Construct this from the same size data view
  constexpr matpack_view(const matpack_data<T, N>& x) noexcept : matpack_view(x.view) {}

  //! Construct this from a smaller view by upping the rank with padded 1-extents to the right
  template <Index M> explicit constexpr matpack_view(const matpack_data<T, M>& x) requires(N > M) : matpack_view(x.data_handle(), upview<N, M>(x.shape())) {}

  //! Construct this by viewing another view with padded 1-extents to the right
  template <Index M> explicit matpack_view(const matpack_view<T, M, constant, strided>& x) requires(N > M) {
    using map_t = typename view_type::mapping_type;
    if constexpr (not strided) {
      view = view_type{x.data_handle(), map_t{upview<N, M>(x.shape())}};
    } else {
      view = view_type{x.unsafe_data_handle(), map_t{upview<N, M>(x.shape()), upview<N, M>(x.strides())}};
    }
  }

  //! Mutable view of x
  explicit constexpr matpack_view(T& x) requires(not constant) : matpack_view(&x, constant_array<N, 1>()) {}

  //! Constant view of x
  explicit constexpr matpack_view(const T& x) : matpack_view(const_cast<T*>(&x), constant_array<N, 1>()) {}

  constexpr matpack_view& operator=(const me_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const ce_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const ms_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const cs_view& x) requires (constant) = delete;

  constexpr matpack_view& operator=(const matpack_data<T, N>& x) requires(not constant) {
    *this = x.view;
    return *this;
  }
  constexpr matpack_view& operator=(const me_view& x) requires (not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::copy(x.elem_begin(), x.elem_end(), elem_begin());
    return *this;
  }
  constexpr matpack_view& operator=(const ce_view& x) requires (not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::copy(x.elem_begin(), x.elem_end(), elem_begin());
    return *this;
  }
  constexpr matpack_view& operator=(const ms_view& x) requires (not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::copy(x.elem_begin(), x.elem_end(), elem_begin());
    return *this;
  }
  constexpr matpack_view& operator=(const cs_view& x) requires (not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::copy(x.elem_begin(), x.elem_end(), elem_begin());
    return *this;
  }

  constexpr matpack_view(matpack_view&&) noexcept = default;

  constexpr matpack_view& operator=(me_view&& x) requires (constant) = delete;
  constexpr matpack_view& operator=(ce_view&& x) requires (constant) = delete;
  constexpr matpack_view& operator=(ms_view&& x) requires (constant) = delete;
  constexpr matpack_view& operator=(cs_view&& x) requires (constant) = delete;
  constexpr matpack_view& operator=(matpack_data<T, N>&& x) requires (constant) = delete;

  constexpr matpack_view& operator=(me_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(ce_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(ms_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(cs_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(matpack_data<T, N>&& x) requires (not constant) { *this = std::move(x.view); return *this; };

  //! Return the total number of elements that can be accessed in this object
  [[nodiscard]] constexpr auto size() const { return static_cast<Index>(view.size()); }

  //! Return the extent of dimenion i
  [[nodiscard]] constexpr auto extent(Index i) const { return view.extent(i); }

  //! Return the stride of dimenion i
  [[nodiscard]] constexpr auto stride(Index i) const { return view.stride(i); }

  //! Return the shape of the object
  [[nodiscard]] constexpr std::array<Index, N> shape() const {
    std::array<Index, N> out;
    for (Index i=0; i<N; i++) out[i] = extent(i);
    return out;
  }

  //! Return the shape of the object skipping the left-most extent
  [[nodiscard]] constexpr std::array<Index, N> inner_shape() const {
    std::array<Index, N-1> out;
    for (Index i=0; i<N-1; i++) out[i] = extent(i+1);
    return out;
  }

  //! Return the strides of the object
  [[nodiscard]] constexpr std::array<Index, N> strides() const {
    std::array<Index, N> out;
    for (Index i=0; i<N; i++) out[i] = stride(i);
    return out;
  }

  //! Return the strides of the object skipping the left-most stride
  [[nodiscard]] constexpr std::array<Index, N> inner_strides() const {
    std::array<Index, N-1> out;
    for (Index i=0; i<N-1; i++) out[i] = stride(i+1);
    return out;
  }

  //! Return a full mdspan mapping of this object
  [[nodiscard]] constexpr auto inner_map() const {
    if constexpr (N == 1) {
      return inner_shape();
    } else {
      using map_t = typename matpack_view<T, N - 1, constant,
                                          strided>::view_type::mapping_type;

      if constexpr (strided) {
        return map_t{inner_shape(), inner_strides()};
      } else {
        return map_t{inner_shape()};
      }
    }
  }

  //! Return the size of the 1-dimensional object
  [[nodiscard]] constexpr Index nelem() const requires(N == 1) { return view.extent(N - 1); }

  //! Return the right-most extent of the multidimensional object
  [[nodiscard]] constexpr Index ncols() const requires(N >= 2) { return view.extent(N - 1); }

  //! Return the right-most extent bar 1 of the multidimensional object
  [[nodiscard]] constexpr Index nrows() const requires(N >= 2) { return view.extent(N - 2); }

  //! Return the right-most extent bar 2 of the multidimensional object
  [[nodiscard]] constexpr Index npages() const requires(N >= 3) { return view.extent(N - 3); }

  //! Return the right-most extent bar 3 of the multidimensional object
  [[nodiscard]] constexpr Index nbooks() const requires(N >= 4) { return view.extent(N - 4); }

  //! Return the right-most extent bar 4 of the multidimensional object
  [[nodiscard]] constexpr Index nshelves() const requires(N >= 5) { return view.extent(N - 5); }

  //! Return the right-most extent bar 5 of the multidimensional object
  [[nodiscard]] constexpr Index nvitrines() const requires(N >= 6) { return view.extent(N - 6); }

  //! Return the right-most extent bar 6 of the multidimensional object
  [[nodiscard]] constexpr Index nlibraries() const requires(N >= 7) { return view.extent(N - 7); }

  /** Returns if the object is empty
   * 
   * @return true if the object views no data
   * @return false otherwise
   */
  [[nodiscard]] constexpr bool empty() const { return view.empty(); }

  template <access_operator... access,
            Index M = num_index<access...>,
            class ret_t = mutable_access<access...>>
  [[nodiscard]] constexpr auto operator()(access... ind) -> ret_t
    requires(sizeof...(access) == N and not constant)
  {
    ARTS_ASSERT(check_index_sizes(view, 0, ind...),
                "Out-of-bounds: ", shape_help<N>(shape()))
    if constexpr (N == M)
      return view(ind...);
    else if constexpr (has_any_range<access...>)
      return ret_t{stdx::submdspan(view, ind...)}.range_adaptor(
          std::forward<access>(ind)...);
    else
      return stdx::submdspan(view, ind...);
  }

  template <access_operator... access,
            Index M = num_index<access...>,
            class ret_t = constant_access<access...>>
  [[nodiscard]] constexpr auto operator()(access... ind) const -> ret_t
    requires(sizeof...(access) == N)
  {
    ARTS_ASSERT(check_index_sizes(view, 0, ind...),
                "Out-of-bounds: ", shape_help<N>(shape()))
    if constexpr (N == M)
      return view(ind...);
    else if constexpr (has_any_range<access...>)
      return ret_t{stdx::submdspan(view, ind...)}.range_adaptor(
          std::forward<access>(ind)...);
    else
      return stdx::submdspan(view, ind...);
  }

  template <access_operator access, Index M = num_index<access>,
            class ret_t = mutable_left_access<access>>
  [[nodiscard]] constexpr auto operator[](access ind) -> ret_t
    requires(not constant)
  {
    ARTS_ASSERT(check_index_sizes(view, 0, ind),
                "Out-of-bounds: ", shape_help<N>(shape()))
    if constexpr (N == M and N == 1)
      return view[ind];
    else if constexpr (has_any_range<access>) {
      ret_t out{view};
      out.template strided_dimension<0>(ind);
      return out;
    } else
      return sub<0, N>(view, ind);
  }

  template <access_operator access, Index M = num_index<access>,
            class ret_t = constant_left_access<access>>
  [[nodiscard]] constexpr auto operator[](access ind) const -> ret_t {
    ARTS_ASSERT(check_index_sizes(view, 0, ind),
                "Out-of-bounds: ", shape_help<N>(shape()))
    if constexpr (N == M and N == 1)
      return view[ind];
    else if constexpr (has_any_range<access>) {
      ret_t out{view};
      out.template strided_dimension<0>(ind);
      return out;
    } else
      return sub<0, N>(view, ind);
  }

  //! The value type of this matpack data is public information
  using value_type = T;

  //! Return the rank of this object
  [[nodiscard]] static constexpr auto rank() { return N; }

  //! Return the constness of this object
  [[nodiscard]] static constexpr auto is_const() { return constant; }

  //! Return a pointer to this object's allocated data that is guaranteed to be continuous until size()
  [[nodiscard]] constexpr auto data_handle() const requires(not strided) { return view.data_handle(); }

  //! Return a pointer to this object's allocated data with no guarantee about continuousness
  [[nodiscard]] constexpr auto unsafe_data_handle() const { return view.data_handle(); }

  using iterator = matpack_mditer<0, false, matpack_view, matpack_view<T, std::max<Index>(N - 1, 1), false, strided>>;
  using const_iterator = matpack_mditer<0, true, const matpack_view, const matpack_view<T, std::max<Index>(N - 1, 1), true, strided>>;

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto begin() requires (not constant) {if constexpr (not strided and N == 1) return data_handle(); else return iterator{*this};}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto end() requires (not constant) {if constexpr (not strided and N == 1) return data_handle()+size(); else return iterator{*this} + extent(0);}

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto begin() const {if constexpr (not strided and N == 1) return data_handle(); else return const_iterator{*this};}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto end() const {if constexpr (not strided and N == 1) return data_handle()+size(); else return const_iterator{*this} + extent(0);}

  //! Iterate over this object by left-most dimension --- return the first of these iterators
  [[nodiscard]] constexpr auto cbegin() const {return begin();}

  //! Iterate over this object by left-most dimension --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto cend() const {return end();}

private:
  //! Helper function to compute the position of an element in the mdspan
  [[nodiscard]] static constexpr std::array<Index, N> mdpos(const std::array<Index, N>& sh, Index ind) {
    std::array<Index, N> pos;
    for (Index i=N-1; i>=0; i--) {
      pos[i] = ind % sh[i];
      ind /= sh[i];
    }
    return pos;
  }

  //! Helper function to take position and stridedness and compute the index it represents
  [[nodiscard]] static constexpr Index mdind(const std::array<Index, N>& pos, const std::array<Index, N>& str) {
      return std::transform_reduce(pos.begin(), pos.end(), str.begin(),
                                   Index{0}, std::plus<>(),
                                   std::multiplies<>());
  }

  //! Finds the true index even if the object is strided
  [[nodiscard]] constexpr Index mdind(Index ind) const {
    if constexpr (view_type::is_always_exhaustive()) {
      return ind;
    } else {
      if constexpr (N == 1) return stride(0) * ind;
      else return mdind(mdpos(shape(), ind), strides());
    }
  }

public:
  //! Access the data elements regardless of the matpack data rank
  [[nodiscard]] constexpr T& elem_at(Index ind) requires(not constant) {
    ARTS_ASSERT(ind < size(), ind, " vs ", size())
    return view.accessor().access(view.data_handle(), mdind(ind));
  }

  //! Access the data elements regardless of the matpack data rank
  [[nodiscard]] constexpr const T& elem_at(Index ind) const {
    ARTS_ASSERT(ind < size(), ind, " vs ", size())
    return view.accessor().access(view.data_handle(), mdind(ind));
  }

  //! Access the data elements regardless of the matpack data rank
  [[nodiscard]] constexpr T& elem_at(const std::array<Index, N>& pos) requires(not constant) {
    ARTS_ASSERT(pos < shape(), shape_help<N>{pos}, " vs ", shape_help<N>{shape()})
    return elem_at(mdind(pos, strides()));
  }

  //! Access the data elements regardless of the matpack data rank
  [[nodiscard]] constexpr const T& elem_at(const std::array<Index, N>& pos) const {
    ARTS_ASSERT(pos < shape(), shape_help<N>{pos}, " vs ", shape_help<N>{shape()})
    return elem_at(mdind(pos, strides()));
  }

  using elem_iterator = matpack_elemwise_mditer<T, false, matpack_view>;
  using const_elem_iterator = matpack_elemwise_mditer<T, true, const matpack_view>;
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_begin() requires (not constant) {if constexpr (not strided) return data_handle(); else return elem_iterator{*this};}

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_end() requires (not constant) {if constexpr (not strided) return data_handle()+size(); else return elem_iterator{*this}+size();}
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_begin() const {if constexpr (not strided) return data_handle(); else return const_elem_iterator{*this};}

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_end() const {if constexpr (not strided) return data_handle()+size(); else return const_elem_iterator{*this}+size();}
  
  //! Iterate over this object element-wise --- return the first of these iterators
  [[nodiscard]] constexpr auto elem_cbegin() const {return elem_begin();}

  //! Iterate over this object element-wise --- return the one-past-the-end of these iterators
  [[nodiscard]] constexpr auto elem_cend() const {return elem_end();}

  friend std::ostream& operator<<(std::ostream& os, const matpack_view& mv) {
    constexpr char extra = N == 1 ? ' ' : '\n';
    bool first = true;
    for (auto&& v: mv) {
      if (not first) os << extra; else first = false;
      os << std::forward<decltype(v)>(v);
    }
    return os;
  }

private:
  //! Helper function to compute the shape and stridedness of a complex view of the data
  constexpr std::pair<std::array<Index, N>, std::array<Index, N>> cmplx_sz() const {
    return {shape(), [s = strides()]() { auto out=s;
              std::transform(out.begin(), out.end(), out.begin(),
                             [](auto &str) { return 2 * str; });
              return out;
            }()};
  }

  //! Helper function to implement real or imaginary view of the complex data
  template <bool impl_constant, Index i>
  constexpr auto cmpl_impl() const {
    using real_t = complex_subtype<T>;
    using strided_t = strided_mdspan<real_t, N>;
    using matpack_t = matpack_view<real_t, N, impl_constant, true>;
    using mapping_type = typename strided_t::mapping_type;
    const auto [sz, st] = cmplx_sz();
    return matpack_t{strided_t{
        &reinterpret_cast<real_t *>(view.data_handle())[i], mapping_type{sz, st}}};
  }

  //! Helper wrapper for real
  template <bool impl_constant> constexpr auto real_impl() const {
    return cmpl_impl<impl_constant, 0>();
  }

  //! Helper wrapper for imag
  template <bool impl_constant> constexpr auto imag_impl() const {
    return cmpl_impl<impl_constant, 1>();
  }

public:
  //! Return a view of the real part of this object
  [[nodiscard]] constexpr auto real() requires(complex_type<T> and not constant) {
    return real_impl<false>();
  }

  //! Return a view of the real part of this object
  [[nodiscard]] constexpr auto real() const requires(complex_type<T>) {
    return real_impl<true>();
  }

  //! Return a view of the imaginary part of this object
  [[nodiscard]] constexpr auto imag() requires(complex_type<T> and not constant) {
    return imag_impl<false>();
  }

  //! Return a view of the imaginary part of this object
  [[nodiscard]] constexpr auto imag() const requires(complex_type<T>) {
    return imag_impl<true>();
  }

  //! Return a view of the diagonal of this object
  [[nodiscard]] matpack_view<T, 1, constant, true> diagonal() const requires(N == 2) {
    ARTS_ASSERT(ncols() == nrows(), ncols(), " vs ", nrows())
    return strided_mdspan<T, 1>{unsafe_data_handle(), {std::array<Index, 1>{nrows()}, std::array<Index, 1>{nrows()+1}}};
  }

  constexpr matpack_view& operator=(std::convertible_to<T> auto x) requires(not constant) {
        if constexpr (N == 1)
      std::fill(elem_begin(), elem_end(), static_cast<T>(x));
    else
      for (auto v : *this)
        v = x;
    return *this;}
  constexpr matpack_view& operator+=(std::convertible_to<T> auto x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(elem_begin(), elem_end(), elem_begin(),
                     [y = static_cast<T>(x)](auto v) { return v + y; });
    else
      for (auto v : *this)
        v += x;
    return *this;
  }
  constexpr matpack_view& operator-=(std::convertible_to<T> auto x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(elem_begin(), elem_end(), elem_begin(),
                     [y = static_cast<T>(x)](auto v) { return v - y; });
    else
      for (auto v : *this)
        v -= x;
    return *this;
  }
  constexpr matpack_view& operator*=(std::convertible_to<T> auto x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(elem_begin(), elem_end(), elem_begin(),
                     [y = static_cast<T>(x)](auto v) { return v * y; });
    else
      for (auto v : *this)
        v *= x;
    return *this;
  }
  constexpr matpack_view& operator/=(std::convertible_to<T> auto x) requires(not constant) {
    if constexpr (N == 1)
      std::transform(elem_begin(), elem_end(), elem_begin(),
                     [y = static_cast<T>(x)](auto v) { return v / y; });
    else
      for (auto v : *this)
        v /= x;
    return *this;
  }

  template <bool c, bool s> constexpr
  matpack_view& operator+=(const matpack_view<T, N, c, s>& x) requires(not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::transform(elem_begin(), elem_end(), x.elem_begin(), elem_begin(), [](auto a, auto b){return a + b;});
    return *this;
  }
  template <bool c, bool s> constexpr
  matpack_view& operator-=(const matpack_view<T, N, c, s>& x) requires(not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::transform(elem_begin(), elem_end(), x.elem_begin(), elem_begin(), [](auto a, auto b){return a - b;});
    return *this;
  }
  template <bool c, bool s> constexpr
  matpack_view& operator*=(const matpack_view<T, N, c, s>& x) requires(not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::transform(elem_begin(), elem_end(), x.elem_begin(), elem_begin(), [](auto a, auto b){return a * b;});
    return *this;
  }
  template <bool c, bool s> constexpr
  matpack_view& operator/=(const matpack_view<T, N, c, s>& x) requires(not constant) {
    ARTS_ASSERT(shape() == x.shape(), shape_help<N>(shape()), " vs ", shape_help<N>(x.shape()))
    std::transform(elem_begin(), elem_end(), x.elem_begin(), elem_begin(), [](auto a, auto b){return a / b;});
    return *this;
  }
  constexpr matpack_view& operator+=(const matpack_data<T, N>& x) requires(not constant) {*this += x.view; return *this;}
  constexpr matpack_view& operator-=(const matpack_data<T, N>& x) requires(not constant) {*this -= x.view; return *this;}
  constexpr matpack_view& operator*=(const matpack_data<T, N>& x) requires(not constant) {*this *= x.view; return *this;}
  constexpr matpack_view& operator/=(const matpack_data<T, N>& x) requires(not constant) {*this /= x.view; return *this;}

  template <bool c, bool s> constexpr bool operator==(const matpack_view<T, N, c, s>& x) const {
    return shape() == x.shape() and std::equal(elem_begin(), elem_end(), x.elem_begin(), x.elem_end());
  }
  template <bool c, bool s> constexpr bool operator!=(const matpack_view<T, N, c, s>& x) const {
    return not (*this == x);
  }
  constexpr bool operator==(const matpack_data<T, N>& x) const {return *this == x.view; }
  constexpr bool operator!=(const matpack_data<T, N>& x) const {return *this != x.view; }

  //! Helper to view a standard array of data as a 1D object
  template <std::size_t M>
  constexpr matpack_view(const std::array<T, M> &x)
    requires(constant and N == 1)
      : matpack_view(const_cast<T *>(x.data()), {static_cast<Index>(M)}) {
    // The cast here should be OK since you are not able to change any values
    // later on as the type is constant
  }

  //! Helper to view a standard vector of data as a 1D object
  constexpr matpack_view(const std::vector<T> &x)
    requires(constant and N == 1)
      : matpack_view(const_cast<T *>(x.data()),
                     {static_cast<Index>(x.size())}) {
    // The cast here should be OK since you are not able to change any values
    // later on as the type is constant
  }

  //! Allow assigning all the values of a compatible object to this view
  template <matpack_convertible<T, N> U>
  constexpr matpack_view& operator=(const U& x) requires(not constant) {
    const auto ext_sh = mdshape(x);
    ARTS_ASSERT(shape() == ext_sh, shape_help<N>(shape()), " vs ", shape_help<N>(ext_sh))
    for (Index i=0; i<size(); i++) elem_at(i) = mdvalue(x, mdpos(ext_sh, i));
    return *this;
  }

  //! Swaps the data around
  constexpr void swap(matpack_view& other) noexcept requires(not constant) { std::swap(view, other.view); }

  //! Unconditionally sets this object's view to another view --- use carefully
  constexpr void set(const matpack_view& other) { view = other.view; }
};

/** Describe the matpack data type
 * 
 * @param m Any matpack_data type
 * @return std::string of the description
 */
template <typename T, Index N, bool constant, bool strided>
std::string describe(const matpack_view<T, N, constant, strided>& m) {
  using namespace matpack;
  if constexpr (constant and strided)
    return var_string("constant and strided matpack_view of rank ", N, " of shape ", shape_help<N>(m.shape()));
  else if constexpr (constant)
    return var_string("constant and exhaustive matpack_view of rank ", N, " of shape ", shape_help<N>(m.shape()));
  else if constexpr (strided)
    return var_string("strided matpack_view of rank ", N, " of shape ", shape_help<N>(m.shape()));
  else
    return var_string("exhaustive matpack_view of rank ", N, " of shape ", shape_help<N>(m.shape()));
}
}  // namespace matpack

//! A ranged access operator
using Range = matpack::matpack_strided_access;

//! A mutable strided view of a vector of Numeric
using VectorView = matpack::matpack_view<Numeric, 1, false, true>;

//! A mutable strided view of a matrix of Numeric
using MatrixView = matpack::matpack_view<Numeric, 2, false, true>;

//! A mutable strided view of a tensor of Numeric of rank 3
using Tensor3View = matpack::matpack_view<Numeric, 3, false, true>;

//! A mutable strided view of a tensor of Numeric of rank 4
using Tensor4View = matpack::matpack_view<Numeric, 4, false, true>;

//! A mutable strided view of a tensor of Numeric of rank 5
using Tensor5View = matpack::matpack_view<Numeric, 5, false, true>;

//! A mutable strided view of a tensor of Numeric of rank 6
using Tensor6View = matpack::matpack_view<Numeric, 6, false, true>;

//! A mutable strided view of a tensor of Numeric of rank 7
using Tensor7View = matpack::matpack_view<Numeric, 7, false, true>;

//! A constant strided view of a vector of Numeric
using ConstVectorView = matpack::matpack_view<Numeric, 1, true, true>;

//! A constant strided view of a matrix of Numeric
using ConstMatrixView = matpack::matpack_view<Numeric, 2, true, true>;

//! A constant strided view of a tensor of Numeric of rank 3
using ConstTensor3View = matpack::matpack_view<Numeric, 3, true, true>;

//! A constant strided view of a tensor of Numeric of rank 4
using ConstTensor4View = matpack::matpack_view<Numeric, 4, true, true>;

//! A constant strided view of a tensor of Numeric of rank 5
using ConstTensor5View = matpack::matpack_view<Numeric, 5, true, true>;

//! A constant strided view of a tensor of Numeric of rank 6
using ConstTensor6View = matpack::matpack_view<Numeric, 6, true, true>;

//! A constant strided view of a tensor of Numeric of rank 7
using ConstTensor7View = matpack::matpack_view<Numeric, 7, true, true>;

//! A mutable continuous view of a vector of Numeric
using ExhaustiveVectorView = matpack::matpack_view<Numeric, 1, false, false>;

//! A mutable continuous view of a matrix of Numeric
using ExhaustiveMatrixView = matpack::matpack_view<Numeric, 2, false, false>;

//! A mutable continuous view of a tensor of Numeric of rank 3
using ExhaustiveTensor3View = matpack::matpack_view<Numeric, 3, false, false>;

//! A mutable continuous view of a tensor of Numeric of rank 4
using ExhaustiveTensor4View = matpack::matpack_view<Numeric, 4, false, false>;

//! A mutable continuous view of a tensor of Numeric of rank 5
using ExhaustiveTensor5View = matpack::matpack_view<Numeric, 5, false, false>;

//! A mutable continuous view of a tensor of Numeric of rank 6
using ExhaustiveTensor6View = matpack::matpack_view<Numeric, 6, false, false>;

//! A mutable continuous view of a tensor of Numeric of rank 7
using ExhaustiveTensor7View = matpack::matpack_view<Numeric, 7, false, false>;

//! A constant continuous view of a vector of Numeric
using ExhaustiveConstVectorView = matpack::matpack_view<Numeric, 1, true, false>;

//! A constant continuous view of a matrix of Numeric
using ExhaustiveConstMatrixView = matpack::matpack_view<Numeric, 2, true, false>;

//! A constant continuous view of a tensor of Numeric of rank 3
using ExhaustiveConstTensor3View = matpack::matpack_view<Numeric, 3, true, false>;

//! A constant continuous view of a tensor of Numeric of rank 4
using ExhaustiveConstTensor4View = matpack::matpack_view<Numeric, 4, true, false>;

//! A constant continuous view of a tensor of Numeric of rank 5
using ExhaustiveConstTensor5View = matpack::matpack_view<Numeric, 5, true, false>;

//! A constant continuous view of a tensor of Numeric of rank 6
using ExhaustiveConstTensor6View = matpack::matpack_view<Numeric, 6, true, false>;

//! A constant continuous view of a tensor of Numeric of rank 7
using ExhaustiveConstTensor7View = matpack::matpack_view<Numeric, 7, true, false>;

//! A mutable strided view of a vector of Complex
using ComplexVectorView = matpack::matpack_view<Complex, 1, false, true>;

//! A mutable strided view of a matrix of Complex
using ComplexMatrixView = matpack::matpack_view<Complex, 2, false, true>;

//! A constant strided view of a vector of Complex
using ConstComplexVectorView = matpack::matpack_view<Complex, 1, true, true>;

//! A constant strided view of a matrix of Complex
using ConstComplexMatrixView = matpack::matpack_view<Complex, 2, true, true>;

//! A mutable continuous view of a vector of Complex
using ExhaustiveComplexVectorView = matpack::matpack_view<Complex, 1, false, true>;

//! A mutable continuous view of a matrix of Complex
using ExhaustiveComplexMatrixView = matpack::matpack_view<Complex, 2, false, true>;

//! A mutable continuous view of a tensor of Complex of rank 3
using ExhaustiveComplexTensor3View = matpack::matpack_view<Complex, 3, false, true>;

//! A constant continuous view of a vector of Complex
using ExhaustiveConstComplexVectorView = matpack::matpack_view<Complex, 1, true, true>;

//! A constant continuous view of a matrix of Complex
using ExhaustiveConstComplexMatrixView = matpack::matpack_view<Complex, 2, true, true>;

//! A constant continuous view of a tensor of Complex of rank 3
using ExhaustiveConstComplexTensor3View = matpack::matpack_view<Complex, 3, true, true>;
