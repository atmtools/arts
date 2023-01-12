#pragma once

#include "matpack_iter.h"

#include <algorithm>
#include <array>
#include <debug.h>
#include <functional>
#include <numeric>
#include <ostream>
#include <string>
#include <type_traits>

namespace matpack {
//! A strided accessor to get the range of a type
struct matpack_strided_access {
  Index offset{0};
  Index extent{-1};
  Index stride{1};
  constexpr matpack_strided_access(Index i0=0, Index len=-1, Index d=1) noexcept : offset(i0), extent(len), stride(d) {}
  constexpr matpack_strided_access(Index i0, Joker, Index d=1) noexcept : offset(i0), stride(d) {}
  constexpr matpack_strided_access(Joker, Index d=1) noexcept : stride(d) {}
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

  constexpr operator Joker() const noexcept {return joker;}
  constexpr operator Joker() noexcept {return joker;}
};

template <access_operator Access, typename... arguments>
[[nodiscard]] consteval Index ranged_access_operator_count() noexcept {
  constexpr bool index = integral<Access>;
  if constexpr (sizeof...(arguments)) {
    return index + ranged_access_operator_count<arguments...>();
  } else {
    return index;
  }
}

template <access_operator... access>
inline constexpr Index num_access_operators = ranged_access_operator_count<access...>();

template <access_operator first, access_operator... access, Index N=sizeof...(access)>
consteval bool any_range() {
  if constexpr (is_matpack_strided_access<first>) return true;
  else if constexpr (N not_eq 0) return any_range<access...>();
  else return false;
};

template <access_operator... access>
inline constexpr bool has_any_range = any_range<access...>();

template <access_operator T, access_operator ... Ts>
[[nodiscard]] constexpr bool check_index_sizes(auto&& arr, Index r, T first, Ts... rest) {
  constexpr bool test = integral<T>;
  if constexpr (test) {
    if (first >= arr.extent(r)) return false;
  }

  if constexpr (sizeof...(Ts) == 0) {
    return true;
  } else {
    return check_index_sizes(arr, r+1, rest...);
  }
}

namespace access_retval {
template <typename T, Index N, bool constant, bool strided,
          access_operator... access>
struct helper {
  static constexpr Index M = num_access_operators<access...>;
  static constexpr bool ranged = has_any_range<access...>;

  using mutable_type =
      std::conditional_t<M == N, T &,
                         matpack_view<T, N - M, constant, ranged or strided>>;

  using constant_type =
      std::conditional_t<M == N, T,
                         matpack_view<T, N - M, true, ranged or strided>>;
};

template <typename T, Index N, bool constant, bool strided, access_operator access>
struct subhelper {
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

template <typename T, Index N, bool constant, bool strided, access_operator... access>
using mutable_access_type = typename access_retval::helper<T, N, constant, strided, access...>::mutable_type;

template <typename T, Index N, bool constant, bool strided, access_operator... access>
using constant_access_type = typename access_retval::helper<T, N, constant, strided, access...>::constant_type;

template <typename T, Index N, bool constant, bool strided, access_operator access>
using mutable_left_access_type = typename access_retval::subhelper<T, N, constant, strided, access>::mutable_type;

template <typename T, Index N, bool constant, bool strided, access_operator access>
using constant_left_access_type = typename access_retval::subhelper<T, N, constant, strided, access>::constant_type;

template <Index N> struct shape_help {
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

template <Index N, Index v> 
constexpr std::array<Index, N> constant_array() {
  std::array<Index, N> x;
  x.fill(v);
  return x;
}

template <Index N, Index M>
std::array<Index, N> upview(const std::array<Index, M> &low) {
  std::array<Index, N> sz = constant_array<N, 1>();
  for (Index i = 0; i < M; i++)
  sz[i] = low[i];
  return sz;
}

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

  template <typename U, Index M> friend class matpack_data;
  template <typename U, Index M, bool c, bool s> friend class matpack_view;
  friend class matpack_elemwise_mditer<T, false, matpack_view>;
  friend class matpack_elemwise_mditer<T, true, const matpack_view>;

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
  constexpr matpack_view(exhaustive_view x) : view(std::move(x)) {}
  constexpr matpack_view(strided_view x) : view(std::move(x)) {}
  constexpr matpack_view(T* ptr, std::array<Index, N> sz) : matpack_view(exhaustive_view{ptr, sz}) {}

  constexpr matpack_view(const me_view& x) : view(x.view) {}
  constexpr matpack_view(const ce_view& x) requires(constant) : view(x.view) {}
  constexpr matpack_view(const ms_view& x) requires(strided) : view(x.view) {}
  constexpr matpack_view(const cs_view& x) requires(constant and strided) : view(x.view) {}
  constexpr matpack_view(const matpack_data<T, N>& x) noexcept : matpack_view(x.view) {}

  template <Index M> explicit constexpr matpack_view(const matpack_data<T, M>& x) requires(N > M) : matpack_view(x.data_handle(), upview<N, M>(x.shape())) {}
  template <Index M> explicit matpack_view(const matpack_view<T, M, constant, strided>& x) requires(N > M) {
    using map_t = typename view_type::mapping_type;
    if constexpr (not strided) {
      view = view_type{x.data_handle(), map_t{upview<N, M>(x.shape())}};
    } else {
      view = view_type{x.unsafe_data_handle(), map_t{upview<N, M>(x.shape()), upview<N, M>(x.strides())}};
    }
  }
  explicit constexpr matpack_view(T& x) requires(not constant) : matpack_view(&x, constant_array<N, 1>()) {}
  explicit constexpr matpack_view(const T& x) : matpack_view(const_cast<T*>(&x), constant_array<N, 1>()) {}

  constexpr matpack_view& operator=(const me_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const ce_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const ms_view& x) requires (constant) = delete;
  constexpr matpack_view& operator=(const cs_view& x) requires (constant) = delete;

  constexpr matpack_view& operator=(const matpack_data<T, N>& x) {
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

  constexpr matpack_view& operator=(me_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(ce_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(ms_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(cs_view&& x) requires (not constant) {*this = x; return *this;};
  constexpr matpack_view& operator=(matpack_data<T, N>&& x) = delete;

  [[nodiscard]] constexpr auto size() const { return static_cast<Index>(view.size()); }
  [[nodiscard]] constexpr auto extent(Index i) const { return view.extent(i); }
  [[nodiscard]] constexpr auto stride(Index i) const { return view.stride(i); }
  [[nodiscard]] constexpr std::array<Index, N> shape() const {
    std::array<Index, N> out;
    for (Index i=0; i<N; i++) out[i] = extent(i);
    return out;
  }
  [[nodiscard]] constexpr std::array<Index, N> inner_shape() const {
    std::array<Index, N-1> out;
    for (Index i=0; i<N-1; i++) out[i] = extent(i+1);
    return out;
  }
  [[nodiscard]] constexpr std::array<Index, N> strides() const {
    std::array<Index, N> out;
    for (Index i=0; i<N; i++) out[i] = stride(i);
    return out;
  }
  [[nodiscard]] constexpr std::array<Index, N> inner_strides() const {
    std::array<Index, N-1> out;
    for (Index i=0; i<N-1; i++) out[i] = stride(i+1);
    return out;
  }
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
  [[nodiscard]] constexpr Index nelem() const requires(N == 1) { return view.extent(N - 1); }
  [[nodiscard]] constexpr Index ncols() const requires(N >= 1) { return view.extent(N - 1); }
  [[nodiscard]] constexpr Index nrows() const requires(N >= 2) { return view.extent(N - 2); }
  [[nodiscard]] constexpr Index npages() const requires(N >= 3) { return view.extent(N - 3); }
  [[nodiscard]] constexpr Index nbooks() const requires(N >= 4) { return view.extent(N - 4); }
  [[nodiscard]] constexpr Index nshelves() const requires(N >= 5) { return view.extent(N - 5); }
  [[nodiscard]] constexpr Index nvitrines() const requires(N >= 6) { return view.extent(N - 6); }
  [[nodiscard]] constexpr Index nlibraries() const requires(N >= 7) { return view.extent(N - 7); }
  [[nodiscard]] constexpr bool empty() const { return view.empty(); }

  template <access_operator... access,
            Index M = num_access_operators<access...>,
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
            Index M = num_access_operators<access...>,
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

  template <access_operator access, Index M = num_access_operators<access>,
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

  template <access_operator access, Index M = num_access_operators<access>,
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

  using value_type = T;
  [[nodiscard]] static constexpr auto rank() { return N; }
  [[nodiscard]] static constexpr auto is_const() { return constant; }
  [[nodiscard]] constexpr auto data_handle() const requires(not strided) { return view.data_handle(); }
  [[nodiscard]] constexpr auto unsafe_data_handle() const { return view.data_handle(); }

  using iterator = matpack_mditer<0, false, matpack_view, matpack_view<T, std::max<Index>(N - 1, 1), false, strided>>;
  using const_iterator = matpack_mditer<0, true, const matpack_view, const matpack_view<T, std::max<Index>(N - 1, 1), true, strided>>;
  [[nodiscard]] constexpr auto begin() requires (not constant) {if constexpr (not strided and N == 1) return data_handle(); else return iterator{*this};}
  [[nodiscard]] constexpr auto end() requires (not constant) {if constexpr (not strided and N == 1) return data_handle()+size(); else return iterator{*this} + extent(0);}
  [[nodiscard]] constexpr auto begin() const {if constexpr (not strided and N == 1) return data_handle(); else return const_iterator{*this};}
  [[nodiscard]] constexpr auto end() const {if constexpr (not strided and N == 1) return data_handle()+size(); else return const_iterator{*this} + extent(0);}
  [[nodiscard]] constexpr auto cbegin() const {return begin();}
  [[nodiscard]] constexpr auto cend() const {return end();}

private:
  [[nodiscard]] static constexpr std::array<Index, N> mdpos(const std::array<Index, N>& sh, Index ind) {
    std::array<Index, N> pos;
    for (Index i=N-1; i>=0; i--) {
      pos[i] = ind % sh[i];
      ind /= sh[i];
    }
    return pos;
  }

  [[nodiscard]] static constexpr Index mdind(const std::array<Index, N>& pos, const std::array<Index, N>& str) {
      return std::transform_reduce(pos.begin(), pos.end(), str.begin(),
                                   Index{0}, std::plus<>(),
                                   std::multiplies<>());
  }

  [[nodiscard]] constexpr Index mdind(Index ind) const {
    if constexpr (view_type::is_always_exhaustive()) {
      return ind;
    } else {
      return mdind(mdpos(shape(), ind), strides());
    }
  }

public:
  [[nodiscard]] constexpr T& elem_at(Index ind) requires(not constant) {
    ARTS_ASSERT(ind < size(), ind, " vs ", size())
    if constexpr (is_always_exhaustive_v<view_type>) {
      return view.accessor().access(view.data_handle(), ind);
    } else {
      return view.accessor().access(view.data_handle(), mdind(ind));
    }
  }

  [[nodiscard]] constexpr T elem_at(Index ind) const {
    ARTS_ASSERT(ind < size(), ind, " vs ", size())
    if constexpr (is_always_exhaustive_v<view_type>) {
      return view.accessor().access(view.data_handle(), ind);
    } else {
      return view.accessor().access(view.data_handle(), mdind(ind));
    }
  }

  [[nodiscard]] constexpr T& elem_at(const std::array<Index, N>& pos) requires(not constant) {
    ARTS_ASSERT(pos < shape(), shape_help<N>{pos}, " vs ", shape_help<N>{shape()})
    return elem_at(mdind(pos), strides());
  }

  [[nodiscard]] constexpr T elem_at(const std::array<Index, N>& pos) const {
    ARTS_ASSERT(pos < shape(), shape_help<N>{pos}, " vs ", shape_help<N>{shape()})
    return elem_at(mdind(pos), strides());
  }

  using elem_iterator = matpack_elemwise_mditer<T, false, matpack_view>;
  using const_elem_iterator = matpack_elemwise_mditer<T, true, const matpack_view>;
  [[nodiscard]] constexpr auto elem_begin() requires (not constant) {if constexpr (not strided) return data_handle(); else return elem_iterator{*this};}
  [[nodiscard]] constexpr auto elem_end() requires (not constant) {if constexpr (not strided) return data_handle()+size(); else return elem_iterator{*this}+size();}
  [[nodiscard]] constexpr auto elem_begin() const {if constexpr (not strided) return data_handle(); else return const_elem_iterator{*this};}
  [[nodiscard]] constexpr auto elem_end() const {if constexpr (not strided) return data_handle()+size(); else return const_elem_iterator{*this}+size();}
  [[nodiscard]] constexpr auto elem_cbegin() const {return elem_begin();}
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
  constexpr std::pair<std::array<Index, N>, std::array<Index, N>> cmplx_sz() const {
    return {shape(), [s = strides()]() { auto out=s;
              std::transform(out.begin(), out.end(), out.begin(),
                             [](auto &str) { return 2 * str; });
              return out;
            }()};
  }

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

  template <bool impl_constant> constexpr auto real_impl() const {
    return cmpl_impl<impl_constant, 0>();
  }

  template <bool impl_constant> constexpr auto imag_impl() const {
    return cmpl_impl<impl_constant, 1>();
  }

public:
  [[nodiscard]] constexpr auto real() requires(complex_type<T> and not constant) {
    return real_impl<false>();
  }

  [[nodiscard]] constexpr auto real() const requires(complex_type<T>) {
    return real_impl<true>();
  }

  [[nodiscard]] constexpr auto imag() requires(complex_type<T> and not constant) {
    return imag_impl<false>();
  }

  [[nodiscard]] constexpr auto imag() const requires(complex_type<T>) {
    return imag_impl<true>();
  }

  constexpr matpack_view& operator=(T x) requires(not constant) {std::fill(elem_begin(), elem_end(), x); return *this;}
  constexpr matpack_view& operator+=(T x) requires(not constant) {std::transform(elem_begin(), elem_end(), elem_begin(), [x](auto v){return v + x;}); return *this;}
  constexpr matpack_view& operator-=(T x) requires(not constant) {std::transform(elem_begin(), elem_end(), elem_begin(), [x](auto v){return v - x;}); return *this;}
  constexpr matpack_view& operator*=(T x) requires(not constant) {std::transform(elem_begin(), elem_end(), elem_begin(), [x](auto v){return v * x;}); return *this;}
  constexpr matpack_view& operator/=(T x) requires(not constant) {std::transform(elem_begin(), elem_end(), elem_begin(), [x](auto v){return v / x;}); return *this;}

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

  template <matpack_convertible<T, N> U>
  constexpr matpack_view& operator=(const U& x) requires(not constant) {
    const auto ext_sh = mdshape(x);
    ARTS_ASSERT(shape() == ext_sh, shape_help<N>(shape()), " vs ", shape_help<N>(ext_sh))
    for (Index i=0; i<size(); i++) elem_at(i) = mdvalue(x, mdpos(ext_sh, i));
    return *this;
  }

  constexpr void swap(matpack_view& other) noexcept requires(not constant) { std::swap(view, other.view); }
};
}  // namespace matpack

template <typename T, Index N, bool constant, bool strided>
std::string describe(const matpack::matpack_view<T, N, constant, strided>& m) {
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
