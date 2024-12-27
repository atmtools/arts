#pragma once

#include <configtypes.h>

#include <experimental/mdspan>
#include <ranges>

#include "matpack_mdspan_common_sizes.h"

namespace stdx = std::experimental;
namespace stdr = std::ranges;
namespace stdv = std::ranges::views;

namespace matpack {
//! A standard type holding a strided multidimensional array
template <class T, Size N>
using mdstrided_t =
    stdx::mdspan<T, stdx::dextents<Index, N>, stdx::layout_stride>;

//! A standard type holding a multidimensional array
template <class T, Size N>
using mdview_t = stdx::mdspan<T, stdx::dextents<Index, N>>;

//! Our data holder
template <typename T, Size N>
class data_t;

//! Our view holder
template <typename T, Size N>
struct view_t;

//! Our strided view holder
template <typename T, Size N>
struct strided_view_t;

//! Our constant-sized data holder
template <typename T, Size... N>
struct cdata_t;

template <typename T>
concept has_value_type =
    requires { typename std::remove_cvref_t<T>::value_type; };

template <has_value_type T>
using value_type = typename std::remove_cvref_t<T>::value_type;

template <typename T>
concept has_is_const = requires { std::remove_cvref_t<T>::is_const; };

template <has_is_const T>
constexpr bool mdmutable = not std::remove_cvref_t<T>::is_const;

template <typename T>
concept has_element_type =
    requires { typename std::remove_cvref_t<T>::element_type; };

template <has_element_type T>
using element_type = typename std::remove_cvref_t<T>::element_type;

template <typename T>
concept has_mapping_type =
    requires { typename std::remove_cvref_t<T>::mapping_type; };

template <has_element_type T>
using mapping_type = typename std::remove_cvref_t<T>::mapping_type;

// Use forward declarations of concepts
template <typename T, typename U, Size N>
static constexpr bool is_data =
    std::same_as<std::remove_cvref_t<T>, data_t<U, N>>;

template <class T, Size N>
concept ranked_data =
    is_data<T, typename std::remove_cvref_t<T>::element_type, N>;

template <class T>
concept any_data = ranked_data<std::remove_cvref_t<T>, rank<T>()>;

template <typename T, typename U, Size N>
static constexpr bool is_strided_view =
    std::same_as<std::remove_cvref_t<T>, strided_view_t<U, N>>;

template <class T, Size N>
concept ranked_strided_view =
    is_strided_view<T, typename std::remove_cvref_t<T>::element_type, N>;

template <class T>
concept any_strided_view = ranked_strided_view<T, rank<T>()>;

template <typename T, typename U, Size N>
static constexpr bool is_view =
    std::same_as<std::remove_cvref_t<T>, view_t<U, N>>;

template <class T, Size N>
concept ranked_view =
    is_view<T, typename std::remove_cvref_t<T>::element_type, N>;

template <class T>
concept any_view = ranked_view<T, rank<T>()>;

template <class T>
concept any_cdata = std::remove_cvref_t<T>::magic_cdata_t_v;

template <typename T, typename U, Index N>
static constexpr bool is_cdata =
    any_cdata<T> and std::same_as<typename T::element_type, U> and
    rank<T>() == N;

template <class T, Index N>
concept ranked_cdata = is_cdata<std::remove_cvref_t<T>, value_type<T>, N>;
template <typename T, Size N>
concept ranked_md = ranked_data<T, N> or ranked_strided_view<T, N> or
                    ranked_view<T, N> or (any_cdata<T> and rank<T>() == N);

template <typename T>
concept any_md =
    any_cdata<T> or any_data<T> or any_strided_view<T> or any_view<T>;

template <typename T, typename U>
concept typed_md = any_md<T> and std::same_as<value_type<T>, U>;

template <typename T, typename U, Size N>
concept exact_md = ranked_md<T, N> and typed_md<T, U>;

template <typename T, typename U, Size N>
concept ranked_convertible_md =
    ranked_md<T, N> and std::convertible_to<value_type<T>, U>;

template <typename T, Size N>
concept mut_ranked_md = ranked_md<T, N> and mdmutable<T>;

template <typename T>
concept mut_any_md = any_md<T> and mdmutable<T>;

template <typename T, typename U>
concept mut_typed_md = typed_md<T, U> and mdmutable<T>;

template <typename T, typename U, Size N>
concept mut_exact_md = exact_md<T, U, N> and mdmutable<T>;

template <typename T, Size N>
concept const_ranked_md = ranked_md<T, N> and not mdmutable<T>;

template <typename T>
concept const_any_md = any_md<T> and not mdmutable<T>;

template <typename T, typename U>
concept const_typed_md = typed_md<T, U> and not mdmutable<T>;

template <typename T, typename U, Size N>
concept const_exact_md = exact_md<T, U, N> and not mdmutable<T>;
}  // namespace matpack
