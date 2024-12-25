#pragma once

#include <configtypes.h>

#include <experimental/mdspan>
#include <ranges>

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
}  // namespace matpack
