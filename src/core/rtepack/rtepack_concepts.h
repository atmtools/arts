#pragma once

#include <matpack.h>
#include <matpack_lazy.h>

namespace rtepack {
using namespace matpack::lazy;

using vec7 = matpack::matpack_constant_data<Numeric, 7>;
using vec4 = matpack::matpack_constant_data<Numeric, 4>;
using mat44 = matpack::matpack_constant_data<Numeric, 4, 4>;

struct muelmat;

template <typename T>
concept lazy_muelmat =
    constexpr_smat_data_like<T> and std::remove_cvref_t<T>::size() == 4;

template <typename T>
concept muelmat_convertible =
    matpack::column_keeper<T> and matpack::row_keeper<T> and
    matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

struct propmat;

template <typename T>
concept lazy_propmat =
    constexpr_vec_data_like<T> and std::remove_cvref_t<T>::size() == 7;

template <typename T>
concept propmat_convertible =
    matpack::column_keeper<T> and matpack::row_keeper<T> and
    matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

struct stokvec;
template <typename T>
concept lazy_stokvec =
    constexpr_vec_data_like<T> and std::remove_cvref_t<T>::size() == 4;

template <typename T>
concept stokvec_convertible =
    matpack::column_keeper<T> and matpack::rank<T>() == 1 and
    matpack::mdvalue_type_compatible<T, Numeric>;
} // namespace rtepack
