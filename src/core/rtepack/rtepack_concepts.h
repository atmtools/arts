#pragma once

#include <matpack.h>

namespace rtepack {
using vec7 = matpack::cdata_t<Numeric, 7>;
using vec4 = matpack::cdata_t<Numeric, 4>;
using mat44 = matpack::cdata_t<Numeric, 4, 4>;
using cmat44 = matpack::cdata_t<Complex, 4, 4>;

struct muelmat;

template <typename T>
concept muelmat_convertible =
    matpack::column_keeper<T> and matpack::row_keeper<T> and
    matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

struct propmat;

template <typename T>
concept propmat_convertible =
    matpack::column_keeper<T> and matpack::row_keeper<T> and
    matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

template <typename T>
concept stokvec_convertible =
    matpack::column_keeper<T> and matpack::rank<T>() == 1 and
    matpack::mdvalue_type_compatible<T, Numeric>;
}  // namespace rtepack
