#pragma once

#include "matpack_mdspan_cdata_t.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_math.h"
#include "matpack_mdspan_strided_view_t.h"
#include "matpack_mdspan_view_t.h"

namespace matpack {
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
