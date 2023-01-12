#pragma once

#include <matpack/matpack_view2.h>

template <typename T, bool constant, bool strided>
matpack::matpack_view<T, 2, constant, true>
transpose(const matpack::matpack_view<T, 2, constant, strided> &x) {
  return matpack::strided_mdspan<Numeric, 2>{
      x.unsafe_data_handle(), {std::array{x.extent(1), x.extent(0)},
                               std::array{x.stride(1), x.stride(0)}}};
}
