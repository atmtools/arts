#pragma once

#include <array>

#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_iter.h"

namespace matpack {
namespace detail {
template <Index I, Index N>
void repeat_set(std::array<Vector, N> &arr, const auto &val, const Index i) {
  std::get<I>(arr)[i] = std::get<I>(val);
}

template <Index... inds, std::size_t N = sizeof...(inds)>
void repeat_set(std::integer_sequence<Index, inds...>,
                std::array<Vector, N> &arr, const auto &val, const Index i) {
  (repeat_set<inds, N>(arr, val, i), ...);
}
} // namespace detail

template <strict_rank_matpack_type<1>... vectors,
          std::size_t N = sizeof...(vectors)>
std::array<Vector, N> repeat(vectors &&...vec)
  requires(N > 0)
{
  std::array<Vector, N> out;
  out.fill(Vector((vec.nelem() * ...)));

  Index i = 0;
  for (auto &&x : elemwise{vec...}) {
    detail::repeat_set(std::make_integer_sequence<Index, N>{}, out, x, i);
    i++;
  }

  return out;
}
}  // namespace matpack
