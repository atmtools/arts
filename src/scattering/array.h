#ifndef __ARTS_SCATTERING_ARRAY__
#define __ARTS_SCATTERING_ARRAY__

#include <array>
#include <utility>

namespace scattering {
namespace detail {
template <typename T,
  size_t M,
  size_t N,
  size_t... IndicesM,
  size_t... IndicesN>
  std::array<T, N + M> concat_impl(
        std::array<T, M> head,
        std::array<T, N> tail,
        std::index_sequence<IndicesM...> /*indices_m*/,
        std::index_sequence<IndicesN...> /*indices_n*/) {
    return std::array<T, N + M>{std::get<IndicesM>(head)...,
            std::get<IndicesN>(tail)...};
}

template <size_t start, typename T, size_t N, size_t... Indices>
std::array<T, sizeof...(Indices)> take_impl(std::array<T, N> array, std::index_sequence<Indices ...> /*indices*/) {
  return std::array<T, sizeof...(Indices)>{std::get<start + Indices>(array)...};
}
}

/** Concatenates array to another.
 *
 * @param head The array to which to concatenate the other.
 * @param tail The array which should be concatenated.
 * @return The array resulting from concatenating head and tail.
 */
template <typename T, size_t M, size_t N>
    std::array<T, N + M> concat(std::array<T, M> head,
                                std::array<T, N> tail) {
    return detail::concat_impl(head,
                               tail,
                               std::make_index_sequence<M>(),
                               std::make_index_sequence<N>());
}

template <typename T, size_t M, size_t N, size_t O, typename ... Ts>
    std::array<T, N + M> concat(std::array<T, M> head,
                                std::array<T, N> tail,
                                std::array<T, N> other_tail,
                                Ts ... ts) {
    return concat(concat(head, tail), other_tail, ts ...);
}



template <typename T, typename... Ts>
std::array<T, sizeof...(Ts) + 1> make_array(const T &t, const Ts &... ts) {
  return {t, ts...};
}

/** Take first M elements from array.
 *
 * @tparam M How many element to take.
 * @param in The input array
 * @return Array containing the M first elements of in.
 */
template <size_t M, typename T, size_t N>
std::array<T, M> take_until(std::array<T, N> in) {
  return detail::take_impl<0>(in, std::make_index_sequence<M>());
}

/** Drop first M elements from array.
*
* @tparam M How many element to drop.
* @param in The input array
* @return Array containing the elements left after dropping the
*         its M first elements.
*/
template <size_t M, typename T, size_t N>
std::array<T, N - M> take_from(std::array<T, N> in) {
  return detail::take_impl<M>(in, std::make_index_sequence<N - M>());
}

}  // namespace scattering

#endif
