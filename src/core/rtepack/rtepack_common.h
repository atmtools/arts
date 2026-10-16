#pragma once

#include <matpack.h>

namespace rtepack {
struct propmat;
struct muelmat;
struct specmat;
struct stokvec;
}  // namespace rtepack

namespace matpack {
template <> struct any_cdata_t<rtepack::propmat> : std::true_type {};
template <> struct any_cdata_t<rtepack::muelmat> : std::true_type {};
template <> struct any_cdata_t<rtepack::specmat> : std::true_type {};
template <> struct any_cdata_t<rtepack::stokvec> : std::true_type {};
}  // namespace matpack

namespace rtepack {
/** The Numeric elements of every fixed-size matrix of x, as a view with one more, trailing, dimension
 *
 * The trailing dimension runs over the elements of each matrix in its storage
 * order: row-major for a muelmat, [F11, F12, F22, F33, F34, F44] for a
 * compact_planar_muelmat.  It is for algorithms that work on one element at a
 * time, such as spherical-harmonic transforms and interpolation, the way
 * matpack's real() and imag() view complex data.
 *
 * @param x A data_t, view_t or strided_view_t of muelmat, compact_planar_muelmat or alike
 */
template <typename X> constexpr auto elements(X &&x) {
  using T           = std::remove_reference_t<decltype(*x.data_handle())>;
  using U           = std::conditional_t<std::is_const_v<T>, const Numeric, Numeric>;
  constexpr Size  N = std::remove_cvref_t<X>::rank();
  constexpr Index K = sizeof(T) / sizeof(Numeric);
  static_assert(K * sizeof(Numeric) == sizeof(T));

  const matpack::strided_view_t<T, N> v{x};
  std::array<Index, N + 1>            shape, strides;
  for (Size i = 0; i < N; i++) {
    shape[i]   = v.extent(i);
    strides[i] = K * v.stride(i);
  }
  shape[N]   = K;
  strides[N] = 1;
  return matpack::strided_view_t<U, N + 1>{
      matpack::mdstrided_t<U, N + 1>{reinterpret_cast<U *>(v.data_handle()), {shape, strides}}};
}
}  // namespace rtepack
