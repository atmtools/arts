#pragma once

#include <matpack.h>

namespace rtepack {
struct propmat;
struct muelmat;
struct specmat;
struct stokvec;
}  // namespace rtepack

namespace matpack {
template <>
struct any_cdata_t<rtepack::propmat> : std::true_type {};
template <>
struct any_cdata_t<rtepack::muelmat> : std::true_type {};
template <>
struct any_cdata_t<rtepack::specmat> : std::true_type {};
template <>
struct any_cdata_t<rtepack::stokvec> : std::true_type {};
}  // namespace matpack
