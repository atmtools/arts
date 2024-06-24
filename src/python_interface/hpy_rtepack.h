#include "hpy_matpack.h"

namespace Python {
template <typename T>
void common_rtepack_interface(py::class_<T>& c) {
  constexpr auto dim = T::data::size()

  matpack_common_interface(c);
}
}  // namespace Python
