#pragma once

#include "python_interface.h"

namespace Python {
namespace {
template <typename FromType, typename ToType>
void implicit_convert_gf(py::class_<ToType>& cls) {
  cls.def("__init__", [](ToType* x, const FromType& v) {
    ToType t1(v);
    new (x) ToType(std::move(t1));
  });
  py::implicitly_convertible<FromType, ToType>();
}
}  // namespace
}  // namespace Python
