#ifndef py_macros_h
#define py_macros_h

#include <python_interface.h>
#include <xml_io.h>

#include <functional>

namespace Python {
namespace py = nanobind;
}  // namespace Python

constexpr Index negative_clamp(const Index i, const Index n) noexcept {
  return (i < 0) ? i + n : i;
}

#define PythonInterfaceIndexItemAccess(Type)                                   \
  def(                                                                         \
      "__len__", [](const Type& x) { return x.size(); }, "Number of elements") \
      .def(                                                                    \
          "__getitem__",                                                       \
          [](Type& x, Index i)                                                 \
              -> std::shared_ptr<std::remove_cvref_t<decltype(x[i])>> {        \
            i = negative_clamp(i, x.size());                                   \
            if (x.size() <= static_cast<Size>(i) or i < 0)                     \
              throw std::out_of_range(                                         \
                  std::format("Bad index access: {}"                           \
                                                                               \
                              " in object of range [0, {})",                   \
                              i,                                               \
                              x.size()));                                      \
            return std::shared_ptr<std::remove_cvref_t<decltype(x[i])>>(       \
                &x[i], [](void*) {});                                          \
          },                                                                   \
          py::rv_policy::reference_internal,                                   \
          py::keep_alive<0, 1>(),                                              \
          "i"_a,                                                               \
          "Get an item of the list")                                           \
      .def(                                                                    \
          "__setitem__",                                                       \
          [](Type& x, Index i, decltype(x[i]) y) {                             \
            i = negative_clamp(i, x.size());                                   \
            if (x.size() <= static_cast<Size>(i) or i < 0)                     \
              throw std::out_of_range(                                         \
                  std::format("Bad index access: {}",                          \
                              " in object of range [0, {})",                   \
                              i,                                               \
                              x.size()));                                      \
            x[i] = std::move(y);                                               \
          },                                                                   \
          "i"_a,                                                               \
          "value"_a,                                                           \
          "Set an item of the list")

#define PythonInterfaceBasicRepresentation(Type) \
  def("__str__", [](const Type& x) {             \
    return std::format("{}", x);                 \
  }).def("__repr__", [](const Type& x) { return std::format("{}", x); })

#define PythonInterfaceCopyValue(Type)  \
  def("__copy__", [](Type& t) -> Type { \
    return t;                           \
  }).def("__deepcopy__", [](Type& t, py::dict&) -> Type { return t; })

#endif
