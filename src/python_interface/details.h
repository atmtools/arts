#ifndef python_interface_pyarts_details_h
#define python_interface_pyarts_details_h

#include <py_auto_interface.h>
#include <pybind11/functional.h>

namespace Python::details {
inline static std::function<py::object(py::object&)> one_arg{[](py::object&) {
  throw std::logic_error("Not implemented");
  return py::none();
}};

inline static std::function<py::object(py::object&, py::object&)> two_args{
    [](py::object&, py::object&) {
      throw std::logic_error("Not implemented");
      return py::none();
    }};

inline static std::function<py::object(py::object&, py::object&, py::object&)>
    three_args{[](py::object&, py::object&, py::object&) {
      throw std::logic_error("Not implemented");
      return py::none();
    }};

inline static std::function<py::object(
    py::object&, py::object&, py::object&, py::object&, py::object&)>
    five_args{[](py::object&, py::object&, py::object&, py::object&, py::object&) {
      throw std::logic_error("Not implemented");
      return py::none();
    }};
}  // namespace Python::details

#endif  // python_interface_pyarts_details_h
