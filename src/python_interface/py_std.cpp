
#include <py_auto_interface.h>

#include <filesystem>

#include "py_macros.h"

namespace Python {
void py_std(py::module_& m) {
  py::class_<std::filesystem::path>(m, "Path")
      .def(py::init<std::string>())
      .doc() = "Wrapper around std::filesystem::path";
  py::implicitly_convertible<std::string, std::filesystem::path>();
}
}  // namespace Python