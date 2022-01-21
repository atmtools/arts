#include "python_interface.h"

namespace Python {
void py_basic(py::module_ &);
void py_matpack(py::module_ &);

PYBIND11_MODULE(pyarts_cpp, m) {
    m.attr("__doc__") = "Contains direct C++ interface for Arts";

    auto classes = m.def_submodule("classes", "Contains internal and external classes");
    py_basic(classes);
    py_matpack(classes);
}
}  // namespace Python
