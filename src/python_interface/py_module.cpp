#include "python_interface.h"

namespace Python {
void py_basic(py::module_ &);
void py_matpack(py::module_ &);
void py_ppath(py::module_ &);
void py_griddedfield(py::module_ &);
void py_time(py::module_ &);
void py_tessem(py::module_&);
void py_rte(py::module_& m);
void py_telsem(py::module_& m);
void py_species(py::module_& m);
void py_workspace(py::module_& m);

PYBIND11_MODULE(pyarts_cpp, m) {
    m.doc() = "Contains direct C++ interface for Arts";

    auto classes = m.def_submodule("classes", "Contains all exposed internal and external classes except Index and Numeric");
    py_basic(classes);
    py_matpack(classes);
    py_ppath(classes);
    py_griddedfield(classes);
    py_time(classes);
    py_tessem(classes);
    py_rte(classes);
    py_telsem(classes);
    py_species(classes);
    py_workspace(classes);
}
}  // namespace Python
