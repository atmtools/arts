#include "python_interface.h"

namespace Python {
void py_basic(py::module_&);
void py_matpack(py::module_&);
void py_ppath(py::module_&);
void py_griddedfield(py::module_&);
void py_time(py::module_&);
void py_tessem(py::module_&);
void py_quantum(py::module_&);
void py_rte(py::module_& m);
void py_telsem(py::module_& m);
void py_species(py::module_& m);
void py_sparse(py::module_& m);
void py_mcantenna(py::module_& m);
void py_scattering(py::module_& m);
void py_spectroscopy(py::module_& m);
void py_jac(py::module_& m);
void py_workspace(py::module_& m);
void py_agenda(py::module_& m);
void py_std(py::module_& m);
void py_global(py::module_& m);

/** Construct a new pybind11 module object to hold all the Arts types and functions
 * 
 * Note: the order of execution mostly does not matter bar for two important things:
 *
 * 1) The auto-generated documentation must know about a type to give the python name
 *
 * 2) The workspace auto-generation should be last, it contains some automatic trans-
 *    lations that would otherwise mess things up
 */
PYBIND11_MODULE(pyarts_cpp, m) {
  m.doc() = "Contains direct C++ interface for Arts";

  auto classes = m.def_submodule(
      "classes",
      "Contains all exposed internal and external Arts classes except Index and Numeric");
  py_std(classes);
  py_basic(classes);
  py_matpack(classes);
  py_griddedfield(classes);
  py_time(classes);
  py_species(classes);
  py_spectroscopy(classes);
  py_ppath(classes);
  py_tessem(classes);
  py_rte(classes);
  py_telsem(classes);
  py_sparse(classes);
  py_mcantenna(classes);
  py_scattering(classes);
  py_jac(classes);
  py_quantum(classes);

  auto workspace = 
      m.def_submodule("workspace",
                      "Contains a way to interactively use the Arts workspace");
  py_agenda(workspace);
  py_global(workspace);

  py_workspace(workspace);  // Must be last, it contains automatic conversion operations
}
}  // namespace Python
