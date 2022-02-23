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
void py_workspace(py::module_& m, py::class_<Workspace>& ws);
void py_agenda(py::module_& m);
void py_std(py::module_& m);
void py_global(py::module_& m);

/** Construct a new pybind11 module object to hold all the Arts types and functions
 * 
 * Note: the order of execution mostly does not matter bar for some important things:
 *
 * 1) The auto-generated documentation must know about a type to give the python name
 *
 * 2) The workspace auto-generation should be last, it contains some automatic trans-
 *    lations that would otherwise mess things up
 * 
 * 3) Implicit conversion can only be defined between two python-defined Arts types
 */
PYBIND11_MODULE(pyarts_cpp, m) {
  m.doc() = "Contains direct C++ interface for Arts";
  auto ws = py::class_<Workspace>(m, "Workspace");

  py_std(m);
  py_basic(m);
  py_matpack(m);
  py_griddedfield(m);
  py_time(m);
  py_species(m);
  py_spectroscopy(m);
  py_ppath(m);
  py_tessem(m);
  py_rte(m);
  py_telsem(m);
  py_sparse(m);
  py_mcantenna(m);
  py_scattering(m);
  py_jac(m);
  py_quantum(m);

  py_agenda(m);
  py_global(m);

  // Must be last, it contains automatic conversion operations
  py_workspace(m, ws);  
}
}  // namespace Python
