#include <workspace_ng.h>
#include <workspace.h>
#include <memory>

#include "py_auto_interface.h"
#include "python_interface.h"

extern Parameters parameters;

void parse_path_from_environment(String envvar, ArrayOfString& paths);

namespace Python {
namespace py = pybind11;

void py_basic(py::module_&);
void py_docserver(py::module_&);
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
void py_workspace(py::module_& m,
                  py::class_<Workspace, std::shared_ptr<Workspace>>& ws,
                  py::class_<WorkspaceVariable>& wsv);
void py_agenda(py::module_& m);
void py_global(py::module_& m);
void py_xsec(py::module_& m);
void py_nlte(py::module_& m);
void py_constants(py::module_& m);
void py_conversions(py::module_& m);
void py_star(py::module_& m);
void py_physics(py::module_& m);
void py_predefined(py::module_& m);
void py_math(py::module_& m);
void py_options(py::module_& m);

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
PYBIND11_MODULE(arts, m) {
  m.doc() = "Contains direct C++ interface for Arts";
  py::class_<Workspace, std::shared_ptr<Workspace>> ws(m, "Pyarts::Workspace");
  py::class_<WorkspaceVariable> wsv(m, "WorkspaceVariable");
  wsv.doc() = "A wrapper around all workspace variables";

  static bool init = true;
  if (init) {
    init = false;

    define_wsv_groups();
    define_wsv_data();
    define_wsv_map();
    define_md_data_raw();
    expand_md_data_raw_to_md_data();
    define_md_map();
    define_md_raw_map();
    define_agenda_data();
    define_agenda_map();
    ARTS_ASSERT(check_agenda_data());
    global_data::workspace_memory_handler.initialize();

    // Set parameters that are know on first execution
#ifdef ARTS_DEFAULT_INCLUDE_DIR
    String arts_default_include_path(ARTS_DEFAULT_INCLUDE_DIR);
    if (arts_default_include_path != "" && !parameters.includepath.nelem()) {
      // Skip delimiters at beginning.
      String::size_type lastPos =
          arts_default_include_path.find_first_not_of(":", 0);
      // Find first "non-delimiter".
      String::size_type pos =
          arts_default_include_path.find_first_of(":", lastPos);

      while (String::npos != pos || String::npos != lastPos) {
        parameters.includepath.push_back(
            arts_default_include_path.substr(lastPos, pos - lastPos));
        lastPos = arts_default_include_path.find_first_not_of(":", pos);
        pos = arts_default_include_path.find_first_of(":", lastPos);
      }
    }
#endif

    parse_path_from_environment("ARTS_INCLUDE_PATH", parameters.includepath);
    parse_path_from_environment("ARTS_DATA_PATH", parameters.datapath);
    parse_path_from_environment("ARTS_CAT_DATA_DIR", parameters.datapath);
    parse_path_from_environment("ARTS_XML_DATA_DIR", parameters.datapath);

    parameters.includepath.insert(parameters.includepath.begin(), ".");
    parameters.datapath.insert(parameters.datapath.begin(), ".");
  }

  py_basic(m);
  py_docserver(m);
  py_matpack(m);
  py_griddedfield(m);
  py_time(m);
  py_species(m);
  py_quantum(m);
  py_spectroscopy(m);
  py_ppath(m);
  py_tessem(m);
  py_rte(m);
  py_telsem(m);
  py_sparse(m);
  py_mcantenna(m);
  py_scattering(m);
  py_jac(m);
  py_xsec(m);
  py_nlte(m);
  py_predefined(m);
  py_star(m);
  py_agenda(m);

  // Must be last, it contains automatic conversion operations
  py_workspace(m, ws, wsv);

  // Extras calling pure internal functions
  py_constants(m);
  py_conversions(m);
  py_physics(m);
  py_global(m);
  py_physics(m);
  py_math(m);
  py_options(m);
}
}  // namespace Python
