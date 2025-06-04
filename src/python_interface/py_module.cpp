#include <nanobind/nanobind.h>
#include <parameters.h>
#include <python_interface.h>

extern Parameters parameters;

void parse_path_from_environment(String envvar, ArrayOfString& paths);

namespace Python {
namespace py = nanobind;

void py_workspace(py::class_<Workspace>& ws);

void py_basic(py::module_& m);
void py_matpack(py::module_& m);
void py_path(py::module_& m);
void py_griddedfield(py::module_& m);
void py_time(py::module_& m);
void py_quantum(py::module_& m);
void py_rtepack(py::module_& m);
void py_species(py::module_& m);
void py_sparse(py::module_& m);
void py_nlte(py::module_& m);
void py_scattering(py::module_& m);
void py_scattering_species(py::module_& m);
void py_psd(py::module_& m);
void py_jac(py::module_& m);
void py_agenda(py::module_& m);
void py_global(py::module_& m);
void py_xsec(py::module_& m);
void py_constants(py::module_& m);
void py_conversions(py::module_& m);
void py_star(py::module_& m);
void py_physics(py::module_& m);
void py_predefined(py::module_& m);
void py_math(py::module_& m);
void py_auto_options(py::module_& m);
void py_hitran(py::module_& m);
void py_atm(py::module_& m);
void py_surf(py::module_& m);
void py_fwd(py::module_& m);
void py_cia(py::module_& m);
void py_operators(py::module_& m);
void py_lbl(py::module_& m);
void py_interp(py::module_& m);
void py_sensor(py::module_& m);
void py_disort(py::module_& m);
void py_igrf(py::module_& m);
void py_zeeman(py::module_& m);
void py_retrieval(py::module_& m);
void py_lookup(py::module_& m);
void py_file(py::module_& m);
void py_auto_agenda_operators(py::module_& m);

/** Construct a new nanobind module object to hold all the Arts types and functions
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
NB_MODULE(arts, m) try {
  m.doc() = "Interface directly to the C++ types, functions and modules via python";
  py::class_<Workspace> ws(m, "CxxWorkspace");

  static bool init = true;
  if (init) {
    init = false;

    // Set parameters that are know on first execution
#ifdef ARTS_DEFAULT_INCLUDE_DIR
    String arts_default_include_path(ARTS_DEFAULT_INCLUDE_DIR);
    if (arts_default_include_path != "" && !parameters.includepath.size()) {
      // Skip delimiters at beginning.
      String::size_type lastPos =
          arts_default_include_path.find_first_not_of(':', 0);
      // Find first "non-delimiter".
      String::size_type pos =
          arts_default_include_path.find_first_of(':', lastPos);

      while (String::npos != pos || String::npos != lastPos) {
        parameters.includepath.push_back(
            arts_default_include_path.substr(lastPos, pos - lastPos));
        lastPos = arts_default_include_path.find_first_not_of(':', pos);
        pos     = arts_default_include_path.find_first_of(':', lastPos);
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

  //! The options names space depends only on static c++ data,
  // so it should be included early for documentation purposes when its
  // data is used by modules below it
  py_auto_options(m);

  py_basic(m);
  py_matpack(m);

  py_interp(m);
  py_griddedfield(m);
  py_time(m);
  py_species(m);
  py_quantum(m);
  py_path(m);
  py_rtepack(m);
  py_sparse(m);
  py_scattering(m);
  py_scattering_species(m);
  py_psd(m);
  py_jac(m);
  py_xsec(m);
  py_predefined(m);
  py_star(m);
  py_agenda(m);
  py_atm(m);
  py_surf(m);
  py_fwd(m);
  py_cia(m);
  py_operators(m);
  py_lbl(m);
  py_sensor(m);
  py_disort(m);
  py_retrieval(m);
  py_lookup(m);
  py_nlte(m);
  py_auto_agenda_operators(m);

  // Must be last, it contains automatic conversion operations
  py_workspace(ws);

  // Extras calling pure internal functions
  py_constants(m);
  py_conversions(m);
  py_physics(m);
  py_global(m);
  py_physics(m);
  py_math(m);
  py_hitran(m);
  py_igrf(m);
  py_zeeman(m);
  py_file(m);

  py::set_leak_warnings(false);
  py::set_implicit_cast_warnings(false);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize module\n{}", e.what()));
}
}  // namespace Python
