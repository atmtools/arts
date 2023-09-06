
#include <arts_omp.h>
#include <parameters.h>
#include <python_interface.h>

extern Parameters parameters;

namespace Python {

void py_global(py::module_& m) try {
  auto global = m.def_submodule("globals");
  global.doc() = "Global settings and data";

  artsclass<Parameters>(global, "parameters")
      .def_readwrite_static(
          "includepath",
          &parameters.includepath,
          py::doc(
              ":class:`~pyarts.arts.ArrayOfString` Automatic include paths"))
      .def_readwrite_static(
          "datapath",
          &parameters.datapath,
          py::doc(":class:`~pyarts.arts.ArrayOfString` Automatic data paths"))
      .def_readwrite_static(
          "numthreads",
          &parameters.numthreads,
          py::doc(
              ":class:`~pyarts.arts.Index` Number of threads allowed to start"))
      .doc() = "Access to static settings data";

  artsclass<WorkspaceGroupRecord>(global, "WorkspaceGroupRecord")
      .def_readonly("file", &WorkspaceGroupRecord::file)
      .def_readonly("desc", &WorkspaceGroupRecord::desc);

  global.def(
      "workspace_groups",
      []() { return internal_workspace_groups(); },
      py::doc("Get a copy of all workspace variables\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of variables"));

  artsclass<WorkspaceVariableRecord>(global, "WorkspaceVariableRecord")
      .def_readonly("default_value", &WorkspaceVariableRecord::default_value)
      .def_readonly("type", &WorkspaceVariableRecord::type)
      .def_readonly("desc", &WorkspaceVariableRecord::desc);

  global.def(
      "workspace_variables",
      []() { return workspace_variables(); },
      py::doc("Get a copy of all workspace variables\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of variables"));

  artsclass<WorkspaceMethodInternalRecord>(global,
                                            "WorkspaceMethodInternalRecord")
      .def_readonly("output", &WorkspaceMethodInternalRecord::out)
      .def_readonly("input", &WorkspaceMethodInternalRecord::in)
      .def_readonly("author", &WorkspaceMethodInternalRecord::author)
      .def_readonly("gout", &WorkspaceMethodInternalRecord::gout)
      .def_readonly("gout_type", &WorkspaceMethodInternalRecord::gout_type)
      .def_readonly("gout_desc", &WorkspaceMethodInternalRecord::gout_desc)
      .def_readonly("gin", &WorkspaceMethodInternalRecord::gin)
      .def_readonly("gin_type", &WorkspaceMethodInternalRecord::gin_type)
      .def_readonly("gin_desc", &WorkspaceMethodInternalRecord::gin_desc)
      .def_readonly("gin_value", &WorkspaceMethodInternalRecord::gin_value)
      .def_readonly("pass_workspace",
                    &WorkspaceMethodInternalRecord::pass_workspace)
      .def_readonly("desc", &WorkspaceMethodInternalRecord::desc);

  global.def(
      "workspace_methods",
      []() { return internal_workspace_methods(); },
      py::doc("Get a copy of all workspace methods\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of methods"));

  artsclass<WorkspaceAgendaInternalRecord>(global,
                                           "WorkspaceAgendaInternalRecord")
      .def_readonly("desc", &WorkspaceAgendaInternalRecord::desc)
      .def_readonly("output", &WorkspaceAgendaInternalRecord::output)
      .def_readonly("input", &WorkspaceAgendaInternalRecord::input)
      .def_readonly("array", &WorkspaceAgendaInternalRecord::array)
      .doc() = "Agenda records used as workspace variables";

  global.def(
      "workspace_agendas",
      []() { return internal_workspace_agendas(); },
      py::doc("Get a copy of all workspace agendas\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of agendas"));

  global.def(
      "get_isotopologues",
      [] { return Species::Isotopologues; },
      py::doc("Get a list of the global isotopologues\n\n"
              "Return\n------\n:class:`list`"
              "\n    List of :class:`pyarts.arts.IsotopeRecord`"));

#ifdef _OPENMP
  global.def("omp_get_max_threads",
             &omp_get_max_threads,
             py::doc("Get maximum number of OpenMP threads\n\n"
                     "Return\n------\n:class:`int`"
                     "\n    Number of maximum threads"));

  global.def(
      "omp_set_num_threads",
      [startup_max = omp_get_max_threads(), one = 1](int threads) {
        omp_set_num_threads(std::clamp<int>(threads, one, startup_max));
      },
      py::doc(R"--(
Sets the maximum number of OpenMP threads

The input is clamped between 1 and n, where n
is the startup value of omp_get_max_threads
(which is also the default value)

Parameters
----------
    threads : int
        Number of threads
)--"),
      py::arg("threads") = omp_get_max_threads());
#endif

} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize global\n", e.what()));
}
}  // namespace Python
