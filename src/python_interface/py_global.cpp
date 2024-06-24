
#include <arts_omp.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>
#include <parameters.h>
#include <python_interface.h>

extern Parameters parameters;

namespace Python {

void py_global(py::module_& m) try {
  auto global  = m.def_submodule("globals");
  global.doc() = "Global settings and data";

  py::class_<Parameters>(global, "parameters")
      .def_rw_static(
          "includepath",
          &parameters.includepath,

          ":class:`~pyarts.arts.ArrayOfString` Automatic include paths")
      .def_rw_static("datapath",
                     &parameters.datapath,
                     ":class:`~pyarts.arts.ArrayOfString` Automatic data paths")
      .def_rw_static(
          "numthreads",
          &parameters.numthreads,

          ":class:`~pyarts.arts.Index` Number of threads allowed to start")
      .doc() = "Access to static settings data";

  py::class_<WorkspaceGroupRecord>(global, "WorkspaceGroupRecord")
      .def_ro("file", &WorkspaceGroupRecord::file)
      .def_ro("desc", &WorkspaceGroupRecord::desc);

  global.def(
      "workspace_groups",
      []() { return internal_workspace_groups(); },
      "Get a copy of all workspace variables\n\n"
      "Return\n------\n:class:`dict`"
      "\n    Map of variables");

  py::class_<WorkspaceVariableRecord>(global, "WorkspaceVariableRecord")
      .def_ro("default_value", &WorkspaceVariableRecord::default_value)
      .def_ro("type", &WorkspaceVariableRecord::type)
      .def_ro("desc", &WorkspaceVariableRecord::desc);

  global.def(
      "workspace_variables",
      []() { return workspace_variables(); },
      "Get a copy of all workspace variables\n\n"
      "Return\n------\n:class:`dict`"
      "\n    Map of variables");

  py::class_<WorkspaceMethodInternalRecord>(global,
                                            "WorkspaceMethodInternalRecord")
      .def_ro("output", &WorkspaceMethodInternalRecord::out)
      .def_ro("input", &WorkspaceMethodInternalRecord::in)
      .def_ro("author", &WorkspaceMethodInternalRecord::author)
      .def_ro("gout", &WorkspaceMethodInternalRecord::gout)
      .def_ro("gout_type", &WorkspaceMethodInternalRecord::gout_type)
      .def_ro("gout_desc", &WorkspaceMethodInternalRecord::gout_desc)
      .def_ro("gin", &WorkspaceMethodInternalRecord::gin)
      .def_ro("gin_type", &WorkspaceMethodInternalRecord::gin_type)
      .def_ro("gin_desc", &WorkspaceMethodInternalRecord::gin_desc)
      .def_ro("gin_value", &WorkspaceMethodInternalRecord::gin_value)
      .def_ro("pass_workspace", &WorkspaceMethodInternalRecord::pass_workspace)
      .def_ro("desc", &WorkspaceMethodInternalRecord::desc);

  global.def(
      "workspace_methods",
      []() { return internal_workspace_methods(); },
      "Get a copy of all workspace methods\n\n"
      "Return\n------\n:class:`dict`"
      "\n    Map of methods");

  py::class_<WorkspaceAgendaInternalRecord>(global,
                                            "WorkspaceAgendaInternalRecord")
      .def_ro("desc", &WorkspaceAgendaInternalRecord::desc)
      .def_ro("output", &WorkspaceAgendaInternalRecord::output)
      .def_ro("input", &WorkspaceAgendaInternalRecord::input)
      .def_ro("array", &WorkspaceAgendaInternalRecord::array)
      .doc() = "Agenda records used as workspace variables";

  global.def(
      "workspace_agendas",
      []() { return internal_workspace_agendas(); },
      "Get a copy of all workspace agendas\n\n"
      "Return\n------\n:class:`dict`"
      "\n    Map of agendas");

  global.def(
      "all_isotopologues",
      []() { return Species::Isotopologues; },
      "List of all valid `~pyarts.arts.SpeciesIsotopeRecord`");

#ifdef _OPENMP
  global.def("omp_get_max_threads",
             &omp_get_max_threads,
             "Get maximum number of OpenMP threads\n\n"
             "Return\n------\n:class:`int`"
             "\n    Number of maximum threads");

  global.def(
      "omp_set_num_threads",
      [startup_max = omp_get_max_threads(), one = 1](int threads) {
        omp_set_num_threads(std::clamp<int>(threads, one, startup_max));
      },
      R"--(
Sets the maximum number of OpenMP threads

The input is clamped between 1 and n, where n
is the startup value of omp_get_max_threads
(which is also the default value)

Parameters
----------
    threads : int
        Number of threads
)--"),
      py::arg("threads") = omp_get_max_threads();
#endif

} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize global\n", e.what()));
}
}  // namespace Python
