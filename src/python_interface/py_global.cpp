
#include <arts_omp.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>
#include <parameters.h>

#include "nanobind/nanobind.h"
#include "python_interface.h"

struct global_data {
  static constexpr bool is_lgpl =
#ifdef ARTS_LGPL
      ARTS_LGPL ? true : false;
#else
      false;
#endif
};

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
      .def_ro("file", &WorkspaceGroupRecord::file, "File path")
      .def_ro("desc", &WorkspaceGroupRecord::desc, "Description")
      .doc() = "Workspace group records";

  global.def(
      "workspace_groups",
      []() { return internal_workspace_groups(); },
      R"(Get a copy of all workspace variables

Return
------
:class:`dict`
    Map of variables)");

  py::class_<WorkspaceVariableRecord>(global, "WorkspaceVariableRecord")
      .def_ro("default_value",
              &WorkspaceVariableRecord::default_value,
              "Default value")
      .def_ro("type", &WorkspaceVariableRecord::type, "Type")
      .def_ro("desc", &WorkspaceVariableRecord::desc, "Description")
      .doc() = "Workspace variable records";

  global.def(
      "workspace_variables",
      []() { return workspace_variables(); },
      R"(Get a copy of all workspace variables

Return
------
:class:`dict`
    Map of variables)");

  py::class_<WorkspaceMethodInternalRecord>(global,
                                            "WorkspaceMethodInternalRecord")
      .def_ro("output", &WorkspaceMethodInternalRecord::out, "Outputs")
      .def_ro("input", &WorkspaceMethodInternalRecord::in, "Inputs")
      .def_ro("author", &WorkspaceMethodInternalRecord::author, "Authors")
      .def_ro("gout", &WorkspaceMethodInternalRecord::gout, "Generic output")
      .def_ro("gout_type",
              &WorkspaceMethodInternalRecord::gout_type,
              "Generic output type")
      .def_ro("gout_desc",
              &WorkspaceMethodInternalRecord::gout_desc,
              "Generic output description")
      .def_ro("gin", &WorkspaceMethodInternalRecord::gin, "Generic input")
      .def_ro("gin_type",
              &WorkspaceMethodInternalRecord::gin_type,
              "Generic input type")
      .def_ro("gin_desc",
              &WorkspaceMethodInternalRecord::gin_desc,
              "Generic input description")
      .def_ro("gin_value",
              &WorkspaceMethodInternalRecord::gin_value,
              "Generic input default value")
      .def_ro("pass_workspace",
              &WorkspaceMethodInternalRecord::pass_workspace,
              "Pass workspace")
      .def_ro("desc", &WorkspaceMethodInternalRecord::desc, "Description")
      .doc() = "Method records used as workspace variables";

  global.def(
      "workspace_methods",
      []() { return internal_workspace_methods(); },
      R"(Get a copy of all workspace methods

Return
------
:class:`dict`
    Map of methods)");

  py::class_<WorkspaceAgendaInternalRecord>(global,
                                            "WorkspaceAgendaInternalRecord")
      .def_ro("desc", &WorkspaceAgendaInternalRecord::desc, "Description")
      .def_ro("output", &WorkspaceAgendaInternalRecord::output, "Outputs")
      .def_ro("input", &WorkspaceAgendaInternalRecord::input, "Inputs")
      .def_ro("array", &WorkspaceAgendaInternalRecord::array, "Is array")
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
      "List of all valid :class:`~pyarts.arts.SpeciesIsotopeRecord`");

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
      "threads"_a = omp_get_max_threads();
#endif

  py::class_<global_data>(global, "data")
      .def(py::init<>())
      .def_ro_static("is_lgpl",
                     &global_data::is_lgpl,
                     "Whether the ARTS library is licensed under the LGPL")
      .def("__repr__", [](const global_data&) {
        return std::format(R"(Global state of ARTS:

is_lgpl: {}
)",
                           global_data::is_lgpl);
      }).doc() = "A set of global data that we might need from ARTS inside pyarts";
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize global\n", e.what()));
}
}  // namespace Python
