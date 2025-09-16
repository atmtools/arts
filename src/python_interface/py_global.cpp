
#include <arts_omp.h>
#include <arts_options.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/chrono.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>
#include <parameters.h>

#include "python_interface.h"

struct global_data {
  static constexpr bool is_lgpl =
#ifdef ARTS_LGPL
      ARTS_LGPL ? true : false;
#else
      false;
#endif

  static constexpr bool has_sht = not is_lgpl and
#ifdef _MSC_VER
                                  false;
#else
                                  true;
#endif

  static constexpr bool has_cdisort =
#ifdef ENABLE_CDISORT
      true;
#else
      false;
#endif

  static constexpr bool has_profiling =
#if ARTS_PROFILING
      true;
#else
      false;
#endif

#ifdef ARTS_SOURCE_DIR
  static constexpr std::string_view arts_source_dir = ARTS_SOURCE_DIR;
#else
  static constexpr std::string_view arts_source_dir = "";
#endif
};

extern Parameters parameters;

namespace Python {

void py_global(py::module_& m) try {
  auto global  = m.def_submodule("globals");
  global.doc() = "Global settings and data";

  py::class_<Parameters>(global, "parameters")
      .def_rw_static("includepath",
                     &parameters.includepath,
                     "Automatic include paths\n\n.. :class:`ArrayOfString`")
      .def_rw_static("datapath",
                     &parameters.datapath,
                     "Automatic data paths\n\n.. :class:`ArrayOfString`")
      .doc() = "Access to static settings data";

  py::class_<WorkspaceGroupRecord>(global, "WorkspaceGroupRecord")
      .def_ro(
          "file", &WorkspaceGroupRecord::file, "File path\n\n.. :class:`str`")
      .def_ro(
          "desc", &WorkspaceGroupRecord::desc, "Description\n\n.. :class:`str`")
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
              "Default value\n\n.. :class:`~pyarts3.arts.Wsv`\n\n.. :class:`None`")
      .def_ro("type", &WorkspaceVariableRecord::type, "Type\n\n.. :class:`str`")
      .def_ro("desc",
              &WorkspaceVariableRecord::desc,
              "Description\n\n.. :class:`str`")
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
      .def_ro("output",
              &WorkspaceMethodInternalRecord::out,
              "Outputs\n\n.. :class:`list[str]`")
      .def_ro("input",
              &WorkspaceMethodInternalRecord::in,
              "Inputs\n\n.. :class:`list[str]`")
      .def_ro("author",
              &WorkspaceMethodInternalRecord::author,
              "Authors\n\n.. :class:`list[str]`")
      .def_ro("gout",
              &WorkspaceMethodInternalRecord::gout,
              "Generic output\n\n.. :class:`list[str]`")
      .def_ro("gout_type",
              &WorkspaceMethodInternalRecord::gout_type,
              "Generic output type\n\n.. :class:`list[str]`")
      .def_ro("gout_desc",
              &WorkspaceMethodInternalRecord::gout_desc,
              "Generic output description\n\n.. :class:`list[str]`")
      .def_ro("gin",
              &WorkspaceMethodInternalRecord::gin,
              "Generic input\n\n.. :class:`list[str]`")
      .def_ro("gin_type",
              &WorkspaceMethodInternalRecord::gin_type,
              "Generic input type\n\n.. :class:`list[str]`")
      .def_ro("gin_desc",
              &WorkspaceMethodInternalRecord::gin_desc,
              "Generic input description\n\n.. :class:`list[str]`")
      .def_ro("gin_value",
              &WorkspaceMethodInternalRecord::gin_value,
              "Generic input default value\n\n.. :class:`list[Wsv | None]`")
      .def_ro("pass_workspace",
              &WorkspaceMethodInternalRecord::pass_workspace,
              "Pass workspace\n\n.. :class:`bool`")
      .def_ro("desc",
              &WorkspaceMethodInternalRecord::desc,
              "Description\n\n.. :class:`str`")
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
      .def_ro("desc",
              &WorkspaceAgendaInternalRecord::desc,
              "Description\n\n.. :class:`str`")
      .def_ro("output",
              &WorkspaceAgendaInternalRecord::output,
              "Outputs\n\n.. :class:`list[str]`")
      .def_ro("input",
              &WorkspaceAgendaInternalRecord::input,
              "Inputs\n\n.. :class:`list[str]`")
      .doc() = "Agenda records used as workspace variables";

  global.def(
      "workspace_agendas",
      []() { return internal_workspace_agendas(); },
      "Get a copy of all workspace agendas\n\n"
      "Return\n------\n:class:`dict`"
      "\n    Map of agendas");

  global.def(
      "workspace_agendas_extra",
      []() { return internal_workspace_agenda_names(); },
      "Get a :class:`dict` of overloaded agendas");

  global.def(
      "option_groups",
      []() -> std::vector<std::string> {
        std::vector<std::string> out(internal_options().size());
        for (Size i = 0; i < out.size(); i++) {
          out[i] = internal_options()[i].name;
        }
        return out;
      },
      "Get a copy of all named options\n\n"
      "Return\n------\n:class:`list`");

  global.def(
      "all_isotopologues",
      []() { return Species::Isotopologues; },
      "List of all valid :class:`~pyarts3.arts.SpeciesIsotopeRecord`");

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
      .def_ro_static(
          "is_lgpl",
          &global_data::is_lgpl,
          "Whether the ARTS library is compiled licensed under the LGPL\n\n.. :class:`bool`")
      .def_ro_static(
          "has_sht",
          &global_data::has_sht,
          "Whether the ARTS library is compiled to be able to use the SHT library\n\n.. :class:`bool`")
      .def_ro_static(
          "has_cdisort",
          &global_data::has_cdisort,
          "Whether the ARTS library is compiled with CDisort support\n\n.. :class:`bool`")
      .def_ro_static(
          "has_profiling",
          &global_data::has_profiling,
          "Whether the ARTS library is compiled with time profiling\n\n.. :class:`bool`")
      .def_ro_static(
          "arts_source_dir",
          &global_data::arts_source_dir,
          "The original ARTS source directory, if available\n\n.. :class:`str`")
      .def("__repr__",
           [](const global_data&) {
             return std::format(R"(Global state of ARTS:

is_lgpl: {}
has_sht: {}
)",
                                global_data::is_lgpl,
                                global_data::has_sht);
           })
      .doc() =
      "A set of global data that we might need from ARTS inside pyarts";

  global.def("time_report",
             &arts::get_report,
             "clear"_a = true,
             R"(Get the time report.

The time report is a :class:`dict` with :class:`int` keys representing threads.

Each entry has another :class:`dict` with :class:`str` key representing the short-name of the C++ method that was timed.

As a thread can call a method multiple times, these results are stored as a :class:`list` of a start and an end :class:`~pyarts3.arts.Time`.

.. note::
    This function is only available if ARTS is compiled with profiling enabled.

    Also be aware that the minimum time is in *native* time units, which
    depends on the platform.

Parameters
----------
    clear : bool
        Clear the report after getting it.  Default: True.

Return
------
See above, :class:`dict`
)");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize global\n{}", e.what()));
}
}  // namespace Python
