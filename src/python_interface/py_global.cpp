
#include <py_auto_interface.h>

#include "isotopologues.h"
#include "py_macros.h"

#include <global_data.h>
#include <parameters.h>
#include <pybind11/pybind11.h>
#include <workspace_global_data.h>

extern Parameters parameters;

namespace Python {

void py_global(py::module_& m) {
  py::class_<Parameters>(m, "parameters")
      .def_readwrite_static("includepath", &parameters.includepath)
      .def_readwrite_static("datapath", &parameters.datapath)
      .def_readwrite_static("numthreads", &parameters.numthreads)
      .doc() = "Access to static settings data";

  m.def(
      "get_md_data",
      []() { return global_data::md_data; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_MdMap",
      []() { return global_data::MdMap; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_md_data_raw",
      []() { return global_data::md_data_raw; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_MdRawMap",
      []() { return global_data::MdRawMap; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_agenda_data",
      []() { return global_data::agenda_data; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_AgendaMap",
      []() { return global_data::AgendaMap; },
      py::doc("Get a copy of the global data variable"));

  auto global = m.def_submodule("global");
  py::class_<GroupRecord>(global, "GroupRecord")
      .def_readwrite("name", &GroupRecord::name)
      .def_readwrite("desc", &GroupRecord::desc);
  py::class_<ArrayOfGroupRecord>(global, "ArrayOfGroupRecord")
      .PythonInterfaceArrayDefault(GroupRecord);

  m.def(
      "get_wsv_groups",
      []() { return global_data::wsv_groups; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_WsvGroupMap",
      []() { return global_data::WsvGroupMap; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_wsv_data",
      []() { return global_data::wsv_data; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_WsvMap",
      []() { return global_data::WsvMap; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_isotopologues",
      [] { return Species::Isotopologues; },
      py::doc("Get a list of the global isotopologues"));

#ifdef _OPENMP
  m.def("omp_get_max_threads",
        &omp_get_max_threads,
        py::doc("Get maximum number of OpenMP threads"));

  m.def(
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
    threads (int): Number of threads
)--"),
      py::arg("threads") = omp_get_max_threads());
#endif

}  // void py_global(py::module_& m)
}  // namespace Python
