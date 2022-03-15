
#include <py_auto_interface.h>

#include "py_macros.h"

#include <global_data.h>
#include <parameters.h>

extern Parameters parameters;
extern String out_basename;

namespace Python {

void py_global(py::module_& m) {
  py::class_<Parameters>(m, "parameters")
      .def_readwrite_static("includepath", &parameters.includepath)
      .def_readwrite_static("datapath", &parameters.datapath)
      .def_readwrite_static("numthreads", &parameters.numthreads)
      .def_readwrite_static("out_basename", &out_basename)
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

  m.def(
      "get_wsv_group_names",
      []() { return global_data::wsv_group_names; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_WsvGroupMap",
      []() { return global_data::WsvGroupMap; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_wsv_data",
      []() { return Workspace::wsv_data; },
      py::doc("Get a copy of the global data variable"));

  m.def(
      "get_WsvMap",
      []() { return Workspace::WsvMap; },
      py::doc("Get a copy of the global data variable"));

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
