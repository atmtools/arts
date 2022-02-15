
#include <global_data.h>
#include <parameters.h>
#include <py_auto_interface.h>

#include "py_macros.h"

extern Parameters parameters;

namespace Python {
void py_global(py::module_& m) {
  py::class_<Parameters>(m, "parameters")
      .def_readwrite_static("includepath", &parameters.includepath)
      .def_readwrite_static("datapath", &parameters.datapath)
      .def_readwrite_static("numthreads", &parameters.numthreads)
      .doc() = "Access to static settings data";

  m.def(
      "md_data",
      []() { return global_data::md_data; },
      "Copy of global attribute");

  m.def(
      "MdMap", []() { return global_data::MdMap; }, "Copy of global attribute");

  m.def(
      "md_data_raw",
      []() { return global_data::md_data_raw; },
      "Copy of global attribute");

  m.def(
      "MdRawMap",
      []() { return global_data::MdRawMap; },
      "Copy of global attribute");

  m.def(
      "agenda_data",
      []() { return global_data::agenda_data; },
      "Copy of global attribute");

  m.def(
      "AgendaMap",
      []() { return global_data::AgendaMap; },
      "Copy of global attribute");

  m.def(
      "wsv_group_names",
      []() { return global_data::wsv_group_names; },
      "Copy of global attribute");

  m.def(
      "WsvGroupMap",
      []() { return global_data::WsvGroupMap; },
      "Copy of global attribute");

#ifdef _OPENMP
  m.def("omp_get_max_threads",
        &omp_get_max_threads,
        "Get maximum number of OpenMP threads");

  m.def(
      "omp_set_num_threads",
      [n = omp_get_max_threads()](int threads) {
        if (threads > 0)
          omp_set_num_threads(std::min(threads, n));
        else
          omp_set_num_threads(n);
      },
      R"--(--

Sets the maximum number of OpenMP threads

With n as the startup maximum number of OpenMP threads:
    If threads is positive the value is set to min(threads, n).
    If threads is zero or negative, the value is set to n.

Parameters:
-----------
    threads (int): Number of threads
)--",
      py::arg("threads") = 0);
}
#endif

}  // namespace Python