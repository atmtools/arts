
#include <parameters.h>

#include <python_interface.h>

#include <arts_omp.h>

extern Parameters parameters;

namespace Python {

void py_global(py::module_& m) try {
  auto global = m.def_submodule("globals");

  py::class_<Parameters>(global, "parameters")
      .def_readwrite_static(
          "includepath",
          &parameters.includepath,
          py::doc(":class:`~pyarts.arts.ArrayOfString` Automatic include paths"))
      .def_readwrite_static(
          "datapath",
          &parameters.datapath,
          py::doc(":class:`~pyarts.arts.ArrayOfString` Automatic data paths"))
      .def_readwrite_static(
          "numthreads",
          &parameters.numthreads,
          py::doc(":class:`~pyarts.arts.Index` Number of threads allowed to start"))
      .doc() = "Access to static settings data";

  global.def(
      "workspace_variables",
      []() { return workspace_variables(); },
      py::doc("Get a copy of all workspace variables\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of variables"));

  global.def(
      "workspace_methods",
      []() { return workspace_methods(); },
      py::doc("Get a copy of all workspace methods\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Map of methods"));

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

} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize global\n", e.what()));
}
}  // namespace Python
