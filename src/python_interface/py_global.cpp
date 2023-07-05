
#include <global_data.h>
#include <parameters.h>
#include <py_auto_interface.h>
#include <pybind11/pybind11.h>
#include <workspace_global_data.h>

#include "isotopologues.h"
#include "py_macros.h"

extern Parameters parameters;
extern String out_basename;

namespace Python {

void py_global(py::module_& m) {
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
      .def_readwrite_static(
          "out_basename",
          &out_basename,
          py::doc(":class:`~pyarts.arts.String` Automatic path for saving data"))
      .doc() = "Access to static settings data";

  global.def(
      "get_md_data",
      []() { return global_data::md_data; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`~pyarts.arts.ArrayOfMdRecord`"
              "\n    List of global method records"));

  global.def(
      "get_MdMap",
      []() { return global_data::MdMap; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Dictionary of global method records' position"));

  global.def(
      "get_md_data_raw",
      []() { return global_data::md_data_raw; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`~pyarts.arts.ArrayOfMdRecord`"
              "\n    List of raw global method records"));

  global.def(
      "get_MdRawMap",
      []() { return global_data::MdRawMap; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Dictionary of global raw method records' position"));

  global.def(
      "get_agenda_data",
      []() { return global_data::agenda_data; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`~pyarts.arts.ArrayOfAgRecord`"
              "\n    List of global agenda records"));

  global.def(
      "get_AgendaMap",
      []() { return global_data::AgendaMap; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Dictionary of global raw agenda records' position"));

  py::class_<GroupRecord>(global, "GroupRecord")
      .def_readwrite(
          "name", &GroupRecord::name, ":class:`~pyarts.arts.String` Group name")
      .def_readwrite("desc",
                     &GroupRecord::desc,
                     ":class:`~pyarts.arts.String` Group description")
      .doc() = "Global workspace data golder - do not use manually";

  py::class_<ArrayOfGroupRecord>(global, "ArrayOfGroupRecord")
      .PythonInterfaceArrayDefault(GroupRecord)
      .doc() = ":class:`list` of :class:`~pyarts.arts.globals.GroupRecord`";

  global.def(
      "get_wsv_groups",
      []() { return global_data::wsv_groups; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`~pyarts.arts.globals.ArrayOfGroupRecord`"
              "\n    List of workspace groups"));

  global.def(
      "get_WsvGroupMap",
      []() { return global_data::WsvGroupMap; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`dict`"
              "\n    Dictionary of workspace groups' position"));

  global.def(
      "get_wsv_data",
      []() { return global_data::wsv_data; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`~pyarts.arts.ArrayOfWsvRecord`"
              "\n    List of workspace variable records"));

  global.def(
      "get_WsvMap",
      []() { return global_data::WsvMap; },
      py::doc("Get a copy of the global data variable\n\n"
              "Return\n------\n:class:`dict`"
              "\n    List of workspace variable records"));

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

}  // void py_global(py::module_& m)
}  // namespace Python
