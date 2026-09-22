#include <nanobind/stl/bind_vector.h>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"
#include "py_griddedfield_helpers.h"

namespace Python {
void py_griddedfield_extra(py::module_& m) {
  py::class_<GeodeticField3> geo3(m, "GeodeticField3");
  gridded_data_interface(geo3);
  generic_interface(geo3);
  implicit_convert_gf<GriddedField3>(geo3);

  py::class_<GeodeticField2> geo2(m, "GeodeticField2");
  gridded_data_interface(geo2);
  generic_interface(geo2);
  implicit_convert_gf<GriddedField2>(geo2);

  py::class_<SortedGriddedField4> sgf4num(m, "SortedGriddedField4");
  gridded_data_interface(sgf4num);
  generic_interface(sgf4num);
  implicit_convert_gf<GriddedField4>(sgf4num);

  py::class_<SortedGriddedField5> sgf5num(m, "SortedGriddedField5");
  gridded_data_interface(sgf5num);
  generic_interface(sgf5num);
  implicit_convert_gf<GriddedField5>(sgf5num);

  py::class_<SortedGriddedField6> sgf6num(m, "SortedGriddedField6");
  gridded_data_interface(sgf6num);
  generic_interface(sgf6num);
  implicit_convert_gf<GriddedField6>(sgf6num);

  py::class_<GriddedSpectralField6> gsf6(m, "GriddedSpectralField6");
  gridded_data_interface(gsf6);
  generic_interface(gsf6);
  implicit_convert_gf<GriddedField6>(gsf6);

  auto a1 = py::bind_vector<ArrayOfGriddedField1, py::rv_policy::reference_internal>(m, "ArrayOfGriddedField1");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfGriddedField2, py::rv_policy::reference_internal>(m, "ArrayOfGriddedField2");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfGriddedField3, py::rv_policy::reference_internal>(m, "ArrayOfGriddedField3");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfGriddedField4, py::rv_policy::reference_internal>(m, "ArrayOfGriddedField4");
  generic_interface(a4);
  vector_interface(a4);
  auto b1 =
      py::bind_vector<ArrayOfArrayOfGriddedField1, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfGriddedField1");
  generic_interface(b1);
  vector_interface(b1);
  auto b2 =
      py::bind_vector<ArrayOfArrayOfGriddedField2, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfGriddedField2");
  generic_interface(b2);
  vector_interface(b2);
  auto b3 =
      py::bind_vector<ArrayOfArrayOfGriddedField3, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfGriddedField3");
  generic_interface(b3);
  vector_interface(b3);
  auto c1 =
      py::bind_vector<ArrayOfGriddedField1Named, py::rv_policy::reference_internal>(m, "ArrayOfGriddedField1Named");
  generic_interface(c1);
  vector_interface(c1);
  auto d2 =
      py::bind_vector<ArrayOfNamedGriddedField2, py::rv_policy::reference_internal>(m, "ArrayOfNamedGriddedField2");
  generic_interface(d2);
  vector_interface(d2);

  auto vsgf1num =
      py::bind_vector<Array<SortedGriddedField1>, py::rv_policy::reference_internal>(m, "ArrayOfSortedGriddedField1");
  vsgf1num.doc() = "A list of :class:`~pyarts3.arts.SortedGriddedField1`";
  generic_interface(vsgf1num);
  vector_interface(vsgf1num);
}
}  // namespace Python
