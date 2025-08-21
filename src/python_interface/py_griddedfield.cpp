#include <debug.h>
#include <nanobind/stl/bind_vector.h>
#include <rtepack.h>

#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "matpack_mdspan_helpers_gridded_data_t.h"
#include "python_interface.h"

static_assert(std::is_nothrow_move_constructible_v<GriddedField3>,
              "Did someone change the gridded field implementation?");

namespace Python {
namespace {
template <typename FromType, typename ToType>
void implicit_convert_gf(py::class_<ToType>& cls) {
  cls.def("__init__", [](ToType* x, const FromType& v) {
    ToType t1(v);
    new (x) ToType(std::move(t1));
  });
  py::implicitly_convertible<FromType, ToType>();
}
}  // namespace

void py_griddedfield(py::module_& m) try {
  py::class_<GriddedField1> gf1(m, "GriddedField1");
  py::class_<GriddedField2> gf2(m, "GriddedField2");
  py::class_<GriddedField3> gf3(m, "GriddedField3");
  py::class_<GriddedField4> gf4(m, "GriddedField4");
  py::class_<GriddedField5> gf5(m, "GriddedField5");
  py::class_<GriddedField6> gf6(m, "GriddedField6");
  gridded_data_interface(gf1);
  generic_interface(gf1);
  gridded_data_interface(gf2);
  generic_interface(gf2);
  gridded_data_interface(gf3);
  generic_interface(gf3);
  gridded_data_interface(gf4);
  generic_interface(gf4);
  gridded_data_interface(gf5);
  generic_interface(gf5);
  gridded_data_interface(gf6);
  generic_interface(gf6);
  implicit_convert_gf<SortedGriddedField1>(gf1);
  implicit_convert_gf<SortedGriddedField2>(gf2);
  implicit_convert_gf<SortedGriddedField3>(gf3);
  implicit_convert_gf<SortedGriddedField4>(gf4);
  implicit_convert_gf<SortedGriddedField5>(gf5);
  implicit_convert_gf<SortedGriddedField6>(gf6);

  py::class_<NamedGriddedField2> ngf2(m, "NamedGriddedField2");
  py::class_<NamedGriddedField3> ngf3(m, "NamedGriddedField3");
  gridded_data_interface(ngf2);
  generic_interface(ngf2);
  gridded_data_interface(ngf3);
  generic_interface(ngf3);

  py::class_<GriddedField1Named> gf1n(m, "GriddedField1Named");
  gridded_data_interface(gf1n);
  generic_interface(gf1n);

  py::class_<ComplexGriddedField2> gf2c(m, "ComplexGriddedField2");
  gridded_data_interface(gf2c);
  generic_interface(gf2c);

  py::class_<StokvecSortedGriddedField1> gf1sv(m, "StokvecSortedGriddedField1");
  gridded_data_interface(gf1sv);
  generic_interface(gf1sv);

  py::class_<StokvecSortedGriddedField2> gf2sv(m, "StokvecSortedGriddedField2");
  gridded_data_interface(gf2sv);
  generic_interface(gf2sv);

  py::class_<StokvecSortedGriddedField3> gf3sv(m, "StokvecSortedGriddedField3");
  gridded_data_interface(gf3sv);
  generic_interface(gf3sv);

  py::class_<StokvecSortedGriddedField4> gf4sv(m, "StokvecSortedGriddedField4");
  gridded_data_interface(gf4sv);
  generic_interface(gf4sv);

  py::class_<StokvecSortedGriddedField5> gf5sv(m, "StokvecSortedGriddedField5");
  gridded_data_interface(gf5sv);
  generic_interface(gf5sv);

  py::class_<StokvecSortedGriddedField6> gf6sv(m, "StokvecSortedGriddedField6");
  gridded_data_interface(gf6sv);
  generic_interface(gf6sv);

  py::class_<SortedGriddedField1> sgf1num(m, "SortedGriddedField1");
  gridded_data_interface(sgf1num);
  generic_interface(sgf1num);
  implicit_convert_gf<GriddedField1>(sgf1num);

  py::class_<SortedGriddedField2> sgf2num(m, "SortedGriddedField2");
  gridded_data_interface(sgf2num);
  generic_interface(sgf2num);
  implicit_convert_gf<GriddedField2>(sgf2num);

  py::class_<SortedGriddedField3> sgf3num(m, "SortedGriddedField3");
  gridded_data_interface(sgf3num);
  generic_interface(sgf3num);
  implicit_convert_gf<GriddedField3>(sgf3num);

  py::class_<GeodeticField3> geo3(m, "GeodeticField3");
  gridded_data_interface(geo3);
  generic_interface(geo3);
  implicit_convert_gf<GriddedField3>(geo3);
  gf3.def(
      "make_geodetic",
      [](const GriddedField3& gf) { return matpack::make_geodetic(gf); },
      "Make the field geodetic");

  py::class_<GeodeticField2> geo2(m, "GeodeticField2");
  gridded_data_interface(geo2);
  generic_interface(geo2);
  implicit_convert_gf<GriddedField2>(geo2);
  gf2.def(
      "make_geodetic",
      [](const GriddedField2& gf) { return matpack::make_geodetic(gf); },
      "Make the field geodetic");

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

  auto a1 =
      py::bind_vector<ArrayOfGriddedField1, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField1");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 =
      py::bind_vector<ArrayOfGriddedField2, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField2");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 =
      py::bind_vector<ArrayOfGriddedField3, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField3");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 =
      py::bind_vector<ArrayOfGriddedField4, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField4");
  generic_interface(a4);
  vector_interface(a4);
  auto b1 = py::bind_vector<ArrayOfArrayOfGriddedField1,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField1");
  generic_interface(b1);
  vector_interface(b1);
  auto b2 = py::bind_vector<ArrayOfArrayOfGriddedField2,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField2");
  generic_interface(b2);
  vector_interface(b2);
  auto b3 = py::bind_vector<ArrayOfArrayOfGriddedField3,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField3");
  generic_interface(b3);
  vector_interface(b3);
  auto c1 = py::bind_vector<ArrayOfGriddedField1Named,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfGriddedField1Named");
  generic_interface(c1);
  vector_interface(c1);
  auto d2 = py::bind_vector<ArrayOfNamedGriddedField2,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfNamedGriddedField2");
  generic_interface(d2);
  vector_interface(d2);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize gridded field\n{}", e.what()));
}
}  // namespace Python
