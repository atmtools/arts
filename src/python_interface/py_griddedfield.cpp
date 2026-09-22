#include <debug.h>
#include <nanobind/stl/bind_vector.h>
#include <rtepack.h>
#include <sensor_meta_info.h>

#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"
#include "py_griddedfield_helpers.h"

static_assert(std::is_nothrow_move_constructible_v<GriddedField3>,
              "Did someone change the gridded field implementation?");

namespace Python {
// Defined in py_griddedfield_named.cpp, py_griddedfield_stokvec.cpp and
// py_griddedfield_extra.cpp respectively.  py_griddedfield.cpp used to bind
// all of these types (and many more) in a single function; that made it one
// of the most expensive translation units in the whole build (both in
// compile time and peak RAM), since nanobind's binding machinery is
// instantiated independently for every bound type.  Splitting the bindings
// across a handful of TUs lets the compiler process them in parallel and
// keeps any single TU's memory footprint down, without changing what gets
// bound.
void py_griddedfield_named(py::module_& m);
void py_griddedfield_stokvec(py::module_& m);
void py_griddedfield_extra(py::module_& m);

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
  gf3.def(
      "make_geodetic", [](const GriddedField3& gf) { return matpack::make_geodetic(gf); }, "Make the field geodetic");
  gf2.def(
      "make_geodetic", [](const GriddedField2& gf) { return matpack::make_geodetic(gf); }, "Make the field geodetic");

  py_griddedfield_named(m);
  py_griddedfield_stokvec(m);
  py_griddedfield_extra(m);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize gridded field\n{}", e.what()));
}
}  // namespace Python
