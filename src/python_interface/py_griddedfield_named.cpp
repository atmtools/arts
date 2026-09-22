#include <debug.h>
#include <sensor_meta_info.h>

#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"
#include "py_griddedfield_helpers.h"

namespace Python {
void py_griddedfield_named(py::module_& m) {
  py::class_<NamedGriddedField2> ngf2(m, "NamedGriddedField2");
  py::class_<NamedGriddedField3> ngf3(m, "NamedGriddedField3");
  gridded_data_interface(ngf2);
  generic_interface(ngf2);
  gridded_data_interface(ngf3);
  generic_interface(ngf3);

  py::class_<ZenGriddedField1> zgf1n(m, "ZenGriddedField1");
  gridded_data_interface(zgf1n);
  generic_interface(zgf1n);

  py::class_<GriddedField1Named> gf1n(m, "GriddedField1Named");
  gridded_data_interface(gf1n);
  generic_interface(gf1n);

  py::class_<ComplexGriddedField2> gf2c(m, "ComplexGriddedField2");
  gridded_data_interface(gf2c);
  generic_interface(gf2c);

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

  py::class_<sensor::CameraGriddedField> cgf(m, "CameraGriddedField");
  cgf.doc() = "A gridded field for camera sensors";
  gridded_data_interface(cgf);
  generic_interface(cgf);
}
}  // namespace Python
