#include <nanobind/stl/bind_vector.h>
#include <rtepack.h>

#include <stdexcept>

#include "debug.h"
#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"

namespace Python {
void py_griddedfield(py::module_& m) try {
  py::class_<GriddedField1> gf1(m, "GriddedField1");
  py::class_<GriddedField2> gf2(m, "GriddedField2");
  py::class_<GriddedField3> gf3(m, "GriddedField3");
  py::class_<GriddedField4> gf4(m, "GriddedField4");
  py::class_<GriddedField5> gf5(m, "GriddedField5");
  py::class_<GriddedField6> gf6(m, "GriddedField6");
  gridded_data_interface(gf1);
  workspace_group_interface(gf1);
  gridded_data_interface(gf2);
  workspace_group_interface(gf2);
  gridded_data_interface(gf3);
  workspace_group_interface(gf3);
  gridded_data_interface(gf4);
  workspace_group_interface(gf4);
  gridded_data_interface(gf5);
  workspace_group_interface(gf5);
  gridded_data_interface(gf6);
  workspace_group_interface(gf6);

  py::class_<NamedGriddedField2> ngf2(m, "NamedGriddedField2");
  py::class_<NamedGriddedField3> ngf3(m, "NamedGriddedField3");
  gridded_data_interface(ngf2);
  workspace_group_interface(ngf2);
  gridded_data_interface(ngf3);
  workspace_group_interface(ngf3);

  py::class_<GriddedField1Named> gf1n(m, "GriddedField1Named");
  gridded_data_interface(gf1n);
  workspace_group_interface(gf1n);

  py::class_<ComplexGriddedField2> gf2c(m, "ComplexGriddedField2");
  gridded_data_interface(gf2c);
  workspace_group_interface(gf2c);

  py::class_<StokvecGriddedField6> gf6sv(m, "StokvecGriddedField6");
  gridded_data_interface(gf6sv);
  workspace_group_interface(gf6sv);

  auto a1 =
      py::bind_vector<ArrayOfGriddedField1, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField1");
  workspace_group_interface(a1);
  vector_interface(a1);
  auto a2 =
      py::bind_vector<ArrayOfGriddedField2, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField2");
  workspace_group_interface(a2);
  vector_interface(a2);
  auto a3 =
      py::bind_vector<ArrayOfGriddedField3, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField3");
  workspace_group_interface(a3);
  vector_interface(a3);
  auto a4 =
      py::bind_vector<ArrayOfGriddedField4, py::rv_policy::reference_internal>(
          m, "ArrayOfGriddedField4");
  workspace_group_interface(a4);
  vector_interface(a4);
  auto b1 = py::bind_vector<ArrayOfArrayOfGriddedField1,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField1");
  workspace_group_interface(b1);
  vector_interface(b1);
  auto b2 = py::bind_vector<ArrayOfArrayOfGriddedField2,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField2");
  workspace_group_interface(b2);
  vector_interface(b2);
  auto b3 = py::bind_vector<ArrayOfArrayOfGriddedField3,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfGriddedField3");
  workspace_group_interface(b3);
  vector_interface(b3);
  auto c1 = py::bind_vector<ArrayOfGriddedField1Named,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfGriddedField1Named");
  workspace_group_interface(c1);
  vector_interface(c1);
  auto d2 = py::bind_vector<ArrayOfNamedGriddedField2,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfNamedGriddedField2");
  workspace_group_interface(d2);
  vector_interface(d2);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize gridded field\n{}", e.what()));
}
}  // namespace Python
