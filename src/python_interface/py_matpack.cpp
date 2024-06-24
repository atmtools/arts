#include <nanobind/stl/bind_vector.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_matpack.h"

namespace Python {
void py_matpack(py::module_& m) try {
  py::class_<Vector2> cv2(m, "Vector2");
  py::class_<Vector3> cv3(m, "Vector3");
  matpack_constant_interface(cv2);
  matpack_constant_interface(cv3);

  py::class_<Range>(m, "Range")
      .def(py::init<Index, Index, Index>(),
           py::arg("offset"),
           py::arg("extent"),
           py::arg("stride") = 1,
           "Valued initialization")
      .def("__getstate__",
           [](const Range& self) {
             return std::make_tuple(self.offset, self.extent, self.stride);
           })
      .def("__setstate__",
           [](Range* r, const std::tuple<Index, Index, Index>& t) {
             new (r) Range{std::get<0>(t), std::get<1>(t), std::get<2>(t)};
           })
      .doc() = "A range, used to select parts of a matpack type";

  py::class_<IndexVector> iv1(m, "IndexVector");
  matpack_interface(iv1);

  py::class_<Vector> v1(m, "Vector");
  py::class_<Matrix> v2(m, "Matrix");
  py::class_<Tensor3> v3(m, "Tensor3");
  py::class_<Tensor4> v4(m, "Tensor4");
  py::class_<Tensor5> v5(m, "Tensor5");
  py::class_<Tensor6> v6(m, "Tensor6");
  py::class_<Tensor7> v7(m, "Tensor7");

  matpack_interface(v1);
  matpack_interface(v2);
  matpack_interface(v3);
  matpack_interface(v4);
  matpack_interface(v5);
  matpack_interface(v6);
  matpack_interface(v7);
  workspace_group_interface(v1);
  workspace_group_interface(v2);
  workspace_group_interface(v3);
  workspace_group_interface(v4);
  workspace_group_interface(v5);
  workspace_group_interface(v6);
  workspace_group_interface(v7);

  auto a1 = py::bind_vector<ArrayOfVector>(m, "ArrayOfVector");
  workspace_group_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfVector>(m, "ArrayOfArrayOfVector");
  workspace_group_interface(a2);
  auto a3 = py::bind_vector<ArrayOfMatrix>(m, "ArrayOfMatrix");
  workspace_group_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfMatrix>(m, "ArrayOfArrayOfMatrix");
  workspace_group_interface(a4);
  auto a5 = py::bind_vector<ArrayOfTensor3>(m, "ArrayOfTensor3");
  workspace_group_interface(a5);
  auto a6 = py::bind_vector<ArrayOfArrayOfTensor3>(m, "ArrayOfArrayOfTensor3");
  workspace_group_interface(a6);
  auto a7 = py::bind_vector<ArrayOfTensor4>(m, "ArrayOfTensor4");
  workspace_group_interface(a7);
  auto a8 = py::bind_vector<ArrayOfTensor5>(m, "ArrayOfTensor5");
  workspace_group_interface(a8);
  auto a9 = py::bind_vector<ArrayOfTensor6>(m, "ArrayOfTensor6");
  workspace_group_interface(a9);
  auto a10 = py::bind_vector<ArrayOfArrayOfTensor6>(m, "ArrayOfArrayOfTensor6");
  workspace_group_interface(a10);
  auto a11 = py::bind_vector<ArrayOfTensor7>(m, "ArrayOfTensor7");
  workspace_group_interface(a11);
  auto a12 = py::bind_vector<ArrayOfVector2>(m, "ArrayOfVector2");
  workspace_group_interface(a12);
  auto a13 = py::bind_vector<ArrayOfVector3>(m, "ArrayOfVector3");
  workspace_group_interface(a13);

  py::class_<Rational> rat(m, "Rational");
  rat.def(py::init<Index, Index>(),
          py::arg("n") = 0,
          py::arg("d") = 1,
          "Default rational")
      .def(py::init<String>(), "From :class:`str`")
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
      // .PythonInterfaceCommonMath(Index)
      // .PythonInterfaceCommonMathSelf
      .def_rw("n", &Rational::numer, ":class:`int` Numerator")
      .def_rw("d", &Rational::denom, ":class:`int` Denominator")
      .def("__getstate__",
           [](const Rational& self) {
             return std::make_tuple(self.numer, self.denom);
           })
      .def("__setstate__", [](Rational* r, const std::tuple<Index, Index>& t) {
        new (r) Range{std::get<0>(t), std::get<1>(t)};
      });
  workspace_group_interface(rat);
  py::implicitly_convertible<Index, Rational>();

  py::class_<ComplexVector> comv1(m, "ComplexVector");
  py::class_<ComplexMatrix> comv2(m, "ComplexMatrix");
  py::class_<ComplexTensor3> comv3(m, "ComplexTensor3");
  matpack_interface(comv1);
  matpack_interface(comv2);
  matpack_interface(comv3);

  py::class_<AscendingGrid> g1(m, "AscendingGrid");
  matpack_grid_interface(g1);

  py::implicitly_convertible<Vector, AscendingGrid>();
  py::implicitly_convertible<AscendingGrid, Vector>();

  auto b1 = py::bind_vector<ArrayOfAscendingGrid>(m, "ArrayOfAscendingGrid");
  workspace_group_interface(b1);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize matpack\n", e.what()));
}
}  // namespace Python
