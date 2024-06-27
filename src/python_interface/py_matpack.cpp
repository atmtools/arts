#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>

#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"

namespace Python {
template <typename Type, typename T>
void common_math_interface(py::class_<T>& cls_) {
  cls_.def(py::self + Type())
      .def(py::self - Type())
      .def(py::self * Type())
      .def(py::self / Type())
      .def(py::self += Type(), py::rv_policy::none)
      .def(py::self -= Type(), py::rv_policy::none)
      .def(py::self *= Type(), py::rv_policy::none)
      .def(py::self /= Type(), py::rv_policy::none)
      .def(py::self == Type())
      .def(py::self != Type())
      .def(py::self <= Type())
      .def(py::self < Type())
      .def(py::self >= Type())
      .def(py::self > Type())
      .def(Type() + py::self)
      .def(Type() - py::self)
      .def(Type() * py::self)
      .def(Type() / py::self)
      .def(Type() == py::self)
      .def(Type() != py::self)
      .def(Type() <= py::self)
      .def(Type() < py::self)
      .def(Type() >= py::self)
      .def(Type() > py::self);
}

template <typename T>
void common_self_math_interface(py::class_<T>& cls_) {
  cls_.def(+py::self)
      .def(-py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self / py::self)
      .def(py::self += py::self, py::rv_policy::none)
      .def(py::self -= py::self, py::rv_policy::none)
      .def(py::self *= py::self, py::rv_policy::none)
      .def(py::self /= py::self, py::rv_policy::none)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self <= py::self)
      .def(py::self < py::self)
      .def(py::self >= py::self)
      .def(py::self > py::self);
}

void py_matpack(py::module_& m) try {
  py::class_<Vector2> cv2(m, "Vector2");
  py::class_<Vector3> cv3(m, "Vector3");
  matpack_constant_interface(cv2);
  workspace_group_interface(cv2);
  matpack_constant_interface(cv3);
  workspace_group_interface(cv3);

  py::class_<Range>(m, "Range")
      .def(py::init<Index, Index, Index>(),
           "offset"_a,
           "extent"_a,
           "stride"_a = 1,
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

  auto a1 = py::bind_vector<ArrayOfVector, py::rv_policy::reference_internal>(
      m, "ArrayOfVector");
  workspace_group_interface(a1);
  vector_interface(a1);
  auto a2 =
      py::bind_vector<ArrayOfArrayOfVector, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfVector");
  workspace_group_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfMatrix, py::rv_policy::reference_internal>(
      m, "ArrayOfMatrix");
  workspace_group_interface(a3);
  vector_interface(a3);
  auto a4 =
      py::bind_vector<ArrayOfArrayOfMatrix, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfMatrix");
  workspace_group_interface(a4);
  vector_interface(a4);
  auto a5 = py::bind_vector<ArrayOfTensor3, py::rv_policy::reference_internal>(
      m, "ArrayOfTensor3");
  workspace_group_interface(a5);
  vector_interface(a5);
  auto a6 =
      py::bind_vector<ArrayOfArrayOfTensor3, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfTensor3");
  workspace_group_interface(a6);
  vector_interface(a6);
  auto a7 = py::bind_vector<ArrayOfTensor4, py::rv_policy::reference_internal>(
      m, "ArrayOfTensor4");
  workspace_group_interface(a7);
  vector_interface(a7);
  auto a8 = py::bind_vector<ArrayOfTensor5, py::rv_policy::reference_internal>(
      m, "ArrayOfTensor5");
  workspace_group_interface(a8);
  vector_interface(a8);
  auto a9 = py::bind_vector<ArrayOfTensor6, py::rv_policy::reference_internal>(
      m, "ArrayOfTensor6");
  workspace_group_interface(a9);
  vector_interface(a9);
  auto a10 =
      py::bind_vector<ArrayOfArrayOfTensor6, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfTensor6");
  workspace_group_interface(a10);
  vector_interface(a10);
  auto a11 = py::bind_vector<ArrayOfTensor7, py::rv_policy::reference_internal>(
      m, "ArrayOfTensor7");
  workspace_group_interface(a11);
  vector_interface(a11);
  auto a12 = py::bind_vector<ArrayOfVector2, py::rv_policy::reference_internal>(
      m, "ArrayOfVector2");
  workspace_group_interface(a12);
  vector_interface(a12);
  auto a13 = py::bind_vector<ArrayOfVector3, py::rv_policy::reference_internal>(
      m, "ArrayOfVector3");
  workspace_group_interface(a13);
  vector_interface(a13);

  py::class_<Rational> rat(m, "Rational");
  common_math_interface<Index>(rat);
  common_self_math_interface(rat);
  rat.def(py::init<Index, Index>(), "n"_a = 0, "d"_a = 1, "Default rational")
      .def(py::init<String>(), "From :class:`str`")
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
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
  workspace_group_interface(g1);

  py::implicitly_convertible<Vector, AscendingGrid>();
  py::implicitly_convertible<AscendingGrid, Vector>();

  auto b1 =
      py::bind_vector<ArrayOfAscendingGrid, py::rv_policy::reference_internal>(
          m, "ArrayOfAscendingGrid");
  workspace_group_interface(b1);
  vector_interface(b1);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize matpack\n", e.what()));
}
}  // namespace Python
