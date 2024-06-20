#include <python_interface.h>

#include "hpy_matpack.h"
#include "hpy_arts.h"

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
          }).def("__setstate__",
          [](Range*r, const std::tuple<Index, Index, Index>& t) {
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
  xml_interface(v1);
  xml_interface(v2);
  xml_interface(v3);
  xml_interface(v4);
  xml_interface(v5);
  xml_interface(v6);
  xml_interface(v7);

  py::class_<Rational>(m, "Rational")
      .def(py::init<Index, Index>(),
           py::arg("n") = 0,
           py::arg("d") = 1,
           "Default rational")
      .def(py::init<String>(),
           "From :class:`str`")
      .def(py::init<Numeric>(),
           "From :class:`float`")
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
      // .PythonInterfaceCommonMath(Index)
      // .PythonInterfaceCommonMathSelf
      .def_rw("n", &Rational::numer, ":class:`int` Numerator")
      .def_rw("d", &Rational::denom, ":class:`int` Denominator")
      .def("__getstate__", 
          [](const Rational& self) {
            return std::make_tuple(self.numer, self.denom);
          }).def("__setstate__",
          [](Rational*r, const std::tuple<Index, Index>& t) {
            new (r) Range{std::get<0>(t), std::get<1>(t)};
          });
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
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize matpack\n", e.what()));
}
}  // namespace Python
