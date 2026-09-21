#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "python_interface.h"

namespace Python {
template <typename Type, typename T> void common_math_interface(py::class_<T>& cls_) {
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

template <typename T> void common_self_math_interface(py::class_<T>& cls_) {
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

void py_matpack_complex(py::module_& m) {
  py::class_<Rational> rat(m, "Rational");
  common_math_interface<Index>(rat);
  common_self_math_interface(rat);
  rat.def(py::init<Index, Index>(), "n"_a = 0, "d"_a = 1)
      .def(py::init_implicit<const std::string_view>())
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
      .def_rw("n", &Rational::numer, "Numerator\n\n.. :class:`~pyarts3.arts.Index`")
      .def_rw("d", &Rational::denom, "Denominator\n\n.. :class:`~pyarts3.arts.Index`");
  generic_interface(rat);
  py::implicitly_convertible<Index, Rational>();

  py::class_<ComplexVector>  comv1(m, "ComplexVector");
  py::class_<ComplexMatrix>  comv2(m, "ComplexMatrix");
  py::class_<ComplexTensor3> comv3(m, "ComplexTensor3");
  py::class_<ComplexTensor4> comv4(m, "ComplexTensor4");
  py::class_<ComplexTensor5> comv5(m, "ComplexTensor5");
  comv1.doc() = "A complex vector";
  comv2.doc() = "A complex matrix";
  comv3.doc() = "A complex tensor3";
  comv4.doc() = "A complex tensor4";
  comv5.doc() = "A complex tensor5";
  matpack_interface(comv1);
  matpack_interface(comv2);
  matpack_interface(comv3);
  matpack_interface(comv4);
  matpack_interface(comv5);
  generic_interface(comv1);
  generic_interface(comv2);
  generic_interface(comv3);
  generic_interface(comv4);
  generic_interface(comv5);
}
}  // namespace Python
