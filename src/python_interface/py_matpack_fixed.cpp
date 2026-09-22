#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "python_interface.h"

namespace Python {
void py_matpack_fixed(py::module_& m) {
  py::class_<Vector2>         cv2(m, "Vector2");
  py::class_<Vector3>         cv3(m, "Vector3");
  py::class_<Vector4>         cv4(m, "Vector4");
  py::class_<Vector7>         cv7(m, "Vector7");
  py::class_<Matrix33>        cm33(m, "Matrix33");
  py::class_<Matrix44>        cm44(m, "Matrix44");
  py::class_<ComplexMatrix44> ccm44(m, "ComplexMatrix44");
  cv2.doc()   = "Fixed size vector of shape [2]";
  cv3.doc()   = "Fixed size vector of shape [3]";
  cv4.doc()   = "Fixed size vector of shape [4]";
  cv7.doc()   = "Fixed size vector of shape [7]";
  cm33.doc()  = "Fixed size matrix of shape [3, 3]";
  cm44.doc()  = "Fixed size matrix of shape [4, 4]";
  ccm44.doc() = "Fixed size matrix of shape [4, 4]";
  matpack_constant_interface(cv2);
  matpack_constant_interface(cv3);
  matpack_constant_interface(cv4);
  matpack_constant_interface(cv7);
  matpack_constant_interface(cm33);
  matpack_constant_interface(cm44);
  matpack_constant_interface(ccm44);
  generic_interface(cv2);
  generic_interface(cv3);
  generic_interface(cv4);
  generic_interface(cv7);
  generic_interface(cm33);
  generic_interface(cm44);
  generic_interface(ccm44);

  py::class_<StridedRange>(m, "StridedRange")
      .def(py::init<Index, Index, Index>(), "offset"_a, "extent"_a, "stride"_a = 1, "Valued initialization")
      .doc() = "A strided range, used to select parts of a matpack type";

  py::class_<Range>(m, "Range")
      .def(py::init<Index, Index>(), "offset"_a, "extent"_a, "Valued initialization")
      .def_ro("offset", &Range::offset, "The first element index.\n\n.. :class:`int`")
      .def_ro("extent", &Range::nelem, "The number of elements.\n\n.. :class:`int`")
      .doc() = "A range, used to select parts of a matpack type";
}
}  // namespace Python
