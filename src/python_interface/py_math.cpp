#include <pybind11/cast.h>
#include <pybind11/pybind11.h>

#include "wigner_functions.h"

namespace Python {
namespace py = pybind11;

void py_math(py::module_& m) {
  auto math = m.def_submodule("math");
  math.doc() = "Contains select mathematics from Arts internal functions";

  math.def(
      "make_wigner_ready",
      &make_wigner_ready,
      py::arg("fastest") = 250,
      py::arg("largest") = 20000000,
      py::arg("symbol_size") = 6,
      py::doc("Initialize a Wigner computation block for wigner3j/wigner6j"));

  math.def("wigner3j",
           &wigner3j,
           py::arg("j1"),
           py::arg("j2"),
           py::arg("j3"),
           py::arg("m1"),
           py::arg("m2"),
           py::arg("m3"),
           py::doc(R"--(Computes the Wigner 3J symbol

         /                \
         |  j1   j2   j3  |
output = |                |
         |  m1   m2   m3  |
         \                /

Note that an appropriately large call to make_wigner_ready(fastest, largest, 3 or 6)
must have been made ahead of time
)--"));

  math.def("wigner6j",
           &wigner6j,
           py::arg("j1"),
           py::arg("j2"),
           py::arg("j3"),
           py::arg("l1"),
           py::arg("l2"),
           py::arg("l3"),
           py::doc(R"--(Computes the Wigner 6J symbol

         /                \
         |  j1   j2   j3  |
output = <                >
         |  l1   l2   l3  |
         \                /

Note that an appropriately large call to make_wigner_ready(fastest, largest, 6)
must have been made ahead of time
)--"));

  math.def(
      "dwigner3j",
      &dwigner3j,
      py::arg("M"),
      py::arg("J1"),
      py::arg("J2"),
      py::arg("J"),
      py::doc(
          R"--(Computes the Wigner 3J symbol using floating point approximation

         /               \
         |  J1   J2   J  |
output = |               |
         |  M   -M    0  |
         \               /

Note that we expect better output from the pure wigner3j function for the same input
)--"));

  math.def(
      "dwigner6j",
      &dwigner6j,
      py::arg("A"),
      py::arg("B"),
      py::arg("C"),
      py::arg("D"),
      py::arg("F"),
      py::doc(R"--(Computes the Wigner 6J symbol

         /             \
         |  A   B   1  |
output = <             >
         |  D   C   F  |
         \             /

Note that we expect better output from the pure wigner3j function for the same input

Also be careful about the order of the symbols bing sent in, they are not entirerly
intuitive
)--"));
}
}  // namespace Python
