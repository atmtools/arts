#include <pybind11/cast.h>
#include <pybind11/pybind11.h>

#include "wigner_functions.h"

namespace Python {
namespace py = pybind11;

void py_math(py::module_& m) try {
  auto math = m.def_submodule("math");
  math.doc() = "Contains select mathematics from Arts internal functions";

  math.def(
      "make_wigner_ready",
      &make_wigner_ready,
      py::arg("fastest") = 250,
      py::arg("largest") = 20000000,
      py::arg("symbol_size") = 6,
      py::doc(R"--(Initialize a Wigner computation block for :func:`wigner3j` or :func:`wigner6j`

Parameters
----------
fastest : int, optional
    The size of the fast table (maybe unused), defaults to 250
largest : int, optional
    The largest symbol combination, defaults to 20000000
symbol_size : int, optional
    The symbol for which the largest symbol and fastest symbol is selected (3 or 6), defaults to 6

Returns
-------
actual_size : int
    The actual state as returned by the library
)--"));

  math.def("wigner3j",
           &wigner3j,
           py::arg("j1"),
           py::arg("j2"),
           py::arg("j3"),
           py::arg("m1"),
           py::arg("m2"),
           py::arg("m3"),
           py::doc(R"--(Computes the Wigner 3J symbol

.. math::
    w_3 = \left(\begin{array}{ccc} j_1&j_2&j_3\\m_1&m_2&m_3\end{array}\right)

Note that an appropriately large call to ``make_wigner_ready(fastest, largest, 3 or 6)``
must have been made ahead of time

Parameters
----------
j1 : int
    As above
j2 : int
    As above
j3 : int
    As above
m1 : int
    As above
m2 : int
    As above
m3 : int
    As above

Returns
-------
w3 : float
    The value
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

.. math::
    w_6 = \left\{\begin{array}{ccc} j_1&j_2&j_3\\l_1&l_2&l_3\end{array}\right\}

Note that an appropriately large call to ``make_wigner_ready(fastest, largest, 6)``
must have been made ahead of time

Parameters
----------
j1 : int
    As above
j2 : int
    As above
j3 : int
    As above
l1 : int
    As above
l2 : int
    As above
l3 : int
    As above

Returns
-------
w3 : float
    The value
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

.. math::
    w_3 = \left(\begin{array}{ccc} J_1&J_2&J\\M&-M&0\end{array}\right)

Note that we expect better output from the pure :func:`wigner3j` function for the same input

Parameters
----------
M : int
    As above
J1 : int
    As above
J2 : int
    As above
J : int
    As above

Returns
-------
w3 : float
    The value
)--"));

  math.def(
      "dwigner6j",
      &dwigner6j,
      py::arg("A"),
      py::arg("B"),
      py::arg("C"),
      py::arg("D"),
      py::arg("F"),
      py::doc(R"--(Computes the Wigner 6J symbol using floating point approximation

.. math::
    w_6 = \left\{\begin{array}{ccc} A&B&1\\D&C&F\end{array}\right\}

Note that we expect better output from the pure :func:`wigner6j` function for the same input

Also be careful about the order of the symbols being sent in, they are not entirerly
intuitive

Parameters
----------
A : int
    As above
B : int
    As above
C : int
    As above
D : int
    As above
F : int
    As above

Returns
-------
w6 : float
    The value
)--"));
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize math\n", e.what()));
}
}  // namespace Python
