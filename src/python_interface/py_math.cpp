#include <fastgl.h>
#include <legendre.h>
#include <matpack.h>
#include <wigner_functions.h>

#include <utility>

#include <nanobind/nanobind.h>

namespace Python {
namespace py = nanobind;
using namespace nanobind::literals;

void py_math(py::module_& m) try {
  auto math = m.def_submodule("math");
  math.doc() = "Contains select mathematics from Arts internal functions";

  math.def(
      "make_wigner_ready",
      &make_wigner_ready,
      py::arg("fastest") = 250,
      py::arg("largest") = 20000000,
      py::arg("symbol_size") = 6,
          R"--(Initialize a Wigner computation block for :func:`wigner3j` or :func:`wigner6j`

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
)--");

  math.def("wigner3j",
           &wigner3j,
           py::arg("j1"),
           py::arg("j2"),
           py::arg("j3"),
           py::arg("m1"),
           py::arg("m2"),
           py::arg("m3"),
           R"--(Computes the Wigner 3J symbol

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
)--");

  math.def("wigner6j",
           &wigner6j,
           py::arg("j1"),
           py::arg("j2"),
           py::arg("j3"),
           py::arg("l1"),
           py::arg("l2"),
           py::arg("l3"),
           R"--(Computes the Wigner 6J symbol

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
)--");

  math.def(
      "dwigner3j",
      &dwigner3j,
      py::arg("M"),
      py::arg("J1"),
      py::arg("J2"),
      py::arg("J"),
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
)--");

  math.def(
      "dwigner6j",
      &dwigner6j,
      py::arg("A"),
      py::arg("B"),
      py::arg("C"),
      py::arg("D"),
      py::arg("F"),
      R"--(Computes the Wigner 6J symbol using floating point approximation

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
)--");

  math.def(
      "leggauss",
      [](Index deg) {
        if (deg < 1) throw std::invalid_argument("Degree must be at least 1");
        auto out = std::make_pair<Vector, Vector>(deg, deg);
        Legendre::GaussLegendre(out.first, out.second);
        return out;
      },
      "deg"_a = Index{0},
      R"(Computes the Gauss-Legendre quadrature

Parameters
----------
deg : int
    The degree of the Gauss-Legendre quadrature

Returns
-------
x : List[float]
    The abscissas
w : List[float]
    The weights
)");

  math.def(
      "pdleggauss",
      [](Index deg) {
        if (deg < 1) throw std::invalid_argument("Degree must be at least 1");
        if (deg % 2) throw std::invalid_argument("Degree must be even");
        auto out = std::make_pair<Vector, Vector>(deg / 2, deg / 2);
        Legendre::PositiveDoubleGaussLegendre(out.first, out.second);
        return out;
      },
      "deg"_a = Index{0},
      R"(Computes the Positive Double Gauss-Legendre quadrature

Parameters
----------
deg : int
    The degree of the Gauss-Legendre quadrature

Returns
-------
x : List[float]
    The abscissas
w : List[float]
    The weights
)");

  math.def(
      "pleggauss",
      [](Index deg) {
        if (deg < 1) throw std::invalid_argument("Degree must be at least 1");
        auto out = std::make_pair<Vector, Vector>(deg, deg);
        Legendre::PositiveGaussLegendre(out.first, out.second);
        return out;
      },
      "deg"_a = Index{0},
      R"(Computes the Positive Gauss-Legendre quadrature

Parameters
----------
deg : int
    The degree of the Gauss-Legendre quadrature

Returns
-------
x : List[float]
    The abscissas
w : List[float]
    The weights
)");

  math.def(
      "schmidt_legendre_polynomial",
      &Legendre::schmidt,
      py::arg("theta"),
      py::arg("nmax"),
      R"(Computes the Positive Gauss-Legendre quadrature

Parameters
----------
theta : Numeric
    The colatitude in radians
nmax : int
    The maximum degree of the polynomial matrix

Returns
-------
Pnm : Matrix
    The Polynominal matrix (nmax+1 x nmax+1; valid the the left of the diagonal only)
dPnm : Matrix
    The Polynominal matrix derivative (nmax+1 x nmax+1; valid the the left of the diagonal only)
)");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize math\n", e.what()));
}
}  // namespace Python
