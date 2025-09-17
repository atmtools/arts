#include <fastgl.h>
#include <legendre.h>
#include <matpack.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <hpy_arts.h>
#include <wigner_functions.h>

#include <utility>

namespace Python {
namespace py = nanobind;
using namespace nanobind::literals;

void py_math(py::module_& m) try {
  auto math  = m.def_submodule("math");
  math.doc() = "Contains select mathematics from Arts internal functions";

  math.def(
      "make_wigner_ready",
      &make_wigner_ready,
      "fastest"_a     = 250,
      "largest"_a     = 20000000,
      "symbol_size"_a = 6,
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
           "j1"_a,
           "j2"_a,
           "j3"_a,
           "m1"_a,
           "m2"_a,
           "m3"_a,
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
           "j1"_a,
           "j2"_a,
           "j3"_a,
           "l1"_a,
           "l2"_a,
           "l3"_a,
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

  math.def("dwigner3j",
           &dwigner3j,
           "M"_a,
           "J1"_a,
           "J2"_a,
           "J"_a,
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

  math.def("dwigner6j",
           &dwigner6j,
           "A"_a,
           "B"_a,
           "C"_a,
           "D"_a,
           "F"_a,
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
        std::pair<Vector, Vector> out{Vector(deg), Vector(deg)};
        Legendre::GaussLegendre(out.first, out.second);
        return out;
      },
      "deg"_a = Index{0},
      R"(Computes the Gauss-Legendre quadrature
)");

  math.def(
      "pdleggauss",
      [](Index deg) {
        if (deg < 1) throw std::invalid_argument("Degree must be at least 1");
        if (deg % 2) throw std::invalid_argument("Degree must be even");
        std::pair<Vector, Vector> out{Vector(deg / 2), Vector(deg / 2)};
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
x : list[float]
    The abscissas
w : list[float]
    The weights
)");

  math.def(
      "pleggauss",
      [](Index deg) {
        if (deg < 1) throw std::invalid_argument("Degree must be at least 1");
        std::pair<Vector, Vector> out{Vector(deg), Vector(deg)};
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
x : list[float]
    The abscissas
w : list[float]
    The weights
)");

  py::class_<Legendre::SchmidtMatrix> sm(math, "SchmidtMatrix");
  generic_interface(sm);
  sm.def(py::init<Size>(), "N"_a, "Default Schmidt matrix")
      .def(
          "__getitem__",
          [](const Legendre::SchmidtMatrix& self, std::pair<Index, Index> idx) {
            if (idx.first < 0 or idx.second < 0 or
                static_cast<Size>(idx.first) >= self.N or
                idx.second >= idx.first)
              throw std::out_of_range(
                  "Index out of range.  Valid range is [0, N) for first index and [0, first index] for second index");

            return self[idx.first, idx.second];
          },
          "idx"_a,
          R"(Get the value at the given index pair)")
      .def(
          "__setitem__",
          [](Legendre::SchmidtMatrix& self,
             std::pair<Index, Index> idx,
             Numeric value) {
            if (idx.first < 0 or idx.second < 0 or
                static_cast<Size>(idx.first) >= self.N or
                idx.second >= idx.first)
              throw std::out_of_range(
                  "Index out of range.  Valid range is [0, N) for first index and [0, first index] for second index");

            self[idx.first, idx.second] = value;
          },
          "idx"_a,
          "value"_a,
          R"(Get the value at the given index pair)")
      .doc() = R"(Special matrix for the Schmidt polynomials

This is a lower triangular N-by-N-matrix missing the 0, 0 entry.
Can only be accessed at indices
``(i, j)``, where ``0 <= i < N`` and ``0 <= j <= i``.
)";

  math.def("schmidt_legendre_polynomial",
           &Legendre::schmidt,
           "theta"_a,
           "nmax"_a,
           R"(Computes the Positive Gauss-Legendre quadrature

Parameters
----------
theta : Numeric
    The colatitude in radians
nmax : int
    The maximum degree of the polynomial matrix

Returns
-------
Pnm : ~pyarts3.arts.math.SchmidtMatrix
    The Polynominal matrix (nmax+1 x nmax+1).
dPnm : ~pyarts3.arts.math.SchmidtMatrix
    The Polynominal matrix derivative (nmax+1 x nmax+1).
)");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize math\n{}", e.what()));
}
}  // namespace Python
